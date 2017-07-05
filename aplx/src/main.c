/* This program will find the minimum value of a n-order rastrigin function.
 * Every core will behave as a cell that carries several chromosomes.
 * Hence, each core/cell will have DEF_N_CHR_PER_CORE Chromosomes. Initially, it will contain DEF_N_CHR_PER_CORE Chromosomes.
 * During simulation, each cell will receive DEF_N_CHR_PER_CORE*(NUM_CORES_USED-1) Chromosomes from its neighbors,
 * so that in total, it always has DEF_N_CHR_PER_CORE*NUM_CORES_USED Chromosomes.
 * After one generation, it will send DEF_N_CHR_PER_CORE Chromosomes to all neighbors.
 * See myAlgorithm.txt for detail.
 *
 *
 * For visualization:
 * I created a simple matlab in ~/Projects/Matlab/try_rastrigin.m
 * */

#include <spin1_api.h>
#include <stdfix.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>	// needed by abs()

#define REAL            accum
#define REAL_CONST(x)   x##k
#define RAND_MAX_UINT	0xFFFFFFFF
#define RAND_MAX_UINT16	0xFFFF

#define PI		(REAL_CONST(3.14159265359))

#define USE_ABS_OBJECTIVE			TRUE

#define DEF_N_CHR_PER_CORE			2		// Number of chromosomes per core
#define DEF_N_GENES_PER_CHR			2		// the dimension of the function (1..xxx) = number of GA parameters encoded in a chromosome
#define NUM_CORES_USED				16		// assuming that all 16 cores are available

// WARNING: Bahaya tanpa () untuk konstanta yang didefinisikan dengan #define:
#define TOTAL_GENES					(DEF_N_GENES_PER_CHR * DEF_N_CHR_PER_CORE * NUM_CORES_USED)
#define TOTAL_CHROMOSOMES			(DEF_N_CHR_PER_CORE * NUM_CORES_USED)
//---------------------------------------------------------------------------------------------
#define MAX_ITER					1
#define MIN_PARAM					-10.0
#define MAX_PARAM					10.0
#define rhoC						0.25	// crossover rate 25%
#define rhoM						0.1		// mutation rate
#define MUTATION_BOUND				16		// since the range of parameter value is limited (between -5.xx to +5.xx), it doesn't make any sense to use the whole bit of stdfix
#define STOPPING_OBJ_VAL			0.000001

uint tik = 0;
uint myCoreID;
uint myCellID;
uint genCntr = 0;
uint prevCellPacketCntr = 0;
uint iter = 0;
uint eliteParent[2] = {0};					// elitism, just an index to certain chromosome in the current population

#define TIMER_TICK_PERIOD		1000000
#define BROADCAST_PRIORITY		1
#define GA_COMPUTE_PRIORITY		2
#define TIMER_PRIORITY			3

/* GA parameters */
uint Chromosomes[DEF_N_CHR_PER_CORE * NUM_CORES_USED][DEF_N_GENES_PER_CHR] = {0};
uint collectedGenes;					// a counter to count how many chromosomes are collected so far
										// if Chromosomes buffer is full, then start the calculation
REAL objVal[TOTAL_CHROMOSOMES] = {0};
REAL fitVal[TOTAL_CHROMOSOMES] = {0};
REAL probVal[TOTAL_CHROMOSOMES] = {0};
REAL cdfVal[TOTAL_CHROMOSOMES] = {0};
REAL bestObjVal[MAX_ITER] = {0};
uint bestChromosomes[MAX_ITER] = {0};

REAL randValSelection[TOTAL_CHROMOSOMES];
uint randValMating[TOTAL_CHROMOSOMES];	// split payload of MCPL into two part: [1st-partner, 2nd-partner]
uint randValMutPoint[TOTAL_CHROMOSOMES];
//future mate-pair for crossover, and future mutation point

/* Use mersenne-twister random generator */
extern void init_genrand(unsigned long s);
extern uint genrand_int32(void);
extern REAL genrand_fixp(REAL minVal, REAL maxVal, uint seed);


// forward declaration
REAL decode(uint gen);
uint encode(REAL gen);
void showMyChromosomes();
void checkMyTurn(uint cellID);
void bcastMyChromosomes(uint arg0, uint arg1);
void printObjVal();
void objEval(uint arg0, uint arg1);
void fitProbEval(uint arg0, uint arg1);
void doSelection(uint arg0, uint arg1);
void doCrossover(uint arg0, uint arg1);
void doMutation(uint arg0, uint arg1);
void bcastRandVal(uint arg0, uint arg1);

// helper functions
void getBinary(uint num, char *buf);
void shuffle(uint *array, uint n);
REAL roundr(REAL inVal);

/*-------------------------------- IMPLEMENTATION ----------------------------------*/

void hTimer(uint arg0, uint arg1)
{
	//io_printf(IO_BUF, "%u\n", genrand_int32());
	tik++;
	/*
	if(tik==10){
		io_printf(IO_STD, "Has collected %d genes so far!\n", collectedGenes); spin1_delay_us(10000);
		spin1_exit(0);
	}
	*/
}

/* Routing strategy for sending a chromosome (broadcast):
 * key     : 0xBCiissxx, where "BC" is broadcast id for masking, ii is the chromosome ID, ss is the source (cell-ID, not core-ID), and xx is the gen-ID
 * payload : particular gen that will be exchanged
 * mask    : 0xFFFFFFFF
 * dest    : 0xFFFFFFC0 -> intra chip only
 * */
void initRouter()
{
	uint e, i, key, mask, route;
	// broadcast MC, but don't send to its own core
	e = rtr_alloc(NUM_CORES_USED);
	if (e == 0)
		rt_error(RTE_ABORT);
	for(i=0; i<NUM_CORES_USED; i++) {
		key = 0xBC000000;
		key |= (i << 8);				// use CellID because easier for indexing for the Chromosomes buffer
		route = 1 << (i + 7);
		route = (~route) & 0x7FFF80;	// send to all cores except my self
		mask = 0xFF00FF00;
		rtr_mc_set(e+i, key, mask, route);
	}
}

/* hMCPL() will collect chromosomes sent by another cells.
 *
 * */
void hMCPL(uint key, uint payload)
{
	// Debugging:
	//genCntr++;
	//io_printf(IO_BUF, "Getting-%2d -> 0x%x:0x%x\n", genCntr, key, payload);

	uint ii, ss, xx;	// key: 0xBCiissxx, "BC" = broadcast id, ii = chromosome ID, ss = cell-ID, xx = gen-ID
	if((key >> 24) == 0xBC) {
		xx = key & 0x000000FF;
		ss = (key & 0x0000FF00) >> 8;	// in this case, used as detector to check my turn to broadcast
		ii = (key & 0x00FF0000) >> 16;
		Chromosomes[ii][xx] = payload;
		collectedGenes++;
		checkMyTurn(ss);
	}

	// Debugging: check if all chromosomes have been collected
	if(collectedGenes == TOTAL_GENES) {
		iter++;

		// if all genes are collected, leadApp will broadcast: random number for selection, future mate-pair for crossover, and future mutation point
		if(leadAp) {
			io_printf(IO_STD, "Cell-%d has received all chromosomes!\n", myCellID); spin1_delay_us(myCellID*1000);
			//showMyChromosomes();
			spin1_schedule_callback(bcastRandVal, 0, 0, GA_COMPUTE_PRIORITY);
		}

		spin1_schedule_callback(objEval, 0, 0, GA_COMPUTE_PRIORITY);
	}
}

/* leadApp will broadcast: random number for selection, future mate-pair for crossover, and future mutation point
 *
 * */
void bcastRandVal(uint arg0, uint arg1)
{

}

/* We have problem with instaneous broadcast: dropped & missing packets.
 * Alternative solution: sequential broadcast
 * checkMyTurn() will check if it is its turn to broadcast by counting how many packets of previous cell
 * have been collected;
 * */
void checkMyTurn(uint cellID)
{
	//Debugging:
	//io_printf(IO_BUF, "Receiving genes from cell-%d\n", cellID);
	if(myCellID-cellID==1) {
		prevCellPacketCntr++;
		if(prevCellPacketCntr==DEF_N_CHR_PER_CORE*DEF_N_GENES_PER_CHR) {
			prevCellPacketCntr = 0;
			spin1_schedule_callback(bcastMyChromosomes, 0, 0, BROADCAST_PRIORITY);
		}
	}
}

/* initMyChromosomes() will generate initial N-chromosomes that will be stored in the Chromosomes buffer
 * (where N is DEF_N_GENES_PER_CHR).
 * At this point, the collectedChromosomes will be equal to DEF_N_GENES_PER_CHR
 * Afterwards, each cell will wait for chromosomes from another cells.
 * */
void initMyChromosomes()
{
	/* generate random value between MIN_PARAM and MAX_PARAM */
	uint h,i, gg;
	REAL g;

	// the following pointers are necessary to avoid "error: cannot convert to a pointer type"
	void *dst, *src;
	src = &g;
	dst = &gg;

	// initialize own chromosomes, don't care with the others...
	for(h=0; h<TOTAL_CHROMOSOMES; h++)
		for(i=0; i<DEF_N_GENES_PER_CHR; i++) {
			g = genrand_fixp(MIN_PARAM, MAX_PARAM, sv->clock_ms);

			spin1_memcpy(dst, src, sizeof(uint));
			// see if we can encode back the real value from integer:
			// spin1_memcpy(idst, dst, sizeof(uint));

			/*the following doesn't work:
			  spin1_memcpy((void *)gg, (void *)g, sizeof(uint));
			  spin1_memcpy((void *)lagi, (void *)gg, sizeof(uint));*/

			Chromosomes[h][i] = gg;
			//io_printf(IO_BUF, "Generating gen = %k -> 0x%x -> %k\n", g, gg, ig);
		}
	//collectedGenes = DEF_N_CHR_PER_CORE * DEF_N_GENES_PER_CHR;
}

void showMyChromosomes()
{
	uint i,j;
	char binBuf[33];
	binBuf[32] = '\0';
	for(i=0; i<DEF_N_CHR_PER_CORE*NUM_CORES_USED; i++) {
		io_printf(IO_BUF, "Chr-%2d = [ ", i);
		for(j=0; j<DEF_N_GENES_PER_CHR; j++)
			io_printf(IO_BUF, "%k ", Chromosomes[i][j]);
		//io_printf(IO_BUF, "] = 0x");
		io_printf(IO_BUF, "] = ");
		for(j=0; j<DEF_N_GENES_PER_CHR; j++)
			//io_printf(IO_BUF, "0x%x", Chromosomes[i][j]);
			io_printf(IO_BUF, "0x%x ", Chromosomes[i][j]);
		io_printf(IO_BUF, " = ");
		//io_printf(IO_BUF, "] = ");
		for(j=0; j<DEF_N_GENES_PER_CHR; j++) {
			getBinary(Chromosomes[i][j], binBuf);
			io_printf(IO_BUF, "%s", binBuf);
		}
		io_printf(IO_BUF, "\n");
		spin1_delay_us(1000);
	}
}

void getBinary(uint num, char *buf)
{
	uint i;
	for(i=0; i<32; i++) {
		if(num & 1 == 1)
			buf[31-i] = '1';
		else
			buf[31-i] = '0';
		num = num >> 1;
	}
}

void bcastMyChromosomes(uint arg0, uint arg1)
{
	// Debugging
	//io_printf(IO_STD, "Cell-%d broadcasts chromosomes...\n", myCellID);
	//showMyChromosomes();
	uint i, j, key, data, idx;
	uint ii, ss, xx;	// key: 0xBCiissxx, "BC" = broadcast id, ii = chromosome ID, ss = cell-ID, xx = gen-ID
	for(i=0; i<DEF_N_CHR_PER_CORE; i++)
		for(j=0; j<DEF_N_GENES_PER_CHR; j++){
			idx = myCellID*DEF_N_CHR_PER_CORE + i;
			xx = j;
			ss = myCellID << 8;
			ii = idx << 16; // chromosome ID is computed by the sender cell
			key = 0xBC000000;
			key |= (ii | ss | xx);
			data = Chromosomes[idx][j];
			spin1_send_mc_packet(key, data, WITH_PAYLOAD);
			//io_printf(IO_BUF, "Sending 0x%x:0x%x\n", key, data);
			//io_printf(IO_STD, "Broadcast my chromosomes:\n");
			//TODO: masih ada packet drop!!!
			//spin1_delay_us(100);
		}
}

REAL roundr(REAL inVal)
{
	uint base = (uint)inVal;
	uint upper = base + 1;
	REAL conver = inVal+REAL_CONST(0.5);
	if((uint)conver == base)
		return (REAL)base;
	else
		return (REAL)upper;
}

/*----------------------------------- GA Stuffs ----------------------------------------*/
// decode() convert from uint to REAL
REAL decode(uint gen)
{
	REAL result;
	void *dst, *src;
	src = &gen;
	dst = &result;
	spin1_memcpy(dst, src, sizeof(uint));
	return result;
}

// encode() convert from REAL to uint
uint encode(REAL gen)
{
	uint result;
	void *dst, *src;
	src = &gen;
	dst = &result;
	spin1_memcpy(dst, src, sizeof(uint));
	return result;
}

void printObjVal()
{
	uint i, j;
	REAL total = 0.0;
	REAL g;
	REAL c0Val = 10.0 * DEF_N_GENES_PER_CHR;
	REAL ciVal = 10.0;
	REAL objVal[TOTAL_CHROMOSOMES];
	for(i=0; i<TOTAL_CHROMOSOMES; i++) {
		objVal[i] = c0Val;
		for(j=0; j<DEF_N_GENES_PER_CHR; j++) {
			g = decode(Chromosomes[i][j]);
			objVal[i] += (g*g - ciVal*cosf(2*PI*g));
			//io_printf(IO_BUF, "g = %k, objVal[%d] = %k\n", g,i,objVal[i]);
		}
		if(USE_ABS_OBJECTIVE)
			//objVal[i] = abs(objVal[i]);
			if(objVal[i] < 0.0)
				objVal[i] *= -1.0;
		io_printf(IO_BUF, "objVal[%d] = %k\n", i, objVal[i]);
	}

}

/* TODO: objEval() should be a *virtual* function that modifiable for another GA experiment!! */
void objEval(uint arg0, uint arg1)
{
	uint i, j;


	REAL total = 0.0;
	REAL g;
	REAL c0Val = 10.0 * DEF_N_GENES_PER_CHR;
	REAL ciVal = 10.0;
	for(i=0; i<TOTAL_CHROMOSOMES; i++) {
		objVal[i] = c0Val;
		for(j=0; j<DEF_N_GENES_PER_CHR; j++) {
			g = decode(Chromosomes[i][j]);
			objVal[i] += (g*g - ciVal*cosf(2*PI*g));
			//io_printf(IO_BUF, "g = %k, objVal[%d] = %k\n", g,i,objVal[i]);
		}
		if(USE_ABS_OBJECTIVE)
			//objVal[i] = abs(objVal[i]);
			if(objVal[i] < 0.0)
				objVal[i] *= -1.0;
		io_printf(IO_BUF, "objVal[%d] = %k\n", i, objVal[i]);
	}
	// find the best objVal and record it
	REAL bestVal[2];
	uint cIdx[2] = {0};
	// let's put the best one in the first order (index-0)
	if(objVal[0] < objVal[1]) {
		bestVal[0] = objVal[0]; bestVal[1] = objVal[1];
		cIdx[0] = 0; cIdx[1] = 1;
	}
	else {
		bestVal[0] = objVal[1]; bestVal[1] = objVal[0];
		cIdx[0] = 1; cIdx[1] = 0;
	}

	for(i=2; i<TOTAL_CHROMOSOMES; i++) {
		if(objVal[i] < bestVal[0]) {
			// shift the record index-0 to index-1
			cIdx[1] = cIdx[0];
			bestVal[1] = bestVal[0];
			// save the record to index-0
			cIdx[0] = i;
			bestVal[0] = objVal[i];
		}
		else if(objVal[i] < bestVal[1]) {
			cIdx[1] = i;
			bestVal[1] = objVal[i];
		}
	}
	io_printf(IO_BUF, "bestVal = %k on chromosome-%d\n", bestVal[0],cIdx[0]);
	bestObjVal[iter] = bestVal[0];
	bestChromosomes[iter] = cIdx[0];
	// finally, put the best two chromosomes as the elite chromosomes
	eliteParent[0] = cIdx[0];
	eliteParent[1] = cIdx[1];
	//io_printf(IO_BUF, "\nElite = [%d, %d]\n\n", eliteParent[0], eliteParent[1]);
}

void fitProbEval(uint arg0, uint arg1)
{
	uint i, j;
	REAL TFitness = 0.0;
	for(i=0; i<TOTAL_CHROMOSOMES; i++) {
		fitVal[i] = 1.0 / objVal[i];
		TFitness += fitVal[i];
		//io_printf(IO_BUF, "fitVal[%d] = %k, TFitness = %k\n\n", i, fitVal[i], TFitness);
	}
	//io_printf(IO_BUF, "\n final TFitness = %k\n\n", TFitness);
	for(i=0; i<TOTAL_CHROMOSOMES; i++) {
		probVal[i] = fitVal[i] / TFitness;
		//io_printf(IO_BUF, "probVal[%d] = %k\n", i, probVal[i]);
	}
	for(i=0; i<TOTAL_CHROMOSOMES; i++) {
		cdfVal[i] = 0;
		for(j=0; j<(i+1); j++) {
			cdfVal[i] += probVal[j];
		}
		//io_printf(IO_BUF, "cdfVal[%d] = %k\n", i, cdfVal[i]);
	}
}

void doSelection(uint arg0, uint arg1)
{
	uint i, j, k;
	REAL R[TOTAL_CHROMOSOMES];
	// 1st: generate n random number within [0.0, 1.0]
	for(i=0; i<TOTAL_CHROMOSOMES; i++)
		R[i] = genrand_fixp(0.0, 1.0, sv->clock_ms);

	/*
	// Debugging:
	io_printf(IO_BUF, "During Selection...\n---------------------\n");
	for(i=0; i<TOTAL_CHROMOSOMES; i++) {
		io_printf(IO_BUF, "R[%d] = %k, cdf[%d] = %k\n", i, R[i], i, cdfVal[i]);
	}
	io_printf(IO_BUF, "Check: if CDF[j] > R[i], then select the chromosomes\n");
	*/

	// 2nd: if CDF > R, then select the chromosomes
	for(i=0; i<TOTAL_CHROMOSOMES; i++)
		for(j=0; j<TOTAL_CHROMOSOMES; j++) {
			if(cdfVal[j] > R[i]) {
				for(k=0; k<DEF_N_GENES_PER_CHR; k++)
					Chromosomes[i][k] = Chromosomes[j][k];
				break; // break the j-loop
			}
		}
}

void doCrossover(uint arg0, uint arg1)
{
	// 1st: select, which chromosome will be a parent
	uint k, l;
	REAL R;
	uint parent[TOTAL_CHROMOSOMES] = {0};
	// implement elitism:
	for(k=0; k<2; k++)
		parent[k] = eliteParent[k];
	uint Tparent = 2;	// total selected individuals as parents
	for(k=0; k<TOTAL_CHROMOSOMES-2; k++) {
		R = genrand_fixp(0.0, 1.0, sv->clock_ms);
		if(R < rhoC && k!=eliteParent[0] && k!=eliteParent[1]) {
			parent[Tparent] = k;
			Tparent++;
		}
	}

	// Just a trick: if random number is too small, then select all as parents
	// due to elitism, this will never happens
	/*
	if(Tparent < 2) {
		Tparent = TOTAL_CHROMOSOMES;
		for(k=0; k<TOTAL_CHROMOSOMES; k++) parent[k] = k;
	}
	*/

	// if the number is odd, remove one
	// if(Tparent % 2 !=0) Tparent--;	// TODO: this is oversimplified, normally, we remove random parent, not only the last one!!!
	if(Tparent % 2 != 0) {
		uint excluded = genrand_int32() % Tparent;	// generate random index number between 0..Tparent to select a single parent
		for(k=0; k<Tparent-1; k++) {				// check for index from 0..(Tparent-1), at the last index might be removed
			if(k>=excluded)
				parent[k] = parent[k+1];
		}
		Tparent--;
	}

	// 2nd: create a mating list with random pair
	uint group1[TOTAL_CHROMOSOMES/2];
	uint group2[TOTAL_CHROMOSOMES/2];

	shuffle(parent, Tparent);
	for(k=0; k<Tparent/2; k++) {
		group1[k] = parent[k];
		group2[k] = parent[Tparent-k-1];
	}

	// 3rd: generate crossover point that will be shared by group1 and group2: hence it is half of Tparent
	uint cp[TOTAL_CHROMOSOMES/2];
	for(k=0; k<Tparent/2; k++)
		cp[k] = 1 + genrand_int32() % (32 * DEF_N_GENES_PER_CHR - 1);


	// Debugging: print parent groups and the crossover point
	/*
	io_printf(IO_BUF, "group1 = [ ");
	for(k=0; k<Tparent/2; k++)
		io_printf(IO_BUF, "%d ", group1[k]);
	io_printf(IO_BUF, "]\n");
	io_printf(IO_BUF, "group2 = [ ");
	for(k=0; k<Tparent/2; k++)
		io_printf(IO_BUF, "%d ", group2[k]);
	io_printf(IO_BUF, "]\n");
	io_printf(IO_BUF, "cp = [ ");
	for(k=0; k<Tparent/2; k++)
		io_printf(IO_BUF, "%d ", cp[k]);
	io_printf(IO_BUF, "]\n");
	*/

	// 4th: cross the parent
	uint children[2][DEF_N_GENES_PER_CHR];
	uint cGenes[2];		// crossed Genes
	uint tGenes[2];		// temporary crossed Genes
	// genPos determine, in which gene-index in a chromosome does the cp exist
	uint genPos, bitPos, mask;
	char binBuf[33];
	binBuf[32] = '\0';

	for(k=0; k<Tparent/2; k++) {
		genPos = (DEF_N_GENES_PER_CHR - 1) - cp[k]/32;
		//io_printf(IO_BUF, "cp[%d] = %d, genPos = %d\n", k, cp[k], genPos);
		for(l=0; l<DEF_N_GENES_PER_CHR; l++) {
			if(l < genPos) {
				children[0][l] = Chromosomes[group1[k]][l];
				children[1][l] = Chromosomes[group2[k]][l];
			}
			// if the cross point is in the current gene-l, then do the logical swap on that particular gene
			else if(l == genPos) {
				tGenes[0] = Chromosomes[group1[k]][l];
				tGenes[1] = Chromosomes[group2[k]][l];
				// bitPos = cp[k] - 32*(DEF_N_GENES_PER_CHR - 1 - l) - (DEF_N_GENES_PER_CHR - 1 - l);
				bitPos = cp[k] - 32*(DEF_N_GENES_PER_CHR - 1 - l);
				/*
				io_printf(IO_BUF, "cp[%d] = %d, genPos = %d, bitPos = %d\n", k, cp[k], genPos, bitPos);
				getBinary(tGenes[0], binBuf);
				io_printf(IO_BUF, "\nBefore crossing:\ncGenes[0] = %s\n", binBuf);
				getBinary(tGenes[1], binBuf);
				io_printf(IO_BUF, "cGenes[1] = %s\n", binBuf);
				*/
				// apply masking and swap
				mask = 0xFFFFFFFF << bitPos;
				cGenes[0] = (tGenes[0] & mask) | (tGenes[1] & (~mask));
				cGenes[1] = (tGenes[1] & mask) | (tGenes[0] & (~mask));
				/*
				getBinary(cGenes[0], binBuf);
				io_printf(IO_BUF, "After crossing:\ncGenes[0] = %s\n", binBuf);
				getBinary(cGenes[1], binBuf);
				io_printf(IO_BUF, "cGenes[1] = %s\n\n", binBuf);
				*/
				// move to children genes
				children[0][l] = cGenes[0];
				children[1][l] = cGenes[1];
			}
			else {	// otherwise, swap the entire gene
				children[0][l] = Chromosomes[group2[k]][l];
				children[1][l] = Chromosomes[group1[k]][l];
			}
		}
		// then replace the parent with the child
		for(l=0; l<DEF_N_GENES_PER_CHR; l++) {
			Chromosomes[group1[k]][l] = children[0][l];
			Chromosomes[group2[k]][l] = children[1][l];
		}
	}

}

void shuffle(uint *array, uint n)
{
	uint i;
	for (i = 0; i < n - 1; i++)
	{
		uint j = i + genrand_int32() / (RAND_MAX_UINT / (n - i) + 1);
		uint t = array[j];
		array[j] = array[i];
		array[i] = t;
	}
}

void doMutation(uint arg0, uint arg1)
{
	uint TBasePair = TOTAL_GENES * 32;	// total base pair in *population*, NOTE: we use stdfix (32 bit representation)
	uint nMutatedGen = (uint)roundr(rhoM * (REAL)TBasePair);
	uint i;
	uint rc, rg, rb;	// random chromosome-, random gen-, and random base-pair-position
	uint newGen;
	// first, generate mutation point
	//io_printf(IO_BUF, "nMutatedGen = %d\n", nMutatedGen);
	for(i=0; i<nMutatedGen; i++){
		rc = genrand_int32() % TOTAL_CHROMOSOMES;	// which chromosome/individual?
		rg = genrand_int32() % DEF_N_CHR_PER_CORE;	// which gene?
		//rb = genrand_int32() % 32;					// which base pair?
		rb = genrand_int32() % MUTATION_BOUND;					// change only the fraction part?
		newGen = Chromosomes[rc][rg];
		newGen ^= (1 << rb);
		Chromosomes[rc][rg] = newGen;
		//io_printf(IO_BUF, "rc[%d] = %d, rg[%d] = %d, rb[%d] = %d\n", i, rc, i, rg, i, rb);
	}
	// then apply to the current population

}

void c_main(void)
{
	myCoreID = spin1_get_core_id();
	myCellID = myCoreID - 1;

	/* Initialize system */
	unsigned long seed = myCoreID * sv->clock_ms;
	init_genrand(seed);

	/* Generate initial population */
	initMyChromosomes();
	//io_printf(IO_BUF, "After initial generation:\n--------------------------\n");

	do {
		io_printf(IO_BUF, "\n\n======================================= Iteration-%d ===========================================\n", iter);
		io_printf(IO_BUF, "\nCurrent chromosomes:\n");
		showMyChromosomes();
		io_printf(IO_BUF, "\nObjective values of current chromosomes:\n");
		objEval(0,0);		// evaluate current chromosomes and record its best objective value
		fitProbEval(0,0);
		doSelection(0,0);
		io_printf(IO_BUF, "\nCurrent chromosomes after Selection:\n--------------------------\n");
		showMyChromosomes();
		io_printf(IO_BUF, "\nObjective values after Selection:\n");
		printObjVal();		// evaluate current chromosomes and record its best objective value
		doCrossover(0,0);
		io_printf(IO_BUF, "\nCurrent chromosomes after Crossover:\n--------------------------\n");
		showMyChromosomes();
		io_printf(IO_BUF, "\nObjective values after Crossover:\n");
		printObjVal();		// evaluate current chromosomes and record its best objective value
		doMutation(0,0);
		io_printf(IO_BUF, "\nCurrent chromosomes after Mutation:\n--------------------------\n");
		showMyChromosomes();
		io_printf(IO_BUF, "\nObjective values after Mutation:\n");
		printObjVal();		// evaluate current chromosomes and record its best objective value
		io_printf(IO_STD, "The %d-generation with bestObjVal = %k on chromosome-%d\n", iter, bestObjVal[iter], bestChromosomes[iter]);
		iter++;
	} while(iter < MAX_ITER || bestObjVal[iter] < STOPPING_OBJ_VAL);

	if(iter < MAX_ITER-1 && bestObjVal[iter] < STOPPING_OBJ_VAL)
		io_printf(IO_STD, "Optimal value is obtained!\n");

	io_printf(IO_BUF, "The final generation:\n--------------------------\n");
	showMyChromosomes();

	/* Setup callbacks */

	/*
	spin1_set_timer_tick(TIMER_TICK_PERIOD);
	spin1_callback_on(TIMER_TICK, hTimer, TIMER_PRIORITY);			// optional, just for debugging
	spin1_callback_on(MCPL_PACKET_RECEIVED, hMCPL, 0);


	if(leadAp) {
		io_printf(IO_BUF, "Initializing router...\n");
		//io_printf(IO_STD, "The leader core is %u\n", myCoreID);
		initRouter();

		// then trigger the population
		io_printf(IO_STD, "Starting iteration-%d\n", iter);
		spin1_schedule_callback(bcastMyChromosomes, 0, 0, BROADCAST_PRIORITY);	// DON'T USE PRIORITY LOWER THAN 1
	}
	else {
		//spin1_delay_us(1000); // give a break for the leader to setup the routing table
		//io_printf(IO_BUF, "Waiting for router initialization...\n");
	}
	*/

	/* Run simulation */
	//spin1_start(SYNC_WAIT);
	spin1_start(SYNC_NOWAIT);
}

