from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import math

import time
import random

def rastriginsfcn(x,y):
    zx = x**2 - (10*math.cos(2*math.pi*x))
    zy = y**2 - (10*math.cos(2*math.pi*y))
    z = 20 + zx + zy
    return z

"""
x=-6:0.01:6;y=-6:0.01:6; z=zeros(length(x),length(y)); 
for i=1:length(x)
    for j=1:length(y), 
         z(i,j)=rastriginsfcn([x(i),y(j)]); 
    end
end
"""

x = np.arange(-5, 5, 0.1)
y = np.arange(-5, 5, 0.1)
Z = np.zeros((len(x), len(y)))
for i in range(len(x)):
    for j in range(len(y)):
        Z[i][j] = rastriginsfcn(x[i], y[j])

fig = plt.figure()
ax = fig.gca(projection='3d')


X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

#plt.show()

gif = plt.figure()
CS = plt.contour(X, Y, Z)

plt.ion()
plt.show()
ax = gif.gca()
for i in range(10):
    d1=[random.random()*10-5]; d2=[random.random()*10-5]
    d3=[random.random()*10-5]; d4=[random.random()*10-5]

    print "generating point [{},{}] and [{},{}]".format(d1,d2,d3,d4)
    ax.plot(d1,d2,'ro')
    ax.plot(d3,d4,'bo')
    plt.draw()
    time.sleep(1)
    gif.canvas.flush_events()
    ax.lines.pop()
    ax.lines.pop()

raw_input("The press enter...")
