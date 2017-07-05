TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    main.c \
    MT19937.c

INCLUDEPATH += \
    /opt/spinnaker_tools_3.1.0/include

DISTFILES += \
    Makefile \
    ../myAlgorithm.txt
