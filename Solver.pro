TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    solver.cpp \
    cell.cpp \
    param.cpp

HEADERS += \
    solver.h \
    cell.h \
    header.h \
    param.h

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../Libraries/OSG/OSG-3.0.1-x86-shared/lib/ -losg
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../Libraries/OSG/OSG-3.0.1-x86-shared/lib/ -losgd

INCLUDEPATH += $$PWD/../../../Libraries/OSG/OSG-3.0.1-x86-shared/include
DEPENDPATH += $$PWD/../../../Libraries/OSG/OSG-3.0.1-x86-shared/include
