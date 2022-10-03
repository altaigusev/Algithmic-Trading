QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    approximation.cc \
    interpolation.cc \
    loader.cc \
    main.cc \
    mainwindow.cc \
    qcustomplot.cpp

HEADERS += \
    approximation.h \
    interpolation.h \
    loader.h \
    mainwindow.h \
    qcustomplot.h \
    s21_matrix_oop.h \
    s21_matrix_oop.inl

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    CVX.csv

ICON = stocks.png
