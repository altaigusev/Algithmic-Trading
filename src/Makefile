.PHONY: all clean check gcov_report dvi install uninstall dist
XX=g++
APP=$(BD)Trading.app

FILES= approximation.cc \
    interpolation.cc \
    loader.cc \
    main.cc \
    mainwindow.cc \
    qcustomplot.cpp

OS := $(shell uname -s)
BD= build/

ifeq ($(OS), Darwin)
	CC += -D OS_MAC
else
	CC += -D OS_LINUX
endif

all: build install

build: clean copy
	cd $(BD); qmake Trading.pro
	make -C $(BD)

copy:
	@mkdir build
	cp *.inl *.pro *.cpp *.cc *.png *.h *.ui $(BD)

install: uninstall build
	@mkdir $(HOME)/Applications/Trading/ 2>/dev/null
	cp -rf $(APP) $(HOME)/Applications/Trading/
	make clean

uninstall: 
	rm -rf $(HOME)/Applications/Trading

dist:
	rm -rf Trading
	mkdir Trading
	cp Makefile *.pro *.cpp *.user Trading/
	rsync -rupE $(DIRS)  Trading/
	tar -cvzf $(HOME)/Desktop/Trading.tgz Trading/
	rm -rf Trading/

dvi:
	open Doc/index.html

clean:
	rm -rf build test *.a *.o *.cfg *.gcda *.gcno *.html *.info *.dSYM report
