BOOST_DIR = /opt/local
BOOST_CFLAGS = -I$(BOOST_DIR)/include
BOOST_LDLIBS = -L$(BOOST_DIR)/lib -lboost_system-mt
EIGEN_CFLAGS = -I/opt/local/include/eigen3

OMPL_DIR = /Users/A.I/Documents/ompl-1.1.0-Source
OMPL_CFLAGS = -I$(OMPL_DIR)/src
OMPL_LDLIBS = -L$(OMPL_DIR)/../build-ompl-1.1.0-Source-Desktop_Qt_5_9_1_clang_64bit-Default/lib -lompl

# to override these variables, call make OMPL_DIR="/path/to/ompl" OMPL_LDLIBS="..."

# for when $PWD is a symlink:
PARENT_DIR = $(shell sh -c 'cd $$PWD/..; pwd')

# -std=c++11 -stdlib=libc++

CXXFLAGS = -ggdb -O2 -I$(PARENT_DIR)/include -Wall -Wno-unused -Wno-overloaded-virtual -Wno-sign-compare -fPIC -std=c++11 -stdlib=libc++  $(BOOST_CFLAGS) $(OMPL_CFLAGS) $(EIGEN_CFLAGS)
LDLIBS = -ggdb -O2 -lpthread -ldl $(BOOST_LDLIBS) $(OMPL_LDLIBS)

.PHONY: clean all install doc

OS = $(shell uname -s)
ifeq ($(OS), Linux)
	CFLAGS += -D__linux
	EXT = so
	INSTALL_DIR ?= $(PARENT_DIR)/..
else
	CFLAGS += -D__APPLE__
	EXT = dylib
	INSTALL_DIR ?= $(PARENT_DIR)/../vrep.app/Contents/MacOS/
endif

all: libv_repExtOMPL.$(EXT) doc

# reference.html: callbacks.xml callbacks.xsl
# saxon -s:$< -a:on -o:$@

v_repExtOMPL.o: stubs.h

stubs.o: stubs.h stubs.cpp

stubs.h: callbacks.xml
	python -m v_repStubsGen -H $@ $<

stubs.cpp: callbacks.xml
	python -m v_repStubsGen -C $@ $<

libv_repExtOMPL.$(EXT): v_repExtOMPL.o stubs.o $(PARENT_DIR)/common/v_repLib.o
	$(CXX) $^ $(LDLIBS) -shared -o $@

clean:
	rm -f libv_repExtOMPL.$(EXT)
	rm -f *.o
	rm -f stubs.cpp stubs.h
	rm -f reference.html

install: all
	cp libv_repExtOMPL.$(EXT) $(INSTALL_DIR)
