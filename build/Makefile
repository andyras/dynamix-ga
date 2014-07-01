# name of main executable
### This line should be just about all you have to change ###
BIN = dynamix-ga
#############################################################

#DEBUG_MAKE = "foo"

# directories for headers, objects and source files
INCDIR = ./include
OBJDIR = ./obj
SRCDIR = ./src
BINDIR = ../bin

# compiler
CPP = $(HOME)/bin/mpich/bin/mpic++
CPP = mpic++

CPPFLAGS = -O3 -Wall -std=c++11 -fopenmp
INCLUDES = -I$(INCDIR)

# this should be first in linking order because it depends on sundials
# dynamix source directory
DYNAMIXDIR = $(HOME)/git/dynamix/dm/build
INCLUDES += -I$(DYNAMIXDIR)/include -I$(HOME)/bin/include
LDFLAGS += -L$(DYNAMIXDIR)/lib -ldynamix -L$(HOME)/bin/lib

# GALib-mpi source directory
GALIBDIR = $(HOME)/git/GALib-mpi
INCLUDES += -I$(GALIBDIR)
#LDFLAGS += -L$(GALIBDIR)/ga-mpi -l/home/andyras/git/GAlib-mpi/ga-mpi/libga-mpi.a
LDFLAGS += -L$(GALIBDIR)/ga-mpi -lga-mpi

LDFLAGS += -lsundials_cvode -lsundials_nvecserial

# basic compiler-dependent flags
ifeq ($(CPP),icpc)
  CPPFLAGS += -fast -xHOST -no-prec-div -mkl -no-multibyte-chars
  LDFLAGS += -mkl
else # g++*
  LDFLAGS += -liomp5
  LDFLAGS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
  ifeq ($(shell hostname),tim.selfip.org)
    LDFLAGS += -lpthread -lm
  endif
endif

ifdef DEBUG_MAKE
  CPPFLAGS += -pg -g -debug
  LDFLAGS += -pg -g -debug
endif

# make will look for .cpp files in $(SRCDIR)
vpath %.cpp $(SRCDIR)

# list of source files
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# object files have same name as .cpp files, but with .o extension
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,obj/%.o,$(SOURCES))

# build the main executable; this should be listed first
$(BIN): $(OBJECTS)
	$(CPP) -o $@ $^ $(LDFLAGS) -fopenmp

# automatic rule for building objects
$(OBJDIR)/%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@ $(INCLUDES)

# all objects depend on object directory
$(OBJECTS): | $(OBJDIR)

.PHONY: clean install uninstall $(BINDIR) $(OBJDIR)

clean:
	rm -rf $(OBJDIR)
	rm -f $(BIN)

install: $(BIN) | $(BINDIR)
	cp $(BIN) $(BINDIR)

uninstall:
	rm -rf $(BINDIR)/$(BIN)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p ../bin
