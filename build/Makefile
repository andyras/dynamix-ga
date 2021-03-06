### This line should be just about all you have to change ######################
# name of main executable
BIN = dynamix-ga
################################################################################

#DEBUG_MAKE = "foo"
#__USE_MKL__ = yes

ifdef DEBUG_MAKE
  CPPFLAGS += -pg -g -debug
  LDFLAGS += -pg -g -debug
endif

ifdef __USE_MKL__
  CPPFLAGS += -D__USE_MKL__
  LDFLAGS += -liomp5
  LDFLAGS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
else
  LDFLAGS += -lblas
endif

# directories for headers, objects and source files
INCDIR = ./include
OBJDIR = ./obj
SRCDIR = ./src
BINDIR = ../bin

# compiler
CPP = mpic++

CPPFLAGS = -O0 -Wall -std=c++11 -fopenmp
INCLUDES = -I$(INCDIR)

################################################################################
# this is needed for included header files from dynamix
CPPFLAGS += -D__BOOST_SERIALIZE__
################################################################################

# this should be first in linking order because it depends on sundials
LDFLAGS += -ldynamix

LDFLAGS += -lga-mpi

LDFLAGS += -lsundials_cvode -lsundials_nvecserial
LDFLAGS += -lboost_mpi -lboost_serialization

# basic compiler-dependent flags
COMPILER = $(word 1, $(shell mpic++ -show))
ifeq ($(COMPILER),icpc)
  CPPFLAGS += -xHOST -no-prec-div -mkl -no-multibyte-chars
  LDFLAGS += -mkl
else # g++*
  ifeq ($(shell hostname),tim.selfip.org)
    LDFLAGS += -lpthread -lm
  endif
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
