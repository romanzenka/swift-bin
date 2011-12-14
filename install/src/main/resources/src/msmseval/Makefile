#makefile for msmsEval

EXECUTABLE = ./bin/msmsEval
LINKCC = $(CXX)

#CXXFLAGS denotes flags for the C++ compiler

CXX = g++
CXXFLAGS = -O3 -march=pentium4
LDFLAGS = -lexpat

SRCS := $(wildcard ./src/*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(wildcard ./src/*.cpp))
#DEPS := $(patsubst %.o,%.d,$(OBJS))


all: $(EXECUTABLE)

	@echo ""
	@echo "By using msmsEval you agree to the terms and conditions set out in the"
	@echo "LICENSE file included with this package. If you have not done so, please"
	@echo "refer to the LICENSE file before use of this application."
	@echo ""

#define the components of the program, and how to link them
#these components are defined as dependencies; that is they must be up-to-date before the code is linked

$(EXECUTABLE): $(DEPS) $(OBJS)
	$(LINKCC) $(CXXFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDFLAGS)
	
#specify the dep files depend on the cpp files

%.d: %.cpp
	$(CXX) -M $(CXXFLAGS) $< > $@
	$(CXX) -M $(CXXFLAGS) $< | sed s/\\.o/.d/ > $@
	


clean:
	-rm $(OBJS)  *~

explain:
	@echo "The following info represents the program:"
	@echo "Final exec name: $(EXECUTABLE)"
	@echo "Source files:       $(SRCS)"
	@echo "Object files:       $(OBJS)"
#	@echo "Dep files:          $(DEPS)"

depend: $(DEPS)
	@echo "Deps are now up-to-date."
 	
#-include $(DEPS)
