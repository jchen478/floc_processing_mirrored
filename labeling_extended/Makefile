SRCS := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS))

COMPILER = icc 
CFLAGS = -g -I$(MKLROOT)/include
MKLLIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl
LFLAGS = -L/$(MKLROOT)/lib/intel64 -lmkl_core $(MKLLIBS)
EXE = labeling_extended 

$(EXE): $(OBJS) 
	$(COMPILER) $(CFLAGS) $(LFLAGS) -o $@ $^

%.o: %.cpp
	$(COMPILER) $(CFLAGS) $(LFLAGS) -c $<

clean:
	-rm -f $(OBJS) $(EXE)

