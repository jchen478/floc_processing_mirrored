SRCS := $(wildcard *.cu)
OBJS := $(patsubst %.cu,%.o,$(SRCS))

COMPILER = nvcc 
NVCCFLAGS = -arch=sm_60  
LFLAGS = -lm 
EXE = clustering

$(EXE): $(OBJS) 
	$(COMPILER) $(NVCCFLAGS) $(LFLAGS) -o $@ $^

%.o: %.cu
	$(COMPILER) $(NVCCFLAGS) $(LFLAGS) -c $<

clean:
	-rm -f $(OBJS) $(EXE)

