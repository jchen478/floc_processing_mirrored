SRCS := $(wildcard *.cu)
OBJS := $(patsubst %.cu,%.o,$(SRCS))

COMPILER = nvcc -g -G 
NVCCFLAGS = -arch=sm_60 -lm 
EXE = sizeChar 

$(EXE): $(OBJS) 
	$(COMPILER) $(NVCCFLAGS) -o $@ $^

%.o: %.cu
	$(COMPILER) $(NVCCFLAGS) -c $<

clean:
	-rm -f $(OBJS) $(EXE)

