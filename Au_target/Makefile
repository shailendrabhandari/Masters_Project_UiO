Target  = runAnalysis
Objects = Event.o 

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

all:$(Target)

runAnalysis: Analysis.C $(Objects)
	c++ -o $@ Analysis.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Event.o: Event.C Event.h
	c++ -o $@ Event.C -c $(ROOTFLAGS)

clean:
	@rm -rf $(Target) *.o
