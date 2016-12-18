CXXFLAGS = -Wall -g # version with the blas lib ...
#CFLAGS = -Wall -g   # version without the blas lib ...
#LDFLAGS =   -lblas -lm # pour Linux 1
#LDFLAGS =  -lopenblas -lm # pour linux 2
PROGS= lap2d lap3d
all:$(PROGS)

mainGC: mainGC.o GC.o 
	$(CXX) -o $@ $^  $(CXXFLAGS) $(LDFLAGS)
#  <=> $(CXX) -o mainGC  mainGC.o GC.o  $(CFLAGS) $(LDFLAGS)
lap2d: lap2d.o GC.o gnuplot.o
	$(CXX) -o $@ $^  $(CXXFLAGS) $(LDFLAGS)
lap3d: lap3d.o GC.o 
	$(CXX) -o $@ $^  $(CXXFLAGS) $(LDFLAGS)
clean:
	-rm *.o  $(PROGS) lap1d *~ *.exe *.txt

# les depandences	
mainGC.o GC.o: GC.hpp
gnuplot.o: gnuplot.hpp 
