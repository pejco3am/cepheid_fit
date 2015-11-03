CXX= g++
CPPFLAGS= -Wall -O3 -fopenmp
LDFLAGS= -Wall -O3 -fopenmp

OBJS=initialize.o data_io.o fitting.o data_managment.o main.o

fit_single: $(OBJS)
	$(CXX) -o fit_single.exe $(OBJS) $(LDFLAGS) cmpfit-1.2/mpfit.o

mpfit:
	$(MAKE) -C cmpfit-1.2

clean:
	rm -f *.o
	rm -f cmpfit-1.2/*.o

.cpp.o:
	$(CXX) -c $(CPPFLAGS) $<	
