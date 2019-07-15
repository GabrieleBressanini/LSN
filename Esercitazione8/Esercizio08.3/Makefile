CC = g++
CFLAGS = -Wall -O3 --std=c++11

qmc1d : qmc1d.o random.o
	$(CC) random.o qmc1d.o -o qmc1d
qmc1d.o : qmc1d.cpp random.h constants.h functions.h
	$(CC) -c qmc1d.cpp -o qmc1d.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o seed.out qmc1d potential.dat kinetic.dat probability.dat

