CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
