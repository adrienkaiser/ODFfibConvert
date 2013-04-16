# target: dependencies
# [tab] system command

all: Fib

main.o: main.cpp zlib.h zconf.h
	g++ -c main.cpp

Fib: main.o
	g++ main.o -o Fib -lz
	rm main.o
