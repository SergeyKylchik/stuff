all: a.out

a.out: main_copy.o functions_copy.o solve_copy.o
	g++ -pthread -O3 -std=c++11 main_copy.o functions_copy.o solve_copy.o -o a.out

main_copy.o: main_copy.cpp
	g++ -pthread -O3 -c main_copy.cpp

functions_copy.o: functions_copy.cpp
	g++ -pthread -O3 -c functions_copy.cpp

solve_copy.o: solve_copy.cpp
	g++ -pthread -O3 -c solve_copy.cpp

clean:
	rm -rf *.o a.out
