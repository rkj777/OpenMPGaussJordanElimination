main: timer.h Lab3IO.h Lab3IO.c Main.c
	gcc -Wall -fopenmp Main.c Lab3IO.c -o Main

clean:
	rm Main
