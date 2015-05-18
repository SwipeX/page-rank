pagerank-serial:
	mpicc pagerank-serial.c -o pagerank-serial
profile:
	mpicc -pg pagerank-serial.c -o pagerank-serial

clean:
	rm -f pagerank-serial