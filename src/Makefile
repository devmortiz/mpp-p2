CC = gcc
MPICC = mpicc
CFLAGS = -O3 -Wall -std=c99 -g

N = 
M = 
N_GEN = 
T_POB = 
Numpro = 3


EXEC = sec
MPI = mpi
OUTFILE = ../output/out
C_FILES = main.c io.c mh.c

sec: $(C_FILES)
	$(CC) $(CFLAGS) $(C_FILES) -o $(EXEC) -DTIME -DDEBUG
mpi: $(C_FILES)
	$(MPICC) $(CFLAGS) $(C_FILES) -o $(MPI) -DTIME -DDEBUG

test_sec:
	./$(EXEC) $(N) $(M) $(N_GEN) $(T_POB) > $(OUTFILE)_$(N)_$(M)_$(N_GEN)_$(T_POB).txt
test_mpi:
	mpirun -np $(Numpro) $(MPI) $(N) $(M) $(N_GEN) $(T_POB) > $(OUTFILE)_$(N)_$(M)_$(N_GEN)_$(T_POB).txt

valgrind_sec:
	valgrind --tool=memcheck --leak-check=yes --track-origins=yes --log-file=../valgrind.output ./$(EXEC) $(N) $(M) $(N_GEN) $(T_POB) > $(OUTFILE)_$(N)_$(M)_$(N_GEN)_$(T_POB).txt

valgrind_mpi:
	valgrind --tool=memcheck --leak-check=yes --track-origins=yes --log-file=../valgrind.output mpirun -np $(Numpro) $(MPI) $(N) $(M) $(N_GEN) $(T_POB) > $(OUTFILE)_$(N)_$(M)_$(N_GEN)_$(T_POB).txt


clean:
	rm -f $(EXEC) $(MPI)
