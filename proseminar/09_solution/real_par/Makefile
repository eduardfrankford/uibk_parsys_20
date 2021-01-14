CC=gcc
CC_FLAGS=-O3 -g3

.PHONEY: all
all: real

SRC=$(wildcard *.c)

real: $(SRC)
	@$(CC) $(CC_FLAGS) $(SRC) -fopenmp -o real -g -pg -lm

.PHONEY: clean
clean:
	@rm real

.PHONEY: run
run: all
	@echo "Parallel:"
	@./real
	@echo
