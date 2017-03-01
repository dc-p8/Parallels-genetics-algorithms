// mpicc -O3 main_mpi.c -o mpi
// mpicc -np 4 ./mpi 10000 50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#define MCW MPI_COMM_WORLD

typedef unsigned long int uli;

typedef struct problem{
	int msize;
	int mrank;
	int nb_bits;
	int size_bitarray;
	int population_size;
	int seed;
	uli *solution;
	uli **bitarrays;
	int *distances;
	int *sorted_indexes;
}problem;

problem p;

void shuffle(uli **one, int nb_bits, int *seed)
{
	for(int i = 0; i < nb_bits; i++)
	{
		short bi = ((*one)[i / 64] & ((0x1UL) << (i % 64))) != 0;
		int k = rand_r(seed) % (nb_bits);
		
		short bk = ((*one)[k / 64] & ((0x1UL) << (k % 64))) != 0;
		if(bk != bi)
		{
			(*one)[k / 64] ^= (0x1UL << (k % 64));
			(*one)[i / 64] ^= (0x1UL << (i % 64));
		}
	}
}


void eval_pop()
{
	for(int i = 0; i < p.population_size; i++)
	{
		p.distances[i] = 0;
		for(int j = 0; j < p.size_bitarray; j++)
		{
			uli diff = p.bitarrays[i][j] ^ p.solution[j];
			for(int k = 0; k < (j == p.size_bitarray - 1 ? p.nb_bits % 64 : 64); k++)
			{
				if( (diff & (0x1UL << (k) )) != 0 )
					p.distances[i]++;
			}
		}
	}
}

int compare(const void *arg_a, const void *arg_b)
{
	int *i1 = (int *)arg_a;
	int *i2 = (int *)arg_b;
	return p.distances[*i1] > p.distances[*i2]; 
}

void gen(int nb_bits, int pop_size)
{
	MPI_Comm_rank(MCW, &p.mrank);
	MPI_Comm_size(MCW, &p.msize);

	p.nb_bits = nb_bits;
	p.size_bitarray = p.nb_bits / 64 + ((p.nb_bits % 64) != 0);
	p.solution = NULL;
	p.bitarrays = NULL;
	p.distances = NULL;
	p.sorted_indexes = NULL;
	p.population_size = pop_size;
	p.seed = p.mrank;

	p.solution = malloc(sizeof(uli) * p.size_bitarray);
	p.distances = malloc(sizeof(int) * p.population_size);
	p.bitarrays = malloc(sizeof(uli *) * p.population_size);
	p.sorted_indexes = malloc(sizeof(int) * p.population_size);
	
	int seed_solution = 0;

	int randn;
	randn = nb_bits / 5;
	for(int i = 0; i < p.size_bitarray; i++)
		p.solution[i] = 0x0UL;
	for(int i = 0; i < randn; i++)
		p.solution[i / 64] = (p.solution[i / 64] | 0x01UL << (i % 64));

	shuffle(&p.solution, p.nb_bits, &seed_solution);

	for(int i = 0; i < p.population_size; i++)
	{
		p.distances[i] = 0;
		p.bitarrays[i] = malloc(sizeof(uli) * p.size_bitarray);

		for(int j = 0; j < p.size_bitarray; j++)
			p.bitarrays[i][j] = 0x0UL;

		randn = (rand_r(&(p.seed)) % p.nb_bits);
		for(int j = 0; j < randn; j++)
			p.bitarrays[i][j/64] |= 0x01UL << (j % 64);

		shuffle(&p.bitarrays[i], p.nb_bits, &(p.seed));

	    p.sorted_indexes[i] = i;
	}
}

void mutate()
{
	for(int i = 1; i < p.population_size; i++)
	{
		int k = rand_r(&p.seed) % p.nb_bits;
		p.bitarrays[p.sorted_indexes[i]][k / 64] ^= (0x1UL << (k % 64));
	}
}


int select_pop()
{
	int best = p.nb_bits;
	int best_i = 0;
	int bests[p.msize];
	MPI_Allgather(&p.distances[p.sorted_indexes[0]], 1, MPI_INT, bests, 1, MPI_INT, MCW);
	for(int i = 0; i < p.msize; i++)
	{
		if(bests[i] < best)
		{
			best_i = i;
			best = bests[i];
		}
	}

	MPI_Bcast(p.bitarrays[p.sorted_indexes[0]], p.size_bitarray, MPI_UNSIGNED_LONG, best_i, MCW);
	
	p.distances[p.sorted_indexes[0]] = best;

	int nb = p.population_size / 4;
	int nb2 = p.population_size  / 10;
	for(int i = 0; i < nb; i++)
	{
		for(int j = 0; j < p.nb_bits; j++)
		{
			if(p.bitarrays[p.sorted_indexes[rand_r(&p.seed) % nb2]][j / 64] & (0x01UL << (j % 64)))
				p.bitarrays[p.sorted_indexes[p.population_size - i - 1]][j / 64] |= (0x01UL << (j % 64));
			else
				p.bitarrays[p.sorted_indexes[p.population_size - i - 1]][j / 64] &= ~(0x01UL << (j % 64));
		}
	}
	return best;
}


int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	if(argc != 3)
	{
		MPI_Finalize();
		return 0;
	}

	int nb_bits = atoi(argv[1]);
	int pop_size_per_process = atoi(argv[2]);

	if(nb_bits < pop_size_per_process || nb_bits > 100000000 || pop_size_per_process < 10)
	{
		MPI_Finalize();
		return 0;
	}

	gen(nb_bits, pop_size_per_process);

	if(!p.mrank)
		printf("nb_bits : %d\npop_size_per_process : %d\nprocessus : %d\n", nb_bits, pop_size_per_process, p.msize);

	int gen = 0;
	while(1)
	{
		eval_pop();
		qsort(p.sorted_indexes, p.population_size, sizeof(int), compare);
		if(!select_pop() || gen++ >= 1000000)
			break;

		mutate();
		if(!(gen % 100))
		{
			printf("Gen : %d. Best : %d\n", gen, p.distances[p.sorted_indexes[0]]);
		}
	}

	FILE *fp = NULL;
	char filename[20];
	sprintf(filename, "log_%d.out", p.mrank);

	fp = fopen(filename, "w");

	fprintf(fp, "PROCESS %d : --- GENERATION : %d\n", p.mrank, gen);
	fprintf(fp, "SOLUTION :\n       | ");
	for(int i = 0; i < p.nb_bits; i++)
		fprintf(fp, "%c", ' ' + ((p.solution[i / 64] & (0x1UL << (i % 64))) != 0));
	
	fprintf(fp, " |\n");
	for(int i = 0; i < p.population_size; i++)
	{
		fprintf(fp, "%4d : | ", i);
		for(int k = 0; k < p.nb_bits; k++)
		{
			fprintf(fp, "%c", ' ' + ((p.bitarrays[p.sorted_indexes[i]][k / 64] & (0x1UL << (k % 64))) != 0));
		}
		fprintf(fp, " | %4d |\n", p.distances[p.sorted_indexes[i]]);
	}
	fprintf(fp, "\n");


	fclose(fp);
	for(int i = 0; i < p.population_size; i++)
		free(p.bitarrays[i]);
	free(p.solution);
	free(p.distances);
	free(p.sorted_indexes);
	MPI_Finalize();
	return 0;
}