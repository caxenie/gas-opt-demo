/*
	Sample demo software for Genetic Algorithms Exercise
	Simple nonlinear function maximisation using GAs.

	To be used as skeleton in the CI class @ TUM WSS 2012-13
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
  #define M_PI 3.1415926535897932384626433832795
#endif

/* upper and lower interval limits for the candidate function */
#define X_MIN 	-1.0 
#define X_MAX 	2.0 	

/* GAs parameters */
#define POPULATION_SIZE 1000  // chromosomes
#define CHROMOSOME_SIZE 22    // bits
#define MAX_EPOCHS        100

/////////////////////////////// PROBLEM SPECIFIC CODE //////////////////////////////////

/* the candidate function to be maximised */
double f(double x){
  return (double)(x*sin(10*x*M_PI)+1.0);
}

/////////////////////////////END OF PROBLEM SPECIFIC CODE /////////////////////////////


//////////////////////////////////// GAS SPECIFIC CODE ///////////////////////////////////////////

/* chromosome abstraction */
struct chromosome{
  /* number of bits in the representation */
  int n;
  /* binary chromosome genes - bits */
  int* b;
  /* converted binary to decimal */
  long x_star;
};

/* chromosomes population */
struct population{
  /* size in chromosomes */
  int size;
  /* chromosomes in the population */
  struct chromosome *c;
};

/* initialize a chromosome */
void init_chromosome(struct chromosome c, int n){
  c.n = n;
  c.b = (int*)calloc(c.n, sizeof(int));
  /* setup the genes randomly */
  for(int i=0; i<c.n; ++i)
    c.b[i] = (((double) rand () / (double) RAND_MAX)<0.5)?0:1;

  /* compute the decimal representation of the genes */
  for(int i=0; i<c.n; ++i)
    c.x_star += (c.b[i] * 2^i);
}

/* encoding function from binary chromosome to a decimal value using
	- the size of the representation
	- a min and max values (interval of values)
	- a conveted base 2 to base 10 value of the chromosome
*/
double encode_chromosome(struct chromosome c){
  return (double)(X_MIN + c.x_star*((X_MAX - X_MIN)/(pow(2, c.n) - 1)));
}


/* compute fitness value for a given chromosome - evaluation function, e */
double compute_chromosome_fitness(struct chromosome c){
  return (double)f(encode_chromosome(c));
}

/* initialize a chromosomes population with given parameters */
void init_population(struct population p, int psize, int csize){
  p.size = psize;
  p.c = (struct chromosome*)calloc(psize, sizeof(struct chromosome));
  for(int i=0; i<p.size; ++i)
    init_chromosome(p.c[i], csize);

}

/* perform xover with 2 chromosomes resulting 2 chromosomes after a xover point */
struct chromosome* compute_crossover(struct chromosome c1, struct chromosome c2, short xover_pt){
  struct chromosome *ret = (struct chromosome*)calloc(2, sizeof(struct chromosome));
  int aux_gene = 0;
  for(int i=0; i<c1.n; ++i){
      if(i>=xover_pt){
          aux_gene = c1.b[i];
          c1.b[i] = c2.b[i];
          c2.b[i] = aux_gene;
        }
    }
  ret[0] = c1;
  ret[1] = c2;
  return ret;
}

/* perform a mutation over one/more genes with a probability equal to the mutation rate */
struct chromosome compute_mutation(struct chromosome c, int mg){
      struct chromosome mc;
      init_chromosome(mc, c.n);
      mg = rand()%c.n;
      (c.b[mg]==0)?1:0;
}

/* select the fittest chromosome in the population */
struct chromosome select_fittest(struct population p)
{
  double  best_val = 0.0;
  int best_idx =0;
  struct chromosome best;
  for(int i=0;i<p.size;i++){
    if(compute_chromosome_fitness(p.c[i])> best_val)
      best_val = compute_chromosome_fitness(p.c[i]);
      best_idx = i;
  }
  best = p.c[best_idx];
  init_population(p, p.size, p.c->n);
  p.c[0] = best;
  return best;
}

/* entry point */
int main(int argc, char* argv[]){

  /* the problem is to determine the maximum of f(x) = x*sin(10*pi*x) + 1
           in the interval xe[-1, 2]
                - analytically we have an infinite number of solutions for this equation
                - 3 types / forms of solutions
                      |	xi = (2*i-1)/20 + ei, i = 1,2....
                      | xi = 0
                      | xi = (2*i+1)/20 + ei, i = 1,2....
                - we split the input interval in 3 x 1000000 equal sized ranges
                - this means we have a 22-bit representation for chromosomes
        */

  // iterator through epochs
  int epoch = 0;
  // initial population of chromosomes
  struct population p;
  // the best / fittest chromosome after evolution
  struct chromosome best;
  // xover output
  struct chromosome *xover_out = (struct chromosome*)calloc(2, sizeof(struct chromosome));
  // initialize population
  init_population(p, POPULATION_SIZE, CHROMOSOME_SIZE);


  // loop to evolve
  while(++epoch < MAX_EPOCHS){
      for(int i=1; i<p.size; ++i){
//          // select
//          select_fittest(p);
//          // crossover
//          xover_out = compute_crossover(p.c[rand()%p.size], p.c[rand()%p.size], 5);
//          // mutate
//          p.compute_mutation()
          // evaluate


        }
  }

  printf("The evolution finished.\n The fittest individual gives a maximum of %lf\n", compute_chromosome_fitness(best));
  return EXIT_SUCCESS;
}
