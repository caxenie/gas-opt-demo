/*
	Sample demo software for Genetic Algorithms Exercise
	Simple nonlinear function maximisation using GAs.

	To be used as skeleton in the CI class @ TUM WSS 2012-13
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/* GAs parameters */
#define REPRESENTATION_SIZE 22		  // bits
#define POPULATION_SIZE     50            // chromosomes
#define MAX_GENERATIONS     150         // number of generations to evolve
#define XOVER_PROB          0.8           // crossover probability
#define MUTATION_PROB       0.15          // mutation probability

/////////////////////////////// PROBLEM SPECIFIC CODE //////////////////////////////////


/* upper and lower interval limits for the candidate function */
#define X_MIN 	-1.0
#define X_MAX 	2.0

/* log data */
FILE *fp;

/* the candidate function to be maximised */
double f(double x){
  return (double)(x*sin(10*x*M_PI)+1.0);
}

/////////////////////////////END OF PROBLEM SPECIFIC CODE /////////////////////////////


//////////////////////////////////// GAS SPECIFIC CODE ///////////////////////////////////////////

/* chromosome abstraction */
struct chromosome{
  /* binary chromosome representation */
  int *b;
  /* converted binary to decimal value encoded */
  double x_star;
  /* the fitness value of the chromosome */
  double fitness;
  /* relative fitness */
  double rfitness;
  /* cumulative fitness */
  double cfitness;
};

/* chromosomes population */
struct population{
  /* size in chromosomes */
  int size;
  /* chromosomes in the population */
  struct chromosome *c;
  /* current evolution generation */
  int gen;
  /* index of the fittest chromosome */
  int best_chromosome_idx;
};

/* random value generator within the problem bounds */
double randomize(double l, double h)
{
  return (((double)(rand()%1000)/1000.0)*(h - l) + l);
}

/* mapping / encoding function from binary to decimal given a binary string */
double encode_binary_chromosome(int *b)
{
  double ret = 0.0f;
  for(int i=0;i<REPRESENTATION_SIZE; ++i){
      ret+=b[i]*pow(2, i);
    }
  return ret;
}

/* initialize a chromosome */
void init_chromosome(struct chromosome *c){
  /* create a binary chromosome representation with random genes */
  c->b = (int*)calloc(REPRESENTATION_SIZE, sizeof(int));
  for(int i=0;i<REPRESENTATION_SIZE;++i){
      c->b[i] = ((rand()%1000/1000.0)<0.5)?0:1;
    }
  /* compute the decimal representation of the genes */
  c->x_star = encode_binary_chromosome(c->b); //randomize(X_MIN, X_MAX);
  /* fitness values */
  c->fitness = 0.0f;
  c->rfitness = 0.0f;
  c->cfitness = 0.0f;
}

/* initialize a chromosomes population with given parameters */
void init_population(struct population *p, int psize){
  p->size = psize;
  p->gen = 0;
  p->best_chromosome_idx = 0;
  p->c = (struct chromosome*)calloc(p->size, sizeof(struct chromosome));
  for(int i=0; i<p->size; i++){
      init_chromosome(&p->c[i]);
    }
}

/* evaluate function, takes a user defined function and computes it for every chromosome */
void evaluate_population(struct population *p)
{
  for(int i=0; i<p->size; i++){
      p->c[i].fitness = f(p->c[i].x_star);
    }
}

/* select the best (fittest) chromosome in the population */
void select_best(struct population *p)
{
  p->best_chromosome_idx = 0;
  for(int i=0; i<p->size; ++i){
      /* the last entry in the population is the best chromosome */
      if(p->c[i].fitness > p->c[POPULATION_SIZE].fitness){
          p->best_chromosome_idx = i;
          p->c[POPULATION_SIZE].fitness = p->c[i].fitness;
        }
    }
  /* found the fittest then copy the genes */
  p->c[POPULATION_SIZE].x_star = p->c[p->best_chromosome_idx].x_star;
}

/* apply elitism so that if the previous best chromosome is better than the
 * current generation best the first will replace the worst chromosome in the
 * current generation.
 */
void apply_elitism(struct population *p)
{
  struct chromosome *best = (struct chromosome*)calloc(1, sizeof(struct chromosome));
  struct chromosome *worst= (struct chromosome*)calloc(1, sizeof(struct chromosome));
  int best_idx = 0, worst_idx = 0;
  init_chromosome(best);
  init_chromosome(worst);
  best->fitness = p->c[0].fitness;
  worst->fitness = p->c[0].fitness;

  for(int i=0;i< p->size-1;++i){
      if(p->c[i].fitness > p->c[i+1].fitness){
          if(p->c[i].fitness >= best->fitness){
              best->fitness = p->c[i].fitness;
              best_idx = i;
            }
          if(p->c[i+1].fitness <= worst->fitness){
              worst->fitness = p->c[i+1].fitness;
              worst_idx = i+1;
            }
        }
      else{
          if(p->c[i].fitness <= worst->fitness){
              worst->fitness = p->c[i].fitness;
              worst_idx = i;
            }
          if(p->c[i+1].fitness >= best->fitness){
              best->fitness = p->c[i+1].fitness;
              best_idx = i+1;
            }
        }
    }
  /* if best chromosome from the new population is better than */
  /* the best chromosome from the previous population, then    */
  /* copy the best from the new population; else replace the   */
  /* worst chromosome from the current population with the     */
  /* best one from the previous generation                     */
  if(best->fitness >= p->c[POPULATION_SIZE].fitness){
      p->c[POPULATION_SIZE].x_star = p->c[best_idx].x_star;
      p->c[POPULATION_SIZE].fitness = p->c[best_idx].fitness;
    }
  else{
      p->c[worst_idx].x_star = p->c[POPULATION_SIZE].x_star;
      p->c[worst_idx].fitness = p->c[POPULATION_SIZE].fitness;
    }
}

/* selection function using the elitist model in which only the
 * best chromosome survives - Winner-Take-All
 */
void apply_selection(struct population *p, struct population *newp)
{
  double sum_fit = 0.0f;
  double prob = 0.0f;
  /* find the global value of the fitness of the population */
  for(int i=0; i< p->size; ++i){
      sum_fit+=p->c[i].fitness;
    }
  /* compute the relative fitness of the population */
  for(int i=0; i<p->size; ++i){
      p->c[i].rfitness = p->c[i].fitness/sum_fit;
    }
  p->c[0].cfitness = p->c[0].rfitness;

  /* compute the cumulative fitness of the population */
  for(int i=1; i< p->size; ++i){
      p->c[i].cfitness = p->c[i-1].cfitness + p->c[i].rfitness;
    }
  /* select the survivors using the cumulative fitness */
  for(int i=0;i<p->size;++i){
      prob = rand()%1000/1000.0;
      if(prob < p->c[0].cfitness)
        newp->c[i] = p->c[0];
      else
        {
          for(int j=0; j<p->size; ++j){
              if(prob>=p->c[j].cfitness && prob<p->c[j+1].cfitness)
                newp->c[i] = p->c[j+1];
            }
        }
    }
  /* one the new population is created copy it back in the working var */
  for(int i=0 ;i<p->size; ++i)
    p->c[i] = newp->c[i];
}

/* apply the single point crossover operator which takes 2 parents */
void apply_crossover(struct population *p)
{
  /* counter of members chosen */
  int cnt = 0;
  /* probability to xover */
  double prob_xover = 0.0f;
  /* the two parent containers init */
  struct chromosome *p1 = (struct chromosome*)calloc(1, sizeof(struct chromosome));
  init_chromosome(p1);
  /* cross over loop */
  for(int i=0; i< p->size; ++i){
      prob_xover = rand()%1000/1000.0;
      if(prob_xover < XOVER_PROB){
          cnt++;
          if(cnt%2==0){
              double tmp;
              tmp = p1->x_star;
              p1->x_star = p->c[i].x_star;
              p->c[i].x_star = tmp;
            }
          else
            {
              p1 = &p->c[i];
            }
        }
    }
}

/* apply mutation - random uniform mutation of the genes */
void apply_mutation(struct population *p)
{
  double prb = 0.0f;
  for(int i=0;i<p->size;++i){
      prb = rand()%1000/1000.0;
      if(prb < MUTATION_PROB){
          p->c[i].x_star = randomize(X_MIN, X_MAX);
        }
    }
}

/* print the state of the current evolution generation */
void report_state(struct population *p)
{
  printf("Generation: %d | Best fitness: %lf\n", p->gen, p->c[POPULATION_SIZE].fitness);
  fprintf(fp, "%d,%lf\n", p->gen, p->c[POPULATION_SIZE].fitness);
}

/* entry point */
int main(int argc, char* argv[]){
  fp = fopen("gas_demo_log.txt","w+");
  srand(time(NULL));
  printf("\n\nSimulation for GAs started...\n\n");
  /* the problem is to determine the maximum of f(x) = x*sin(10*pi*x) + 1
           in the interval x=[-1, 2]
                - analytically we have an infinite number of solutions for this equation
                - 3 types / forms of solutions
                      |	xi = (2*i-1)/20 + ei, i = 1,2....
                      | xi = 0
                      | xi = (2*i+1)/20 + ei, i = 1,2....
                - we split the input interval in 3 x 1000000 equal sized ranges
                - this means we have a 22-bit representation for chromosomes
        */

  struct population *p = (struct population*)calloc(1, sizeof(struct population));
  struct population *newp = (struct population*)calloc(1, sizeof(struct population));
  init_population(p, POPULATION_SIZE);
  init_population(newp, POPULATION_SIZE);
  evaluate_population(p);
  select_best(p);
  report_state(p);
  while(p->gen < MAX_GENERATIONS ){
      p->gen++;
      apply_selection(p, newp);
      apply_crossover(p);
      apply_mutation(p);
      report_state(p);
      evaluate_population(p);
      apply_elitism(p);
    }
  printf("\nEvolution is completed...\n\n");
  printf("\nBest chromosome: %lf | Best fitness: %lf\n\n", p->c[POPULATION_SIZE].x_star, p->c[POPULATION_SIZE].fitness);
  printf("\nSimulation ended.\n\n");
  fclose(fp);
  return EXIT_SUCCESS;
}
