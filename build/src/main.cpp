#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <ga-mpi/ga.h>
#include <ga-mpi/std_stream.h>
#include <mpi.h>
#include <sys/wait.h>
#include <string>

#include <dynamix.hpp>

#include "objective.hpp"
#include "initializer.hpp"
#include "output.hpp"
#include "parser.hpp"
#include "gaparams.hpp"

// #define DEBUG

// declare Initializer also
void Initializer(GAGenome &);

int mpi_tasks, mpi_rank;

int main(int argc, char **argv)
{
  // MPI init
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // See if we've been given a seed to use (for testing purposes).  When you
  // specify a random seed, the evolution will be exactly the same each time
  // you use that seed number
  unsigned int seed = 0;
  for(int i=1 ; i<argc ; i++)
    if(strcmp(argv[i++],"seed") == 0)
      seed = atoi(argv[i]);

  // Declare variables for the GA parameters and set them to some default values.
  int popsize  = 4; // Population
  int ngen     = 200; // Generations
  float pmut   = 0.25;
  float pcross = 0.65;
  float pconv = 1.01; // convergence

  GAParams gp;

  // popsize / mpi_tasks must be an integer
  popsize = mpi_tasks * int((double)popsize/(double)mpi_tasks+0.999);

  // Create the phenotype for two variables.  The number of bits you can use to
  // represent any number is limited by the type of computer you are using.
  // For this case we use 10 bits for each var, ranging the square domain [0,5*PI]x[0,5*PI]
  ///GABin2DecPhenotype map;
  ///GABin2DecPhenotype map;
  ///map.add(10, 0.0, 5.0 * M_PI);
  ///map.add(10, 0.0, 5.0 * M_PI);

  // Create the template genome using the phenotype map we just made.
  ///GABin2DecGenome genome(map, objective);
  GA1DArrayGenome<double> genome(3, dynamixObjective);
  if (gp.objectiveType.compare("single") == 0) {
    GA1DArrayGenome<double> genome(3, dynamixObjective);
  }
  else if (gp.objectiveType.compare("double") == 0) {
    GA1DArrayGenome<double> genome(3, dualObjective);
  }
  else {
    std::cout << "WARNING [" << __FUNCTION__ << "]: " << "objective type" <<
      gp.objectiveType << "not recognized." << std::endl;
    exit(-1);
  }

  // define own initializer, can do the same for mutator and comparator
  if (gp.variables.compare("g1g2g1_c") == 0) {
    genome.initializer(::gammasInitializer);
  }
  else if (gp.variables.compare("wavepacket") == 0) {
    genome.initializer(::wavepacketInitializer);
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      gp.variables << "not recognized." << std::endl;
    exit(-1);
  }

  omp_set_num_threads(1);
  mkl_set_num_threads(1);

  // Now create the GA using the genome and run it. We'll use sigma truncation
  // scaling so that we can handle negative objective scores.
  GASimpleGA ga(genome); // TODO change to steady-state
  GALinearScaling scaling;
  ga.minimize();    // by default we want to minimize the objective
  ga.populationSize(popsize);
  ga.nGenerations(ngen);
  ga.pMutation(pmut);
  ga.pCrossover(pcross);
  ga.pConvergence(pconv);
  ga.terminator(GAGeneticAlgorithm::TerminateUponConvergence);
  ga.scaling(scaling);
  if(mpi_rank == 0)
    ga.scoreFilename("evolution.txt");
  else
    ga.scoreFilename("/dev/null");
  ga.scoreFrequency(1);
  ga.flushFrequency(1);
  ga.selectScores(GAStatistics::AllScores);
  // Pass MPI data to the GA class
  ga.mpi_rank(mpi_rank);
  ga.mpi_tasks(mpi_tasks);
  ga.evolve(seed);

  omp_set_num_threads(1);
  mkl_set_num_threads(1);

  // Dump the GA results to file
  if(mpi_rank == 0)
  {
    GA1DArrayGenome<double> &bestGenome = (GA1DArrayGenome<double> &)ga.population().best();

    std::cout << "Best: Gene 0: " << bestGenome.score();
    for (int ii = 0; ii < bestGenome.length(); ii++) {
      std::cout << " Gene " << ii << ": " << bestGenome.gene(ii);
    }
    std::cout << " Score: " << ga.population().best().score();
    std::cout << std::endl;
    }

  MPI_Finalize();

  return 0;
}