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

int mpi_tasks, mpi_rank;

int main(int argc, char **argv)
{
  // MPI init
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int pid = getpid();

  // See if we've been given a seed to use (for testing purposes).  When you
  // specify a random seed, the evolution will be exactly the same each time
  // you use that seed number
  unsigned int seed = 0;
  for(int i=1 ; i<argc ; i++)
    if(strcmp(argv[i++],"seed") == 0)
      seed = atoi(argv[i]);

  // Declare variables for the GA parameters and set them to some default values.
  int ngen = 200; // Generations

  GAParams gp;

  // each thread needs to read this file separately.
  int pRank = 0;
  while (pRank < mpi_tasks) {
    if (mpi_rank == pRank) {
      assignGAParams("./ins/ga.in", &gp);
    }
    pRank++;
    MPI_Barrier(MPI_COMM_WORLD);
  }

  void * userData = &gp;
#ifdef DEBUG
  std::cout << "[" << pid << ":" << mpi_rank <<  "] userData is " << userData << std::endl;
#endif

#ifdef DEBUG
  std::cout << "[" << pid << ":" << mpi_rank <<  "] Using " << gp.objectiveType << " objective function." << std::endl;
#endif
  if (gp.objectiveType.compare("single") == 0) {
    gp.objectiveFn = singleObjective;
  }
  else if (gp.objectiveType.compare("double") == 0) {
    gp.objectiveFn = doubleObjective;
  }
  else {
    std::cout << "WARNING [" << __FUNCTION__ << "]: " << "objective type" <<
      gp.objectiveType << "not recognized." << std::endl;
    exit(-1);
  }

  unsigned int genomeLength = 0;
  if (gp.initializer.compare("g1g2g1_c") == 0) {
    gp.initializerFn = ::gammasInitializer;
    genomeLength = 3;
  }
  else if (gp.initializer.compare("wavepacket") == 0) {
    gp.initializerFn = ::wavepacketInitializer;
    genomeLength = 2;
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      gp.initializer << "not recognized." << std::endl;
    exit(-1);
  }

  GA1DArrayGenome<double> genome(genomeLength, gp.objectiveFn, userData);
  genome.initializer(gp.initializerFn);

  // initialize structures for best genomes/scores
  double initVal = 0.0;
  if (gp.minmax.compare("min") == 0) {
    initVal = INFINITY;
  }
  else if (gp.minmax.compare("max") == 0) {
    initVal = -INFINITY;
  }

  gp.bestScore = initVal;
  gp.bestGenome.resize(genomeLength);
  // Now create the GA using the genome and run it. We'll use sigma truncation
  // scaling so that we can handle negative objective scores.
  GASimpleGA ga(genome); // TODO change to steady-state
  ga.userData(userData);
  GALinearScaling scaling;
  if (gp.minmax.compare("min") == 0) {
    ga.minimize();
    ga.pConvergence(1.0 + gp.convergence);
  }
  else if (gp.minmax.compare("max") == 0) {
    ga.maximize();
    ga.pConvergence(1.0 - gp.convergence);
  }
  else {
    std::cerr << "ERROR: unrecognized minmax: << " << gp.minmax << std::endl;
    exit(0);
  }
  ga.populationSize(gp.popsize);
  ga.nGenerations(ngen);
  ga.pMutation(gp.pMut);
  ga.pCrossover(gp.pCross);
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

  std::vector<double> bestScores;
  std::vector<double> bestGenomes;

  if (mpi_rank == 0) {
    bestScores.resize(mpi_tasks, 0.0);
    bestGenomes.resize(mpi_tasks*genomeLength, 0.0);
  }

  // wait for everyone to finish ///////////////////////////////////////////////
  MPI_Barrier(MPI_COMM_WORLD);

  // gather best scores ////////////////////////////////////////////////////////
  MPI_Gather(&gp.bestScore, 1, MPI_DOUBLE,
             &(bestScores[0]), 1, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  // gather best genomes ///////////////////////////////////////////////////////
  MPI_Gather(&(gp.bestGenome[0]), genomeLength, MPI_DOUBLE,
             &(bestGenomes[0]), genomeLength, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  // Dump the GA results to file ///////////////////////////////////////////////
  if(mpi_rank == 0)
  {
    // find rank of best score
    auto bestScore = std::min_element(bestScores.begin(), bestScores.end());
    if (gp.minmax.compare("max") == 0) {
      bestScore = std::max_element(bestScores.begin(), bestScores.end());
    }
    unsigned int bestIdx = std::distance(bestScores.begin(), bestScore);

    std::cout << "[" << pid << ":" << mpi_rank << "] Best score: " <<
      bestScores[bestIdx] << " Genome:";
    for (unsigned int ii = 0; ii < genomeLength; ii++) {
      std::cout << " " << bestGenomes[bestIdx*genomeLength + ii];
    }
    std::cout << std::endl;

    std::cout << "[" << pid << ":" << mpi_rank << "] Scores:";
    for (unsigned int ii = 0; ii < bestScores.size(); ii++) {
      std::cout << " " << bestScores[ii];
    }
    std::cout << std::endl;
  }

  MPI_Finalize();

  return 0;
}
