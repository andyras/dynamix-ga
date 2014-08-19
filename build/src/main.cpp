#include <boost/mpi.hpp>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <dynamix.hpp>
#include <ga-mpi/ga.h>
#include <ga-mpi/std_stream.h>
#include <params.hpp>

#include "gaparams.hpp"
#include "initializer.hpp"
#include "mutator.hpp"
#include "objective.hpp"
#include "output.hpp"
#include "parser.hpp"

// #define DEBUG

int mpi_tasks, mpi_rank;

int main(int argc, char **argv)
{
  // MPI init
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

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

  if (world.rank() == 0) {
    assignGAParams("./ins/ga.in", &gp);
    assignParams("./ins/parameters.in", &(gp.p));
    initialize(&(gp.p));
    for (int ii = 1; ii < world.size(); ii++) {
      world.send(ii, 0, gp);
    }
  }
  else {
    world.recv(0, 0, gp);
  }

  void * userData = &gp;
#ifdef DEBUG
  std::cout << "[" << pid << ":" << mpi_rank <<  "] userData is " << userData << std::endl;
#endif

#ifdef DEBUG
  std::cout << "[" << pid << ":" << mpi_rank <<  "] Using " << gp.objectiveType << " objective function." << std::endl;
#endif

  unsigned int genomeLength = 0;
  if (gp.initializer.compare("g1g2g1_c") == 0) {
    genomeLength = 3;
    gp.lb.resize(genomeLength);
    gp.ub.resize(genomeLength);
    for (unsigned int ii = 0; ii < 3; ii++) {
      gp.lb[ii] = 0.001;
      gp.ub[ii] = 0.05;
    }
  }
  else if (gp.initializer.compare("wavepacket") == 0) {
    genomeLength = 2;
    gp.lb.resize(genomeLength);
    gp.ub.resize(genomeLength);
    for (unsigned int ii = 0; ii < 2; ii++) {
      gp.lb[ii] = 0.0001; // this needs to be ~> the interlevel spacing
      gp.ub[ii] = 0.01; // it would be nice to have these initialized by reading parameters.in
    }
  }
  else if (gp.initializer.compare("wavepacketGammas") == 0) {
    genomeLength = 5;
    gp.lb.resize(genomeLength);
    gp.ub.resize(genomeLength);
    for (unsigned int ii = 0; ii < 3; ii++) {
      gp.lb[ii] = 0.001;
      gp.ub[ii] = 0.05;
    }
    for (unsigned int ii = 3; ii < (3+2); ii++) {
      gp.lb[ii] = 0.0001;
      gp.ub[ii] = 0.01;
    }
  }
  else if (gp.initializer.compare("torsion") == 0) {
    genomeLength = 3;
    gp.lb.resize(genomeLength);
    gp.ub.resize(genomeLength);
    // torsionSin2V0
    gp.lb[0] = 0.00001;
    gp.ub[0] = 0.01;
    // torsionSin2V1
    gp.lb[1] = 0.00001;
    gp.ub[1] = 0.01;
    // Eb
    gp.lb[2] = gp.p.kBandEdge;
    gp.ub[2] = gp.p.kBandTop;
  }
  else {
    std::cout << "ERROR [" << __FUNCTION__ << "]: " << "variable set" <<
      gp.initializer << "not recognized." << std::endl;
    exit(-1);
  }

  GA1DArrayGenome<double> genome(genomeLength, getObjectiveType(&gp), userData);
  genome.initializer(getInitializer(&gp));

  genome.mutator(getMutator(&gp));

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

    // create new genome and evaluate objective (in order to create outputs)
    GA1DArrayGenome<double> bg(genomeLength, getObjectiveType(&gp), userData);
    // fill genome with best parameters
    for (int ii = 0; ii < bg.length(); ii++) {
      bg.gene(ii, bestGenomes[bestIdx*genomeLength + ii]);
    }
    // call objective function
    (* getObjectiveType(&gp))(bg);
  }

  MPI_Finalize();

  return 0;
}
