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
#include <propagate.hpp>

#include "objective.hpp"
#include "initializer.hpp"

#define DEBUG

// This are the declaration of the objective functions which are defined later.
float objective(GAGenome &);
float dynamixObjective(GAGenome &);
int dynamixMain (int argc, char * argv[]);

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
  int popsize  = 40; // Population
  int ngen     = 200; // Generations
  float pmut   = 0.05;
  float pcross = 0.65;

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
  //GA1DArrayGenome<double> genome(2, objective);
  GA1DArrayGenome<double> genome(3, dynamixObjective);
  // define own initializer, can do the same for mutator and comparator
  genome.initializer(::Initializer);

  omp_set_num_threads(2);
  mkl_set_num_threads(2);

  // Now create the GA using the genome and run it. We'll use sigma truncation
  // scaling so that we can handle negative objective scores.
  GASimpleGA ga(genome); // TODO change to steady-state
  GALinearScaling scaling;
  ga.minimize();    // by default we want to minimize the objective
  ga.populationSize(popsize);
  ga.nGenerations(ngen);
  ga.pMutation(pmut);
  ga.pCrossover(pcross);
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
    genome = ga.statistics().bestIndividual();
    printf("GA result:\n");
    printf("x = %f, y = %f\n",
  genome.gene(0), genome.gene(1));
  }

  MPI_Finalize();

  return 0;
}

float dynamixObjective(GAGenome &c) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  // get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

#ifdef DEBUG
  std::cout << "rank: " << rank << " of " << size << "\n";
#endif

  // output value
  float output = 1.0;

  // Struct of parameters //////////////////////////////////////////////////////

  Params p;
  // TODO XXX FUGLINESS HARD-CODING
  std::string jobPrefix = "/home/andyras/git/dynamix-ga/";
  p.inputFile = jobPrefix + "ins/parameters.in";
  p.cEnergiesInput = jobPrefix + "ins/c_energies.in";
  p.bEnergiesInput = jobPrefix + "ins/b_energies.in";
  p.VNoBridgeInput = jobPrefix + "ins/Vnobridge.in";
  p.VBridgeInput = jobPrefix + "ins/Vbridge.in";

  // assign parameters from input file /////////////////////////////////////////

  assignParams(p.inputFile.c_str(), &p);

  initialize(&p);

  // assign GA parameters
  // make sure this function matches the objective you are using
  init_wavepacket(c, &p);

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p.nproc);
  mkl_set_num_threads(p.nproc);

  // Make plot files ///////////////////////////////////////////////////////////

  makePlots(&p);

  // only do propagation if not just making plots //////////////////////////////

  if (! p.justPlots) {
    propagate(&p);
  }

  // calculate value of objective function /////////////////////////////////////
#ifdef DEBUG
  std::cout << "Calculating value of objective function..." << std::endl;
#endif
  output = obj_tcpeak(&p);

  std::cout << "whoo" << std::endl;

  return output;
}

void Initializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;

  genome.gene(0, GARandomFloat(0.0,0.01));
  genome.gene(1, GARandomFloat(0.0,0.01));
  genome.gene(2, GARandomFloat(0.0,0.01));

  return;
}
