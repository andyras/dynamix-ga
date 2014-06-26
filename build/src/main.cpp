#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <ga-mpi/ga.h>
#include <ga-mpi/std_stream.h>
#include "mpi.h"
#include <sys/wait.h>
#include <string>

#include "dynamix.hpp"
#include "propagate.hpp"

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
  int popsize  = 2; // Population
  int ngen     = 2; // Generations
  float pmut   = 0.03;
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

float objective(GAGenome &c)
{
  ///GABin2DecGenome &genome = (GABin2DecGenome &)c;
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  float x, y, error;

  x = genome.gene(0);
  y = genome.gene(1);

  // Function with local minima. The lowest is located at (5/2*PI, 5/2*PI)
  error = ((1.-sin(x)*sin(y))+sqrt((x-M_PI*2.5)*(x-M_PI*2.5)+(y-M_PI*2.5)*(y-M_PI*2.5))/10.0)/2.5;

  return error;
}

float dynamixObjective(GAGenome &c) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;
  // get MPI rank and size
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // variables and output
  double g1 = genome.gene(0);
  double g2 = genome.gene(1);
  double g1_c = genome.gene(2);
  float output = 1.0;

  // array of strings for commands + arguments
  char * args [10];
  for (int ii = 0; ii < 10; ii++) {
    args[ii] = new char [1000];
  }

  std::string arg = "/extra/scratch/foo";

  // names for original directory and job directory
  // TODO pass as cmd line params
  std::string tmpDir ("/extra/scratch/dynamix");
  // job directory name
  char jobFmt [] = "%s%.12e_%.12e_%.12e_%d_%d";
  char jobDir [1000];
  sprintf(jobDir, jobFmt, tmpDir.c_str(), g1, g2, g1_c, rank, size);
  std::cout << jobDir << std::endl;
  std::string dynDir ("/home/andyras/git/dynamix-ga/");
  std::string dynInsDir (dynDir + "ins/");
  std::string dynamix (dynDir + "bin/dynamix-ga");
  std::string changeParam (dynDir + "tools/changeParam.py");
  char sciFmt [] = "%.12e";
  char g1Str [100];
  char g2Str [100];
  char g1_cStr [100];
  sprintf(g1Str, sciFmt, g1);
  sprintf(g2Str, sciFmt, g2);
  sprintf(g1_cStr, sciFmt, g1_c);

  std::string jobInsDir (jobDir + std::string("/ins/"));
  std::string jobOutsDir (jobDir + std::string("/outs/"));

  // construct arguments
  for (int ii = 0; ii < 10; ii++) {
    delete [] args[ii];
    args[ii] = new char [1000];
  }

  // for debugging
  jobInsDir = "./ins/";
  jobOutsDir = "/tmp/";

  strncpy(args[0], dynamix.c_str(), dynamix.length() + 1);
  strncpy(args[1], "-i", strlen("-i") + 1);
  strncpy(args[2], jobInsDir.c_str(), jobInsDir.length() + 1);
  strncpy(args[3], "-o", strlen("-o") + 1);
  strncpy(args[4], jobOutsDir.c_str(), jobOutsDir.length() + 1);
  args[5] = NULL;

  fprintf(stdout, "COMMAND:");
  for (int ii = 0; ii < 10; ii++) fprintf(stdout, " %s", args[ii]);
  fprintf(stdout, "\n");

  dynamixMain(5, args);

  // ---- check for success ---- //
  //
  // ---- read in outputs ---- //
  //
  // ---- calculate objective ---- //
  //
  // ---- remove job directory ---- //

  // clean up
  for (int ii = 0; ii < 10; ii++) {
    delete [] args[ii];
  }

  return output;
}

void Initializer(GAGenome &g) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;

  // there are two genes
  genome.gene(0, GARandomFloat(0.0,5*M_PI));
  genome.gene(1, GARandomFloat(0.0,5*M_PI));
  genome.gene(2, GARandomFloat(0.0,5*M_PI));

  return;
}

int dynamixMain (int argc, char * argv[]) {
  // Struct of parameters //////////////////////////////////////////////////////

  Params p;

  // process command line flags ////////////////////////////////////////////////

  opterr = 0;
  int c;
  std::string insDir;
  /* process command line options */
  while ((c = getopt(argc, argv, "i:o:")) != -1) {
    switch (c) {
      case 'i':
        // check that it ends in a slash
        insDir = optarg;
        if (strcmp(&(insDir.at(insDir.length() - 1)), "/")) {
          std::cerr << "ERROR: option -i requires argument ("
                    << insDir << ") to have a trailing slash (/)." << std::endl;
          return 1;
        }
        else {
          // ---- assign input files ---- //
          p.inputFile = insDir + "parameters.in";
          p.cEnergiesInput = insDir + "c_energies.in";
          p.bEnergiesInput = insDir + "b_energies.in";
          p.VNoBridgeInput = insDir + "Vnobridge.in";
          p.VBridgeInput = insDir + "Vbridge.in";
        }
        break;
      case 'o':
        p.outputDir = optarg;
        break;
      case '?':
        if (optopt == 'i') {
          fprintf(stderr, "Option -%c requires a directory argument.\n", optopt);
        }
        else if (isprint(optopt)) {
          fprintf(stderr, "Unknown option -%c.\n", optopt);
        }
        else {
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        }
        return 1;
      default:
        continue;
    }
  }

  // assign parameters from input file /////////////////////////////////////////

  assignParams(p.inputFile.c_str(), &p);

  initialize(&p);

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p.nproc);
  mkl_set_num_threads(p.nproc);

  // Make plot files ///////////////////////////////////////////////////////////

  makePlots(&p);

  // only do propagation if not just making plots //////////////////////////////

  if (! p.justPlots) {
    propagate(&p);
  }

  std::cout << "whoo" << std::endl;

  return 0;
}