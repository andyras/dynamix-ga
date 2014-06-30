#include <dynamix.hpp>

#include "initializer.hpp"

#define DEBUG

void init_gammas(GAGenome &c, Params * p) {
  // This function initializes the parameters which are to be changed in an
  // optimization over relaxation constants.

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;

  // variables to be changed
  double g1 = genome.gene(0);
  double g2 = genome.gene(1);
  double g1_c = genome.gene(2);

  // assign parameters from GA /////////////////////////////////////////////////

#ifdef DEBUG
  std::cout << "Value of gamma1 is: " << g1 << std::endl;
  std::cout << "Value of gamma1 is: " << g2 << std::endl;
  std::cout << "Value of gamma1 is: " << g1_c << std::endl;
#endif

  p->gamma1 = g1;
  p->gamma2 = g2;
  p->gamma1_c = g1_c;

  initialize(p);

  return;
}

void init_wavepacket(GAGenome &c, Params * p) {
  // This function initializes the parameters which are to be changed in an
  // optimization over relaxation constants.

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;

  // variables to be changed
  double bulkGaussSigma = genome.gene(0);
  double bulkGaussMu = genome.gene(1);

  // assign parameters from GA /////////////////////////////////////////////////

#ifdef DEBUG
  std::cout << "Value of bulkGaussSigma is: " << bulkGaussSigma << std::endl;
  std::cout << "Value of bulkGaussMu is: " << bulkGaussMu << std::endl;
#endif

  p->bulkGaussSigma = bulkGaussSigma;
  p->bulkGaussMu = bulkGaussMu;

  initialize(p);

  return;
}