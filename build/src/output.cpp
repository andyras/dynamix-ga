#include "output.hpp"

void print1DGenes(GAGenome & g) {
  /* This function prints the genes contained in a 1D array genome. */

  // cast as 1D array
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;

  // print first guy
  std::cout << "Gene 0: " << genome.gene(0);
  for (int ii = 1; ii < genome.length(); ii++) {
    std::cout << " Gene " << ii << ": " << genome.gene(ii);
  }
  std::cout << std::endl;

  return;
}