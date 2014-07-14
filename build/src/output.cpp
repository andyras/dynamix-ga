#include "output.hpp"

void print1DGenes(GAGenome & g) {
  /* This function prints the genes contained in a 1D array genome. */
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int pid = getpid();

  // cast as 1D array
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;

  // print first guy
  std::cout << "[" << pid << ":" << rank << "] Gene 0: " << genome.gene(0);
  for (int ii = 1; ii < genome.length(); ii++) {
    std::cout << " Gene " << ii << ": " << genome.gene(ii);
  }
  std::cout << std::endl;

  return;
}