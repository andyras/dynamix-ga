#include "initializeParams.hpp"

void init_gammas(GAGenome &c, Params * p) {
  // This function initializes the parameters which are to be changed in an
  // optimization over relaxation constants.

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;

  // variables to be changed
  double g1 = genome.gene(0);
  double g2 = genome.gene(1);
  double g1_c = genome.gene(2);

  // assign parameters from GA /////////////////////////////////////////////////

  p->gamma1 = g1;
  p->gamma2 = g2;
  p->gamma1_c = g1_c;

  initWavefunction(p);

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

  p->bulkGaussSigma = bulkGaussSigma;
  p->bulkGaussMu = bulkGaussMu;

  initWavefunction(p);

  return;
}

void init_wavepacketGammas(GAGenome &c, Params * p) {
  // This function initializes the parameters which are to be changed in an
  // optimization over relaxation constants and wavepacket parameters.

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;

  // variables to be changed
  double g1 = genome.gene(0);
  double g2 = genome.gene(1);
  double g1_c = genome.gene(2);
  double bulkGaussSigma = genome.gene(3);
  double bulkGaussMu = genome.gene(4);

  // assign parameters from GA /////////////////////////////////////////////////

  p->gamma1 = g1;
  p->gamma2 = g2;
  p->gamma1_c = g1_c;
  p->bulkGaussSigma = bulkGaussSigma;
  p->bulkGaussMu = bulkGaussMu;

  initWavefunction(p);

  return;
}

void init_torsion(GAGenome &c, Params * p) {
  // this function initializes variables relevant to having a torsional coupling

  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)c;

  p->torsionSin2V0 = genome.gene(0); // constant component of torsional coupling
  p->torsionSin2V1 = genome.gene(1); // time-dependent component of torsional coupling
  p->energies[p->Ib] = genome.gene(2);    // first bridge site energy
  p->energies[p->Ib+1] = genome.gene(2);  // second bridge site energy

  return;
}