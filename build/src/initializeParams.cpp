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

  p->torsionCouplingV0 = genome.gene(0); // constant component of torsional coupling
  p->torsionCouplingV1 = genome.gene(1); // time-dependent component of torsional coupling
  p->energies[p->Ib] = genome.gene(2);    // first bridge site energy
  p->energies[p->Ib+1] = genome.gene(2);  // second bridge site energy

  return;
}

void initParamsToChange(GAGenome &g, Params * p) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;
  GAParams * gp = (GAParams *) genome.userData();

  std::string name; // parameter name

  // flags to check if acceptor QC needs to be generated
  bool isNaChanged = false;
  bool isDaChanged = false;

  for (unsigned int ii = 0; ii < gp->paramsToChange.size(); ii++) {
    paramTuple ptc = gp->paramsToChange[ii];
    name = std::get<0>(ptc);
    if (name.compare("gamma1") == 0) {
      p->gamma1 = genome.gene(ii);
    }
    else if (name.compare("gamma2") == 0) {
      p->gamma2 = genome.gene(ii);
    }
    else if (name.compare("gamma1_c") == 0) {
      p->gamma1_c = genome.gene(ii);
    }
    else if (name.compare("bulkGaussSigma") == 0) {
      p->bulkGaussSigma = genome.gene(ii);
    }
    else if (name.compare("bulkGaussMu") == 0) {
      p->bulkGaussMu = genome.gene(ii);
    }
    else if (name.compare("torsionCouplingV0") == 0) {
      p->torsionCouplingV0 = genome.gene(ii);
    }
    else if (name.compare("torsionCouplingV1") == 0) {
      p->torsionCouplingV1 = genome.gene(ii);
    }
    else if (name.compare("torsionCouplingOmega") == 0) {
      p->torsionCouplingOmega = genome.gene(ii);
    }
    else if (name.compare("torsionCouplingPhi") == 0) {
      p->torsionCouplingPhi = genome.gene(ii);
    }
    else if (name.compare("Nb") == 0) {
      p->Nb = genome.gene(ii);
      p->b_energies.resize(p->Nb, 0.0);
    }
    else if (name.compare("Na") == 0) {
      p->Nc = genome.gene(ii);
      p->c_energies.resize(p->Nc, 0.0);
    }
    else if (name.compare("Da") == 0) {
      gp->Da = genome.gene(ii);
    }
    else if (name.compare("Ea") == 0) {
      gp->Ea = genome.gene(ii);
    }
    else {
      std::cerr << "ERROR: unrecognized parameter name '" << name << "'" << std::endl;
      _exit(-1);
    }

    // check if Na was set
    if (isNaChanged || isDaChanged) {
      for (int ii = 0; ii < p->Nc; ii++) {
        p->c_energies[ii] = gp->Ea + ii*gp->Da;
      }
    }
    // create acceptor QC according to GAParams
  }

  return;
}