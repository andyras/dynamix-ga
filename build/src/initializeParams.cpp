#include "initializeParams.hpp"

// #define DEBUG

void initParamsToChange(GAGenome &g, Params * p) {
  GA1DArrayGenome<double> &genome = (GA1DArrayGenome<double> &)g;
  GAParams * gp = (GAParams *) genome.userData();

  std::string name; // parameter name

  // flags to check if acceptor QC needs to be generated
  bool isNaChanged = false;
  bool isDaChanged = false;

#ifdef DEBUG
  std::cout << "Changing " << gp->paramsToChange.size() << " parameters." << std::endl;
#endif

  for (unsigned int ii = 0; ii < gp->paramsToChange.size(); ii++) {
    paramTuple ptc = gp->paramsToChange[ii];
    name = std::get<0>(ptc);

#ifdef DEBUG
    std::cout << "Parameter: " << std::get<0>(ptc) <<
      " lower bound: " << std::get<1>(ptc) <<
      " upper bound: " << std::get<2>(ptc) << std::endl;
#endif

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
      p->Nb = static_cast<int>(genome.gene(ii));
      p->b_energies.resize(p->Nb, 0.0);
    }
    else if (name.compare("Na") == 0) {
      p->Nc = static_cast<int>(genome.gene(ii));
      p->c_energies.resize(p->Nc, 0.0);
      isNaChanged = true;
    }
    else if (name.compare("Da") == 0) {
      gp->Da = genome.gene(ii);
      isDaChanged = true;
    }
    else if (name.compare("Ea") == 0) {
      gp->Ea = genome.gene(ii);
    }
    else if (name.compare("Eb") == 0) {
      for (unsigned int jj = 0; jj < p->b_energies.size(); jj++) {
        p->b_energies[jj] = genome.gene(ii);
      }
    }
    else {
      std::cerr << "ERROR: unrecognized parameter name '" << name << "'" << std::endl;
      _exit(-1);
    }

    // check if Na or Da were set
    if (isNaChanged || isDaChanged) {
      for (int ii = 0; ii < p->Nc; ii++) {
        p->c_energies[ii] = gp->Ea + ii*gp->Da;
      }
    }
    // create acceptor QC according to GAParams
  }

  return;
}