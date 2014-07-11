#ifndef __OBJECTIVE__
#define __OBJECTIVE__

#include <sys/types.h>
#include <unistd.h>

#include "params.hpp"

double obj_tcpeak(Params * p);

double obj_Pcavg(Params * p);

double obj_Pcavg_after_peak(Params * p);

double obj_maxFinal(Params * p);

#endif
