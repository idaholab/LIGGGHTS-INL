#ifndef LAMMPS_RCB_H
#define LAMMPS_RCB_H

#include <mpi.h>
#include "pointers.h"
#include "lammps.h"

namespace LAMMPS_NS {

class RCB : protected Pointers {
 public:

  RCB(class LAMMPS *lmp):Pointers(lmp) {}
  virtual   ~RCB() {}
  void compute(int, int, double **, double *, double *, double *){}
  void compute_old(int, int, double **, double *, double *, double *){}
  void invert(int sortflag = 0){}
  bigint memory_usage(){return 0;}

  int noriginal;              
  int nfinal;                 
  int nkeep;                  
                              
  int *recvproc;              
  int *recvindex;             
                              
  double *lo,*hi;             
  double cut;                 
  int **sendproc;              
  int **sendindex;             
  int cutdim;                 
};

}

#endif
