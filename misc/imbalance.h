#ifndef LMP_IMBALANCE_H
#define LMP_IMBALANCE_H

#include <stdio.h>
#include "pointers.h"

namespace LAMMPS_NS {

class Imbalance : protected Pointers {
 public:
  Imbalance(class LAMMPS *):Pointers(lmp) {}
  virtual ~Imbalance() {};
  virtual int options(int, char **) {}
  virtual void init(int) {} 
  virtual void compute(double *) {}
  virtual void info(FILE *) {}
};

typedef Imbalance ImbalanceGroup;
typedef Imbalance ImbalanceVar;
typedef Imbalance ImbalanceTime;
typedef Imbalance ImbalanceNeigh;
typedef Imbalance ImbalanceStore;

}

#endif
