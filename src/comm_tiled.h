#ifndef LMP_COMM_TILED_H
#define LMP_COMM_TILED_H

#include "comm.h"

namespace LAMMPS_NS {

class CommTiled : public Comm {
 public:
  CommTiled(class LAMMPS *lmp) :Comm(lmp) {}
  CommTiled(class LAMMPS *lmp, class Comm *cm): Comm(lmp){}
  CommTiled(Comm *main_comm, ParallelBase *) : Comm(*main_comm) {}
  virtual ~CommTiled(){}

  void init(){}
  void setup(){}                        
  virtual void forward_comm(int dummy = 0){}   
  virtual void reverse_comm(){}                
  virtual void exchange(){}                    
  virtual void borders(){}                     

  virtual void forward_comm_pair(class Pair *){} 
  virtual void reverse_comm_pair(class Pair *){} 
  virtual void forward_comm_fix(class Fix *){}   
  virtual void reverse_comm_fix(class Fix *){}   
  virtual void forward_comm_variable_fix(class Fix *){}

  virtual void forward_comm_compute(class Compute *){} 
  virtual void reverse_comm_compute(class Compute *){} 
  virtual void forward_comm_dump(class Dump *){}   
  virtual void reverse_comm_dump(class Dump *){}   

  virtual void forward_comm_array(int, double **){} 
  void coord2proc_setup(){}
  int coord2proc(double *, int &, int &, int &){return 0;}
  bigint memory_usage(){return 0;}

  void evaluate_box_drop(int dim, double *lo, double *hi, int &indexme)  {  }
  void evaluate_box_other(int dim, int dir, int proc, double *lo, double *hi)  {  }
  int evaluate_box_touch(int proc, int dim, int dir)  { return 0; }
  int evaluate_point_drop(int dim, double *x)  { return 0; }
};

}

#endif
