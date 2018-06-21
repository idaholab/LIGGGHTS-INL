#include "property_type_base.h" 
 #include "lammps.h" 
 namespace LAMMPS_NS{ 
 class PropertyTypeLookup : public PropertyTypeBase 
 { 
 public: 
 PropertyTypeLookup(LAMMPS * const , const char * const arg[], const int, const int, const int, const bool, const bool) {} 
  double compute_scalar(const double)  { return 0;} 
     double compute_vector(const double, const int i)  { return 0;} 
     double compute_array(const double, const int, const int)  { return 0;} 
     double memory_usage() { return 0;} 
 }; 
 }
