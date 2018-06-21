#ifndef PROPERTY_TYPE_CONSTANT_H
#define PROPERTY_TYPE_CONSTANT_H

#include "property_type_base.h"
#include "fix_property_global.h"
#include "lammps.h"

namespace LAMMPS_NS
{

class PropertyTypeConstant : public PropertyTypeBase
{
public:
    PropertyTypeConstant(LAMMPS *const lmp, char *const *const arg, const int narg, const int ncols, const int data_style, const double is_symmetric);
    ~PropertyTypeConstant();

    double compute_scalar(const double x);
    double compute_vector(const double x, const int i);
    double compute_array(const double x, const int i, const int j);
    double memory_usage();
    void grow(const int len1, const int len2);
    double* get_values() const { return values; }
    double* get_values_modified() { return values_recomputed; }
    double const* const* get_array() { return array; }
    double const* const* get_array_modified() { return array_recomputed; }
    void vector_modify(const int i, const double val);
    void array_modify(const int i, const int j, const double val);
    void new_array(const int l1, const int l2);

private:
    LAMMPS *lmp_;
    int nvalues;
    int nvalues_new_array;
    int nrows_, ncols_;
    double *values;            // original values to be stored in this fix
    double *values_recomputed; // values to be stored in this fix, recomputed by eg another fix
    double **array;
    double **array_recomputed;
    int data_style_;
};

}

#endif
