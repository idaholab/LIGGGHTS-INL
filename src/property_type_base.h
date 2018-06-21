#ifndef PROPERTY_TYPE_BASE_H
#define PROPERTY_TYPE_BASE_H

namespace LAMMPS_NS
{

class PropertyTypeBase
{
public:
    PropertyTypeBase() {}
    virtual ~PropertyTypeBase() {}

    virtual double compute_scalar(const double x) = 0;
    virtual double compute_vector(const double x, const int i) = 0;
    virtual double compute_array(const double x, const int i, const int j) = 0;
    virtual double memory_usage() = 0;

    virtual void grow(const int len1, const int len2) {}
};

}

#endif
