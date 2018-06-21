#include "property_type_constant.h"
#include "fix_property_global.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

PropertyTypeConstant::PropertyTypeConstant(LAMMPS *const lmp, char *const *const arg, const int narg, const int ncols, const int data_style, const double is_symmetric) :
    lmp_(lmp),
    nvalues(narg),
    nvalues_new_array(0),
    nrows_(-1),
    ncols_(ncols),
    values(NULL),
    values_recomputed(NULL),
    array(NULL),
    array_recomputed(NULL),
    data_style_(data_style)
{
    values = (double*) lmp_->memory->smalloc(nvalues*sizeof(double),"PropertyTypeConstant:values");
    values_recomputed = (double*) lmp_->memory->smalloc(nvalues*sizeof(double),"PropertyTypeConstant:values");

    for (int j = 0; j < nvalues; j++)
        values[j] = lmp_->force->numeric(FLERR,arg[j]);

    if(data_style_ == FIXPROPERTY_GLOBAL_MATRIX)
    {
        nrows_ = nvalues/ncols_;
        array = (double**)lmp_->memory->smalloc(nrows_*sizeof(double**),"PropertyTypeConstant:array");
        array_recomputed = (double**)lmp_->memory->smalloc(nrows_*sizeof(double**),"PropertyTypeConstant:array_recomputed");
        for(int i = 0; i < nrows_; i++)
            array[i] = &values[i*ncols_];
        for(int i = 0; i < nrows_; i++)
            array_recomputed[i] = &values_recomputed[i*ncols_];
    }

    // error check if matrix is symmetric (if required)
    if(is_symmetric)
    {
        if(nrows_ != ncols_)
            lmp_->error->all(FLERR, "per-atomtype property matrix must be symmetric, i.e. N atom types "
                                    "require you to define N columns and N rows with N*N total values");

        int sflag = true;
        for(int i = 0; i < nrows_; i++)
            for(int j = 0; j < ncols_; j++)
                if(array[i][j] != array[j][i])
                    sflag = false;

        if(!sflag)
            lmp_->error->all(FLERR, "per-atomtype property matrix must be symmetric");
    }
}

PropertyTypeConstant::~PropertyTypeConstant()
{
    lmp_->memory->sfree(values);
    lmp_->memory->sfree(values_recomputed);

    if(array)
        lmp_->memory->sfree(array);
    if(array_recomputed)
        lmp_->memory->sfree(array_recomputed);
}

double PropertyTypeConstant::compute_scalar(const double x)
{
    return values[0];
}

double PropertyTypeConstant::compute_vector(const double x, const int i)
{
    if (i >= nvalues)
        lmp_->error->one(FLERR, "Trying to access vector, but index out of bounds");
    return values[i];
}

double PropertyTypeConstant::compute_array(const double x, const int i, const int j)
{
    if (i >= nrows_)
        lmp_->error->one(FLERR, "Trying to access matrix, but row index out of bounds");
    if (j >= ncols_)
        lmp_->error->one(FLERR, "Trying to access matrix, but column index out of bounds");

    return array[i][j];
}

double PropertyTypeConstant::memory_usage()
{
    return nvalues*sizeof(double);
}

void PropertyTypeConstant::grow(const int len1, const int len2)
{
    if(data_style_ == FIXPROPERTY_GLOBAL_SCALAR)
        lmp_->error->one(FLERR, "Can not grow global property of type scalar");
    else if(data_style_ == FIXPROPERTY_GLOBAL_VECTOR && len1 > nvalues)
        lmp_->memory->grow(values,len1,"PropertyTypeConstant:values");
    else if(data_style_ == FIXPROPERTY_GLOBAL_MATRIX && len1*len2 > nvalues)
    {
        values = (double*) lmp_->memory->srealloc(values,len1*len2*sizeof(double),"PropertyTypeConstant:values");
        nrows_ = len1;
        ncols_ = len2;
        nvalues = len1*len2;
        array = (double**) lmp_->memory->srealloc(array, nrows_*sizeof(double**),"PropertyTypeConstant:array");
        for(int i = 0; i < nrows_; i++)
            array[i] = &values[i*ncols_];
    }
}

void PropertyTypeConstant::vector_modify(const int i, const double val)
{
    if (i >= nvalues)
        lmp_->error->one(FLERR, "Trying to access vector, but index out of bounds");
    values_recomputed[i] = val;
}

void PropertyTypeConstant::array_modify(const int i, const int j, const double val)
{
    if (i >= nrows_)
        lmp_->error->one(FLERR, "Trying to access matrix, but row index out of bounds");
    if (j >= ncols_)
        lmp_->error->one(FLERR, "Trying to access matrix, but column index out of bounds");

    array_recomputed[i][j] = val;
}

void PropertyTypeConstant::new_array(const int l1, const int l2)
{
    nrows_ = l1;
    ncols_ = l2;
    nvalues_new_array = l1*l2;

    lmp_->memory->create(array,nrows_,ncols_,"PropertyTypeConstant:array");
    lmp_->memory->create(array_recomputed,nrows_,ncols_,"PropertyTypeConstant:array_recomputed");
}
