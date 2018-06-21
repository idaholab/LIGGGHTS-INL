/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include "fix_property_global.h"
#include "property_type_constant.h"
#include "property_type_lookup.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

FixPropertyGlobal::FixPropertyGlobal(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    variablename(NULL),
    data_style(FIXPROPERTY_GLOBAL_SCALAR),
    data_type_(FIXPROPERTY_GLOBAL_TYPE_CONSTANT),
    is_symmetric(false),
    is_atomtype_bound(false),
    filename(NULL),
    grpname(NULL),
    me(-1),
    property_(NULL)
{
    //Check args
    if (narg < 6)
        error->all(FLERR,"Illegal fix property/global command, not enough arguments");

    if(strcmp(style,"custom_property/global") == 0)
    {
        int len = strlen("property/global") + 1;
        delete []style;
        style = new char[len];
        strcpy(style,"property/global");
    }

    //Read args
    int n = strlen(arg[3]) + 1;
    variablename = new char[n];
    strcpy(variablename,arg[3]);

    if (strcmp(arg[4],"scalar") == 0)
        data_style = FIXPROPERTY_GLOBAL_SCALAR;
    else if (strcmp(arg[4],"vector") == 0)
        data_style = FIXPROPERTY_GLOBAL_VECTOR;
    else if (strcmp(arg[4],"peratomtype") == 0 || strcmp(arg[4],"atomtype") == 0)
    {
        data_style = FIXPROPERTY_GLOBAL_VECTOR;
        is_atomtype_bound = true;
    }
    else if (strcmp(arg[4],"matrix") == 0)
        data_style = FIXPROPERTY_GLOBAL_MATRIX;
    else if (strcmp(arg[4],"peratomtypepair") == 0 || strcmp(arg[4],"atomtypepair") == 0)
    {
        data_style = FIXPROPERTY_GLOBAL_MATRIX;
        is_symmetric = true;
        is_atomtype_bound = true;
    }
    else
        error->fix_error(FLERR,this,"Unknown style. Valid styles are scalar, vector, atomtype/peratomtype, matrix, or atomtypepair/peratomtypepair");

    int darg = 0;
    if (data_style == FIXPROPERTY_GLOBAL_MATRIX)
    {
        size_array_cols = force->inumeric(FLERR,arg[5]);
        darg += 1;
    }

    if (strcmp(arg[5+darg], "type") == 0)
    {
        if (narg < 7+darg)
            error->all(FLERR,"Illegal fix property/global command, not enough arguments");
        if (strcmp(arg[6+darg], "constant") == 0)
            data_type_ = FIXPROPERTY_GLOBAL_TYPE_CONSTANT;
        else if (strcmp(arg[6+darg], "lookup") == 0)
            data_type_ = FIXPROPERTY_GLOBAL_TYPE_LOOKUP;
        else
            error->all(FLERR, "Illegal property/global type");
        darg += 2;
    }

    bool warn_out_of_bounds(false);
    if(strcmp(arg[5+darg], "warn_out_of_bounds") == 0)
    {
        if(data_type_ != FIXPROPERTY_GLOBAL_TYPE_LOOKUP)
            error->fix_error(FLERR,this,"warn_out_of_bounds only valid for property/global type lookup");

        if(strcmp(arg[6+darg],"on") == 0)
            warn_out_of_bounds = true;
        else if(strcmp(arg[6+darg],"off") == 0)
            warn_out_of_bounds = false;
        else
            error->all(FLERR,"illegal value for warn_out_of_bounds");
        darg += 2;
    }

    //assign values
    nvalues = narg - 5 - darg;

    if (data_style == FIXPROPERTY_GLOBAL_MATRIX)
    {
        if (fmod(static_cast<double>(nvalues),size_array_cols) != 0.)
            error->fix_error(FLERR,this,"the number of default values must be a multiple of nCols.");
        size_array_rows = static_cast<int>(static_cast<double>(nvalues)/size_array_cols);
    }
    
    if (data_type_ == FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
    {
        property_ = static_cast<PropertyTypeBase*>(new PropertyTypeConstant(lmp, &arg[5+darg], nvalues, size_array_cols, data_style, is_symmetric));

        if (data_style == FIXPROPERTY_GLOBAL_SCALAR)
            scalar_flag = 1;
        else if (data_style==FIXPROPERTY_GLOBAL_VECTOR)
        {
            vector_flag = 1;
            size_vector = nvalues;
        }
        else if (data_style == FIXPROPERTY_GLOBAL_MATRIX)
            array_flag = 1;
    }
    else if (data_type_ == FIXPROPERTY_GLOBAL_TYPE_LOOKUP)
        property_ = static_cast<PropertyTypeBase*>(new PropertyTypeLookup(lmp, &arg[5+darg], nvalues, size_array_cols, data_style, is_symmetric, warn_out_of_bounds));

    extvector=0; 

    filename = 0;
    grpname = 0;

    //check if there is already a fix that tries to register a property with the same name
    
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((modify->fix[ifix]) && (strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyGlobal*)(modify->fix[ifix]))->variablename,variablename)==0) )
            error->fix_error(FLERR,this,"There is already a fix that registers a variable of the same name");
}

/* ---------------------------------------------------------------------- */

FixPropertyGlobal::~FixPropertyGlobal()
{
    // delete locally stored arrays
    delete[] variablename;

    if(filename)
        delete[] filename;
    if(grpname)
        delete[] grpname;

    delete property_;
}

/* ---------------------------------------------------------------------- */

void FixPropertyGlobal::pre_delete(bool unfixflag)
{
    if(filename)
        write();
}

/* ---------------------------------------------------------------------- */

Fix* FixPropertyGlobal::check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag)
{
    char errmsg[400];

    if(strcmp(varname,variablename) == 0)
    {
        if(strcmp(svmstyle,"scalar") == 0)
            len1 = 1;

        // check variable style
        if(
            (strcmp(svmstyle,"scalar") == 0 && data_style != FIXPROPERTY_GLOBAL_SCALAR) ||
            ((strcmp(svmstyle,"vector") == 0 || strcmp(svmstyle,"peratomtype") == 0) && data_style != FIXPROPERTY_GLOBAL_VECTOR) ||
            ((strcmp(svmstyle,"matrix") == 0 || strcmp(svmstyle,"peratomtypepair") == 0) && data_style != FIXPROPERTY_GLOBAL_MATRIX)
        )
        {
            if(errflag)
            {
                sprintf(errmsg,"%s style required for fix property/global variable %s for usage with %s",svmstyle,varname,caller);
                error->fix_error(FLERR,this,errmsg);
            }
            else return NULL;
        }

        // check length
        if((nvalues < len1) && ((data_style != FIXPROPERTY_GLOBAL_MATRIX) || ((data_style == FIXPROPERTY_GLOBAL_MATRIX) && (size_array_cols < len2))))
        {
            if(errflag)
            {
                sprintf(errmsg,"Length not sufficient for variable %s for usage with %s",varname,caller);
                error->fix_error(FLERR,this,errmsg);
            }
            else return NULL;
        }

        // success
        return static_cast<Fix*>(this);
    }
    return NULL;
}

/* ---------------------------------------------------------------------- */

void FixPropertyGlobal::init()
{
    me = comm->me;

    char errmsg[300];
    int ntypes = atom->ntypes;

    if(FIXPROPERTY_GLOBAL_VECTOR == data_style && is_atomtype_bound && nvalues != ntypes)
    {
        
        sprintf(errmsg,"Fix property/global: Length not correct for variable %s, length should be equal to %d (= number of atom types)",
                variablename,ntypes);
        error->fix_error(FLERR,this,errmsg);
    }
    if(FIXPROPERTY_GLOBAL_MATRIX == data_style && is_atomtype_bound && nvalues != ntypes*ntypes)
    {
        sprintf(errmsg,"Fix property/global: Length not correct for variable %s, length should be equal to %d ( = number of atom types * number of atom types)",
                variablename,ntypes*ntypes);
        error->fix_error(FLERR,this,errmsg);
    }
}

/* ---------------------------------------------------------------------- */

void FixPropertyGlobal::grow(int len1, int len2)
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->fix_error(FLERR, this, "Cannot grow containers of non-constant properties");

    property_->grow(len1, len2);
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_scalar()
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->fix_error(FLERR, this, "Cannot compute scalar of non-constant property");
    return property_->compute_scalar(0.);
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_vector(int i)
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->fix_error(FLERR, this, "Cannot compute vector of non-constant property");
    return property_->compute_vector(0., i);
}

void FixPropertyGlobal::vector_modify(int i,double val)
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->fix_error(FLERR, this, "Cannot modify vector of non-constant property");
    static_cast<PropertyTypeConstant*>(property_)->vector_modify(i, val);
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_array(int i, int j) //i is row, j is column
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->fix_error(FLERR, this, "Cannot compute array of non-constant property");
    return property_->compute_array(0., i, j);
}

void FixPropertyGlobal::array_modify(int i, int j,double val) //i is row, j is column
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->fix_error(FLERR, this, "Cannot modify array of non-constant property");
    static_cast<PropertyTypeConstant*>(property_)->array_modify(i, j, val);
}

/* ---------------------------------------------------------------------- */

int FixPropertyGlobal::setmask()
{
    int mask = 0;
    return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPropertyGlobal::memory_usage()
{
    double bytes = property_->memory_usage();
    return bytes;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyGlobal::new_array(int l1,int l2)
{
    
    if (data_style == FIXPROPERTY_GLOBAL_MATRIX)
        error->fix_error(FLERR,this,"Can not allocate extra array for matrix style");
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->fix_error(FLERR, this, "Can not create new array for non-constant properties");
    array_flag = 1;
    static_cast<PropertyTypeConstant*>(property_)->new_array(l1, l2);
}

/* ----------------------------------------------------------------------
   write out command
------------------------------------------------------------------------- */

void FixPropertyGlobal::write()
{
    
    if(0 != me)
        return;

    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->one(FLERR, "Fix property/global write is only allowed for type constant properties");

    FILE *file = fopen(filename,"w");

    if(!file)
        error->one(FLERR,"Fix property/global cannot open file");

    // fix id group style variablename
    fprintf(file,"fix %s %s %s %s ",id,grpname,style,variablename);

    // datatype
    const char *datatype;
    if(data_style == 0)
        datatype = "scalar";
    if(data_style == 1)
        datatype = "vector";
    if(data_style == 2 && is_symmetric)
        datatype = "atomtypepair";
    else if(data_style == 2)
        datatype = "matrix";
    fprintf(file,"%s ",datatype);

    // size_array_cols if required
    if(2 == data_style)
        fprintf(file,"%d ",size_array_cols);

    // values
    for(int i = 0; i < nvalues; i++)
        fprintf(file,"%f ", property_->compute_vector(0.,i));

    fprintf(file,"\n");
    fclose(file);
}

/* ---------------------------------------------------------------------- */

int FixPropertyGlobal::modify_param(int narg, char **arg)
{
    if (strcmp(arg[0],"file") == 0)
    {
        if (narg < 2)
            error->fix_error(FLERR,this,"not enough arguments for fix_modify 'file'");

        filename = new char[strlen(arg[1])+1];
        strcpy(filename,arg[1]);
        grpname = new char[strlen(group->names[igroup])+1];
        strcpy(grpname,group->names[igroup]);
        return 2;
    }

    return 0;
}

double FixPropertyGlobal::operator() (const double x)
{
    return property_->compute_scalar(x);
}

double FixPropertyGlobal::operator() (const double x, const int i)
{
    
    return property_->compute_vector(x, i-1);
}

double FixPropertyGlobal::operator() (const double x, const int i, const int j)
{
    
    return property_->compute_array(x, i-1, j-1);
}

double* FixPropertyGlobal::get_values() const
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->one(FLERR, "get_values called for non-constant type of property/global");
    return static_cast<PropertyTypeConstant*>(property_)->get_values();
}

double* FixPropertyGlobal::get_values_modified() const
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->one(FLERR, "get_values_modified called for non-constant type of property/global");
    return static_cast<PropertyTypeConstant*>(property_)->get_values_modified();
}

double const* const* FixPropertyGlobal::get_array() const
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->one(FLERR, "get_array called for non-constant type of property/global");
    return static_cast<PropertyTypeConstant*>(property_)->get_array();
}

double const* const* FixPropertyGlobal::get_array_modified() const
{
    if (data_type_ != FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
        error->one(FLERR, "get_array_modified called for non-constant type of property/global");
    return static_cast<PropertyTypeConstant*>(property_)->get_array_modified();
}
