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

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include <cmath>
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_force.h"
#include "fix_property_atom.h"
#include "fix_multisphere.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingForce::FixCfdCouplingForce(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg),
    CAddRhoFluid_(0.0),
    fix_dragforce_(0),
    fix_dragforce_implicit_(0),
    fix_dragforce_total_(0),
    fix_hdtorque_(0),
    fix_hdtorque_implicit_(0),
    fix_coupling_(0),
    use_torque_(false),
    dragforce_implicit_(true)
{
    iarg = 3;

    //check if style is not /implicit
    if (strcmp(style,"couple/cfd/force/implicit") == 0) error->fix_error(FLERR,this,"old fix style couple/cfd/force/implicit used, which is no longer supported.\nPlease use the regular couple/cfd/force and the fix/nve/cfd_cn/ class integrators for implicit drag handling.");

    // register props, which are always used -> x v r
    //std::string propName, bool push, bool pull, enum type:
    //{SCALAR,VECTOR,VECTOR2D,QUATERNION,SCALARMULTISPHERE,VECTORMULTISPHERE,SCALARGLOB,VECTORGLOB,MATRIXGLOB}

    registerProp("x",true,false,VECTOR);
    registerProp("v",true,false,VECTOR);
    registerProp("radius",true,false,SCALAR);
    registerProp("dragforce",false,true,VECTOR);

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;
        if(strcmp(arg[iarg],"force") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'force'");
            iarg++;
            if(strcmp(arg[iarg],"implicit") == 0) dragforce_implicit_ = true;
            else if(strcmp(arg[iarg],"explicit") == 0) dragforce_implicit_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'implicit' or 'explicit' after 'force'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"torque") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'torque'");
            iarg++;
            if(strcmp(arg[iarg],"implicit") == 0)
            {
                registerProp("KslRotation",false,true,VECTOR);
                registerProp("hdtorque",false,true,VECTOR);
                registerProp("hdtorque_implicit",false,true,VECTOR);
            }
            else if(strcmp(arg[iarg],"explicit") == 0)
            {
                 registerProp("hdtorque",false,true,VECTOR);
            }
            else
                error->fix_error(FLERR,this,"expecting 'implicit' or 'explicit' after 'torque'");
            iarg++;
            use_torque_=true;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_superquadric") == 0)
        {
            if(narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'transfer_superquadric'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
            {
                registerProp("hdtorque",false,true,VECTOR);
                registerProp("KslExtra",false,true,VECTOR);
                registerProp("ex",false,true,VECTOR);
                registerProp("volume",true,false,SCALAR);
                registerProp("area",true,false,SCALAR);
                registerProp("shape",true,false,VECTOR);
                registerProp("blockiness",true,false,VECTOR2D);
                registerProp("quaternion",true,false,QUATERNION);
            }
            else if(strcmp(arg[iarg],"no") == 0) { }
            else
              error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_superquadric'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_ellipsoid") == 0)
        {
            if(narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'transfer_ellipsoid'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
            {
                registerProp("hdtorque",false,true,VECTOR);
                registerProp("KslExtra",false,true,VECTOR);
                registerProp("ex",false,true,VECTOR);
                registerProp("shape",true,false,VECTOR);

            }
            else if(strcmp(arg[iarg],"no") == 0) { }
            else
              error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_ellipsoid'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_stochastic") == 0)
        {
            if(narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'transfer_stochastic'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
            {
                registerProp("dispersionTime",false,true,SCALAR);
                registerProp("dispersionVel",false,true,VECTOR);

            }
            else if(strcmp(arg[iarg],"no") == 0) { }
            else
              error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_stochastic'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_property") == 0)
        {
            std::string name,stype;
            cfdCoupleType type;
            if(narg < iarg+5)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_type'");
            iarg++;
            if(strcmp(arg[iarg],"name"))
                error->fix_error(FLERR,this,"expecting 'name' after 'transfer_property'");
            iarg++;
            name = arg[iarg];
            iarg++;
            if(strcmp(arg[iarg],"type"))
                error->fix_error(FLERR,this,"expecting 'type' after property name");
            iarg++;
            stype=arg[iarg];
            iarg++;

            //{SCALAR,VECTOR,VECTOR2D,QUATERNION,SCALARMULTISPHERE,VECTORMULTISPHERE,SCALARGLOB,VECTORGLOB,MATRIXGLOB}
            if (stype.compare("scalar-atom") == 0) type=SCALAR;
            else if (stype.compare("vector-atom") == 0) type=VECTOR;
            else if (stype.compare("vector2d-atom") == 0) type=VECTOR2D;
            else if (stype.compare("quaternion-atom") == 0) type=QUATERNION;
            else if (stype.compare("scalar-multisphere") == 0) type=SCALARMULTISPHERE;
            else if (stype.compare("vector-multisphere") == 0) type=VECTORMULTISPHERE;
            else if (stype.compare("scalar-global") == 0) type=SCALARGLOB;
            else if (stype.compare("vector-global") == 0) type=VECTORGLOB;
            else if (stype.compare("matrix-global") == 0) type=MATRIXGLOB;
            else error->fix_error(FLERR,this,"unknown property type");

            registerProp(name,true,false,type);

            hasargs = true;
        }
        else if (strcmp(arg[iarg],"CAddRhoFluid") == 0)
        {
            if(narg < iarg+3)
                error->fix_error(FLERR,this,"not enough arguments for 'CAddRhoFluid'. You must specify the added mass coefficient AND the fluid density");
            iarg++;
            double CAdd = atof(arg[iarg]);
            iarg++;
            double fluidDensity = atof(arg[iarg]);
            if (comm->me == 0 && screen) fprintf(screen,"fix cfd/coupling/force will consider added mass with CAdd = %g, fluidDensity: %g\n",CAdd, fluidDensity);
            CAddRhoFluid_ = CAdd*fluidDensity;
            iarg++;
            hasargs = true;
        }
        else  if (strcmp(this->style,"couple/cfd/force") == 0) error->all(FLERR,"Illegal fix cfd coupling force command");
    }
    if (dragforce_implicit_)
    {
        registerProp("Ksl",false,true,SCALAR);
        registerProp("uf",false,true,VECTOR);
        registerProp("dragforce_implicit",false,false,VECTOR);
    }

    if (force->typeSpecificCG()) registerProp("type",true,false,SCALAR);

    // flags for vector output
    vector_flag = 1;
    size_vector = 6;
    global_freq = 1;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingForce::~FixCfdCouplingForce()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::post_create()
{
    //loop all mapped entries and check if fixes exist, if not create
    //coupleParams shall have fix property pointer, bool push, bool pull, cfdCoupleType type
    //{SCALAR,VECTOR,VECTOR2D,QUATERNION,SCALARMULTISPHERE,VECTORMULTISPHERE,SCALARGLOB,VECTORGLOB,MATRIXGLOB}
    std::map<std::string, coupleParams >::iterator it;
    for ( it = coupleList_.begin(); it != coupleList_.end(); it++ )
    {
        //look for fix name
        it->second.fix_property_ = modify->find_fix_id((it->first).c_str());
        bool isMS = ( it->second.type == SCALARMULTISPHERE ||  it->second.type == VECTORMULTISPHERE );
        if (!it->second.fix_property_ && !isMS ) //not found, create fix, if not multisphere
        {
            const char* fixarg[12];
            fixarg[0]=(it->first).c_str();
            fixarg[1]="all";
            if (it->second.type <= QUATERNION) fixarg[2]="property/atom";
            else error->fix_error(FLERR,this,"not yet implemented ;)");
            fixarg[3]=(it->first).c_str();
            if (it->second.type == SCALAR) fixarg[4]="scalar";
            else if (it->second.type == VECTOR) fixarg[4]="vector";
            else if (it->second.type == VECTOR2D) fixarg[4]="vector";
            else if (it->second.type == QUATERNION) fixarg[4]="vector";
            else error->fix_error(FLERR,this,"not yet implemented ;)");
            fixarg[5]="no";     // restart
            fixarg[6]="no";     // communicate ghost
            fixarg[7]="no";     // communicate rev
            fixarg[8]="0.";
            fixarg[9]="0.";
            fixarg[10]="0.";
            fixarg[11]="0.";
            
            int narg;
            if (it->second.type == SCALAR || it->second.type == SCALARMULTISPHERE || it->second.type == SCALARGLOB) narg = 9;
            else if (it->second.type == VECTOR2D) narg = 10;
            else if (it->second.type == VECTOR || it->second.type == VECTORMULTISPHERE || it->second.type == VECTORMULTISPHERE) narg = 11;
            else narg = 12;

            it->second.fix_property_ = modify->add_fix_property_atom(narg,const_cast<char**>(fixarg),style);
        }
        else //multisphere treatment
        {
            class FixMultisphere *fix_multisphere = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
            if(!fix_multisphere) error->fix_error(FLERR,this,"this fix with these properties needs a fix of style multisphere defined before this fix in the input script");
            int n_body,len1,len2;
            // check if property exists, custom value tracker throws error if it exists
            if(! fix_multisphere->extract_ms((it->first).c_str(),len1,len2))
            {
                n_body = fix_multisphere->data().n_body();
                if ( it->second.type == SCALARMULTISPHERE )
                {
                    fix_multisphere->data().prop().addElementProperty< ScalarContainer<double> >((it->first).c_str(),"comm_exchange_borders","frame_invariant","restart_no",1,n_body)->setAllToZero();
                    fix_multisphere->data().prop().getElementProperty< ScalarContainer<double> >((it->first).c_str())->setDefaultValue(0.0);
                }
                else
                {
                    fix_multisphere->data().prop().addElementProperty< VectorContainer<double,3> >((it->first).c_str(),"comm_exchange_borders","frame_invariant","restart_no",1,n_body)->setAllToZero();
                    fix_multisphere->data().prop().getElementProperty< VectorContainer<double,3> >((it->first).c_str())->setDefaultValue(0.0);
                }
            }
        }
    }

    fix_dragforce_ = (FixPropertyAtom*)(coupleList_.find("dragforce")->second.fix_property_);
    if (coupleList_.count("dragforce_implicit") > 0)
        fix_dragforce_implicit_= (FixPropertyAtom*)(coupleList_.find("dragforce_implicit")->second.fix_property_);
    if (use_torque_)
        fix_hdtorque_ = (FixPropertyAtom*)(coupleList_.find("hdtorque")->second.fix_property_);
    if (coupleList_.count("hdtorque_implicit") > 0)
        fix_hdtorque_implicit_= (FixPropertyAtom*)(coupleList_.find("hdtorque_implicit")->second.fix_property_);
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if (fix_dragforce_)
            modify->delete_fix("dragforce");
        if (fix_dragforce_implicit_)
            modify->delete_fix("dragforce_implicit");
        if (fix_hdtorque_)
            modify->delete_fix("hdtorque");
        if (fix_hdtorque_implicit_)
            modify->delete_fix("hdtorque_implicit");
    }
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingForce::setmask()
{
    int mask = 0;
    mask |= POST_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");
    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/force needs a fix of type couple/cfd");

    //loop all mapped entries and add push / pull properties
    //coupleParams shall have fix property pointer, bool push, bool pull, cfdCoupleType type
    //{SCALAR,VECTOR,VECTOR2D,QUATERNION,SCALARMULTISPHERE,VECTORMULTISPHERE,SCALARGLOB,VECTORGLOB,MATRIXGLOB}
    const char* cfdCoupleTypeNames[9] = {"scalar-atom","vector-atom","vector2D-atom","quaternion-atom","scalar-multisphere","vector-multisphere","scalar-global","vector-global","matrix-global"};
    std::map<std::string, coupleParams >::iterator it;
    for ( it = coupleList_.begin(); it != coupleList_.end(); it++ )
    {
        if (it->second.push)
        {
              
              fix_coupling_->add_push_property(it->first.c_str(),cfdCoupleTypeNames[it->second.type]);
        }
        if (it->second.pull)
        {
              
              fix_coupling_->add_pull_property(it->first.c_str(),cfdCoupleTypeNames[it->second.type]);
        }
    }

    vectorZeroize3D(dragforce_total);
    vectorZeroize3D(hdtorque_total);

    if (strcmp(update->integrate_style,"respa") == 0)
       error->fix_error(FLERR,this,"'run_style respa' not supported.");
    fix_dragforce_total_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dragforce_total","property/atom","vector",0,0,style,false));
    if (!fix_dragforce_total_)
    {
          const char* fixarg[11];
          fixarg[0]="dragforce_total";
          fixarg[1]="all";
          fixarg[2]="property/atom";
          fixarg[3]="dragforce_total";
          fixarg[4]="vector";
          fixarg[5]="no";    // restart
          fixarg[6]="no";     // communicate ghost
          fixarg[7]="no";     // communicate rev
          fixarg[8]="0.";
          fixarg[9]="0.";
          fixarg[10]="0.";
          fix_dragforce_total_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }

    //check if being used with cfd_cn class integrator
    if (dragforce_implicit_)
    {
        if ( ! modify->find_fix_style("nve/cfd_cn",0)) error->warning(FLERR,"Trying to use implicit drag formulation without integrator style nve/cfd_cn! The drag will be ignored.");
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::setup(int vflag)
{
    if (strstr(update->integrate_style,"verlet"))
        post_force(vflag);
    else
        error->fix_error(FLERR,this,"only 'run_style verlet' supported.");
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::post_force(int)
{
    //explicit treatment of forces
    double **f = atom->f;
    double **torque = atom->torque;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double **dragforce = fix_dragforce_->array_atom;
    double **hdtorque;
    if(use_torque_) hdtorque = fix_hdtorque_->array_atom;
    double **dragforce_implicit;
    if (fix_dragforce_implicit_) dragforce_implicit = fix_dragforce_implicit_->array_atom;
    double **hdtorque_implicit;
    if (fix_hdtorque_implicit_) hdtorque_implicit = fix_hdtorque_implicit_->array_atom;
    double **dragforce_sum;
    if (fix_dragforce_total_) dragforce_sum = fix_dragforce_total_->array_atom;

    fix_dragforce_total_->set_all(0.0, true);
    vectorZeroize3D(dragforce_total);
    vectorZeroize3D(hdtorque_total);

    // add dragforce to force vector
    
    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            vectorAdd3D(f[i], dragforce[i], f[i]);
            if(use_torque_) vectorAdd3D(torque[i], hdtorque[i], torque[i]);
            vectorAdd3D(dragforce_sum[i], dragforce[i], dragforce_sum[i]);
            if (fix_dragforce_implicit_) vectorAdd3D(dragforce_sum[i], dragforce_implicit[i], dragforce_sum[i]);
            vectorAdd3D(dragforce_total, dragforce_sum[i], dragforce_total);
            if(use_torque_)
            {
                vectorAdd3D(hdtorque_total,hdtorque[i],hdtorque_total);
                if (fix_hdtorque_implicit_) vectorAdd3D(hdtorque_total,hdtorque_implicit[i],hdtorque_total);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForce::registerProp(std::string propName, bool push, bool pull, cfdCoupleType type)
{
    std::map<std::string, coupleParams >::const_iterator it = coupleList_.find(propName);
    if ( it == coupleList_.end())
    {
        if(comm->me == 0 && screen && push)
            fprintf(screen,"fix/couple/cfd/force: adding a property to push %s \n",propName.c_str());
        if(comm->me == 0 && screen && pull)
            fprintf(screen,"fix/couple/cfd/force: adding a property to pull %s \n",propName.c_str());
        coupleParams newParams = {NULL, push, pull, type};
        coupleList_.insert(std::pair<std::string, coupleParams>(propName, newParams));
    }
}

/* ----------------------------------------------------------------------
   return components of total force on fix group
------------------------------------------------------------------------- */

double FixCfdCouplingForce::compute_vector(int n)
{
  if(n < 3)
  {
    double dragtotal = dragforce_total[n];
    MPI_Sum_Scalar(dragtotal,world);
    return dragtotal;
  }

  double hdtorque = hdtorque_total[n-3];
  MPI_Sum_Scalar(hdtorque,world);
  return hdtorque;
}
