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

    Andreas Aigner (DCS Computing GmbH, Linz)

    Copyright 2018-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#include "physicsheatconduction.h"

#include "lammps.h"
#include "error.h"
#include "force.h"
#include "pair_gran.h"
#include "pair_gran_proxy.h"
#include "math_extra_liggghts.h"

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "superquadric.h"
#endif

#include "global_properties.h"

#include <cstring>

using namespace LAMMPS_NS;

#define SMALL_FIX_HEAT_GRAN 1e-10
#define SMALL_DEVISION 1e-20

PhysicsHeatConduction::PhysicsHeatConduction(LAMMPS *lmp) :
    Pointers(lmp),
    currentComputeHtPropertiesPP(NULL),
    currentComputeHtPropertiesPW(NULL),
    pair_gran_(NULL),
    fixTemp_(NULL),
    temp_(NULL),
    fixConductivity_(NULL),
    conductivity_(NULL),
    areaCorrectionFlag_(false),
    areaCalculationMode_(CONDUCTION_CONTACT_AREA_OVERLAP),
    areaShapeSquare_(true),
    fixYOrig_(NULL),
    Yeff_(NULL),
    nu_(NULL),
    deltanRatio_(NULL),
    fixedContactArea_(0.0),
    iradiusOffset_(-1),
    cpOffset_(-1),
    alpha1Offset_(-1),
    alpha2Offset_(-1),
    initialized_(false)
{

}

PhysicsHeatConduction::~PhysicsHeatConduction()
{
    if (conductivity_)
      delete []conductivity_;
}

void PhysicsHeatConduction::parseArgs(int &iarg, int narg, char **arg)
{
    bool hasargs = true;
    while (iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"contact_area") == 0)
        {

            if(strcmp(arg[iarg+1],"overlap") == 0)
                areaCalculationMode_ =  CONDUCTION_CONTACT_AREA_OVERLAP;
            else if(strcmp(arg[iarg+1],"projection") == 0)
                areaCalculationMode_ =  CONDUCTION_CONTACT_AREA_PROJECTION;
            else if(strcmp(arg[iarg+1],"constant") == 0)
            {
                if (iarg+3 > narg)
                    error->all(FLERR,"not enough arguments for keyword 'contact_area constant'");
                areaCalculationMode_ =  CONDUCTION_CONTACT_AREA_CONSTANT;
                fixedContactArea_ = force->numeric(FLERR,arg[iarg+2]);
                if (fixedContactArea_ <= 0.)
                    error->all(FLERR,"'contact_area constant' value must be > 0");
                iarg++;
            }
            else if (strcmp(arg[iarg+1],"superquadric") == 0)
                areaCalculationMode_ = CONDUCTION_CONTACT_AREA_SUPERQUADRIC;
            else if (strcmp(arg[iarg+1],"convex") == 0)
                areaCalculationMode_ = CONDUCTION_CONTACT_AREA_CONVEX;
            else
                error->all(FLERR,"expecting 'overlap', 'projection', 'constant', 'superquadric' or 'convex' after 'contact_area'");
            iarg += 2;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"area_shape") == 0)
        {
            if (iarg+2 > narg)
                error->all(FLERR,"not enough arguments for keyword 'area_shape'");
            if(strcmp(arg[iarg+1],"square") == 0)
                areaShapeSquare_ = true;
            else if(strcmp(arg[iarg+1],"circle") == 0)
                areaShapeSquare_ = false;
            else
                error->all(FLERR,"expecting 'square' or 'circle' after 'area_shape'");
            iarg += 2;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"area_correction") == 0)
        {
            if (iarg+2 > narg)
                error->all(FLERR,"not enough arguments for keyword 'area_correction'");
            if(strcmp(arg[iarg+1],"yes") == 0)
                areaCorrectionFlag_ = true;
            else if(strcmp(arg[iarg+1],"no") == 0)
                areaCorrectionFlag_ = false;
            else
                error->all(FLERR,"expecting 'yes' or 'no' after 'area_correction'");
            iarg += 2;
            hasargs = true;
        }
    }

    if( areaCorrectionFlag_ &&
            !(areaCalculationMode_ == CONDUCTION_CONTACT_AREA_OVERLAP ||
              areaCalculationMode_ == CONDUCTION_CONTACT_AREA_CONVEX ||
              areaCalculationMode_ == CONDUCTION_CONTACT_AREA_SUPERQUADRIC) )
        error->all(FLERR,"can use 'area_correction' only for 'contact_area = overlap | convex | superquadric'");
}

void PhysicsHeatConduction::init()
{
    pair_gran_ = static_cast<PairGran*>(force->pair_match("gran", 0));
    if (!pair_gran_)
        error->all(FLERR,"No pair style gran");

    // orginally in post_create
    fixTemp_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,"physics heatconduction"));

    const int max_type = atom->get_properties()->max_type();
    fixConductivity_ =
            static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",max_type,0,"physics heatconduction"));

    if (conductivity_) delete []conductivity_;
    conductivity_ = new double[max_type];

    // pre-calculate conductivity for possible contact material combinations
    for(int i=1;i< max_type+1; i++)
    {
        conductivity_[i-1] = fixConductivity_->compute_vector(i-1);
        if(conductivity_[i-1] < 0.)
            error->all(FLERR,"Fix heat/gran/conduction: Thermal conductivity must not be < 0");
    }

    switch (areaCalculationMode_)
    {
        case CONDUCTION_CONTACT_AREA_CONSTANT:
            currentComputeHtPropertiesPP = &PhysicsHeatConduction::proxyComputeHtPropertiesPP<CONDUCTION_CONTACT_AREA_CONSTANT>;
            currentComputeHtPropertiesPW = &PhysicsHeatConduction::proxyComputeHtPropertiesPW<CONDUCTION_CONTACT_AREA_CONSTANT>;
            break;
        case CONDUCTION_CONTACT_AREA_CONVEX:
            currentComputeHtPropertiesPP = &PhysicsHeatConduction::proxyComputeHtPropertiesPP<CONDUCTION_CONTACT_AREA_CONVEX>;
            currentComputeHtPropertiesPW = &PhysicsHeatConduction::proxyComputeHtPropertiesPW<CONDUCTION_CONTACT_AREA_CONVEX>;
            break;
        case CONDUCTION_CONTACT_AREA_OVERLAP:
            currentComputeHtPropertiesPP = &PhysicsHeatConduction::proxyComputeHtPropertiesPP<CONDUCTION_CONTACT_AREA_OVERLAP>;
            currentComputeHtPropertiesPW = &PhysicsHeatConduction::proxyComputeHtPropertiesPW<CONDUCTION_CONTACT_AREA_OVERLAP>;
            break;
        case CONDUCTION_CONTACT_AREA_PROJECTION:
            currentComputeHtPropertiesPP = &PhysicsHeatConduction::proxyComputeHtPropertiesPP<CONDUCTION_CONTACT_AREA_PROJECTION>;
            currentComputeHtPropertiesPW = &PhysicsHeatConduction::proxyComputeHtPropertiesPW<CONDUCTION_CONTACT_AREA_PROJECTION>;
            break;
        case CONDUCTION_CONTACT_AREA_SUPERQUADRIC:
            currentComputeHtPropertiesPP = &PhysicsHeatConduction::proxyComputeHtPropertiesPP<CONDUCTION_CONTACT_AREA_SUPERQUADRIC>;
            currentComputeHtPropertiesPW = &PhysicsHeatConduction::proxyComputeHtPropertiesPW<CONDUCTION_CONTACT_AREA_SUPERQUADRIC>;
            break;
    }

    // calculate heat transfer correction

    if(areaCorrectionFlag_) // PHYS TODO
    {
        if(!force->pair_match("gran",0))
            error->all(FLERR,"area correction only works with using granular pair styles");

        initAreaCorrection();
    }

    updatePtrs();

    initialized_ = true;
}

void PhysicsHeatConduction::beginPass()
{
    updatePtrs();

    // find offset
    iradiusOffset_ = static_cast<PairGranProxy*>(pair_gran_)->get_history_offset("intersection_radius");
    if (areaCalculationMode_ == CONDUCTION_CONTACT_AREA_CONVEX && iradiusOffset_ < 0)
        error->all(FLERR,"Internal error: need surface model convexhull/manifold");

    cpOffset_ = static_cast<PairGranProxy*>(pair_gran_)->get_history_offset("contact_point_offset");
    alpha1Offset_ = static_cast<PairGranProxy*>(pair_gran_)->get_history_offset("alpha1_offset");
    alpha2Offset_ = static_cast<PairGranProxy*>(pair_gran_)->get_history_offset("alpha2_offset");
    if (areaCalculationMode_ == CONDUCTION_CONTACT_AREA_SUPERQUADRIC && ( cpOffset_ < 0 || alpha1Offset_ < 0 || alpha2Offset_ < 0) )
        error->all(FLERR,"Internal error: need surface model superquadric");
}

double PhysicsHeatConduction::computeHtPropertiesPP(int i, int j, double r, double radsum, double* all_contact_hist, int hist_offset, double &contactArea)
{
    if (!initialized_)
        error->all(FLERR,"not initialized");
    return (this->*currentComputeHtPropertiesPP)(i,j,r,radsum,all_contact_hist,hist_offset, contactArea);
}

double PhysicsHeatConduction::computeHtPropertiesPW(int ip, double delta_n, double *contact_history, int atom_type_wall_, double Temp_wall, double &contactArea)
{
    if (!initialized_)
        error->all(FLERR,"not initialized");
    return (this->*currentComputeHtPropertiesPW)(ip,delta_n,contact_history,atom_type_wall_,Temp_wall,contactArea);
}

double PhysicsHeatConduction::getDeltanRatio(int itype, int jtype, double ival, double jval)
{
    // no correction
    if (!areaCorrectionFlag_)
        return 1.0;

    // in case of constant, just return pre-calculated value
    if (deltanRatio_)
        return  deltanRatio_[itype-1][jtype-1];

    // else calculate the correct value
    const double expo = 1./pair_gran_->stressStrainExponent();
    const double Yeff_orig_ij = 1./((1.-pow(nu_[itype],2.))/(*fixYOrig_)(ival,itype)+(1.-pow(nu_[jtype],2.))/(*fixYOrig_)(jval,jtype));
    const double ratio = pow(Yeff_[itype][jtype]/Yeff_orig_ij,expo);
    
    return ratio;
}

/* ----------------------------------------------------------------------
   precalculate properties for area correction
------------------------------------------------------------------------- */

void PhysicsHeatConduction::initAreaCorrection()
{

    const int max_type = atom->get_properties()->max_type();

    // get the original Young's modulus - don't use registry here, since it may be a lookup table
    fixYOrig_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,"physics heatconduction"));

    // we need at least an effective Young's modulus, Poisson's ratio
    PropertyRegistry *registry = &force->registry;
    registry->registerProperty("heatPoissonsRatio", &MODEL_PARAMS::createPoissonsRatio);
    registry->connect("heatPoissonsRatio",nu_,"physics heatconduction");

    registry->registerProperty("heatYeff", &MODEL_PARAMS::createYeff);
    registry->connect("heatYeff",Yeff_,"physics heatconduction");

    // pre-calculate correction factor for constant type
    if (fixYOrig_->get_data_type() == FIXPROPERTY_GLOBAL_TYPE_CONSTANT)
    {
        const double expo = 1./pair_gran_->stressStrainExponent();
        const double *Y_orig = fixYOrig_->get_values();

        // allocate a new array within youngsModulusOriginal
        fixYOrig_->new_array(max_type,max_type);

        // feed deltan_ratio into this array
        for(int i = 1; i < max_type+1; i++)
        {
            for(int j = 1; j < max_type+1; j++)
            {
                // care different index for registered properties and direct data access
                const double Yeff_ij      = Yeff_[i][j];
                const double Yeff_orig_ij = 1./((1.-pow(nu_[i],2.))/Y_orig[i-1]+(1.-pow(nu_[j],2.))/Y_orig[j-1]);
                const double ratio = pow(Yeff_ij/Yeff_orig_ij,expo);
                
                fixYOrig_->array_modify(i-1,j-1,ratio);
            }
        }
        // get reference to deltan_ratio
        deltanRatio_ = fixYOrig_->get_array_modified();
    }
}

void PhysicsHeatConduction::updatePtrs()
{
    temp_ = fixTemp_->vector_atom;
}

template <int CONTACTAREA>
double PhysicsHeatConduction::proxyComputeHtPropertiesPP(int i, int j, double r, double radsum, double* all_contact_hist, int hist_offset, double &contactArea)
{
    int const * const type = atom->type;
    double const * const radius = atom->radius;

    const double radi = radius[i];
    const double radj = radius[j];

    const double deltan_ratio = getDeltanRatio(type[i],type[j],temp_[i],temp_[j]);

    // calculate cross sectional area of smaller bounding sphere
    // it is used as a upper limit and for scaling in case of convex simulation
    const double rmin = std::min(radi,radj);
    const double Amin = rmin*rmin*M_PI;

    contactArea = 0.0;
    bool checkContactArea = true; // by default we will check the contact area
    if(CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP)
    {
        
        if(areaCorrectionFlag_)
        {
            double delta_n = radsum - r;
            delta_n *= deltan_ratio;
            r = radsum - delta_n;
        }

        if (r < fmax(radi, radj)) // one sphere is inside the other
        {
            // set contact area to area of smaller sphere
            contactArea = fmin(radi,radj);
            contactArea *= contactArea * M_PI;
        }
        else
            //contact area of the two spheres
            contactArea = MathExtraLiggghts::contactAreaTwoSpheres(radi,radj,r);
    }
    else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONSTANT)
    {
        contactArea = fixedContactArea_;
        checkContactArea = false; // disable the contact area check for constant contact area
    }
    else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_PROJECTION)
    {
        double rmax = std::max(radi,radj);
        contactArea = M_PI*rmax*rmax;
        checkContactArea = false; // disable the contact area check for projected contact area
    }
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_SUPERQUADRIC)
    {
        double *const *const x = atom->x;
        double *const *const quat = atom->quaternion;
        double *const *const shape = atom->shape;
        double *const *const blockiness = atom->blockiness;

        Superquadric parti = Superquadric(x[i], quat[i], shape[i], blockiness[i]);
        Superquadric partj = Superquadric(x[j], quat[j], shape[j], blockiness[j]);

        const double * const contact_history = &all_contact_hist[hist_offset];
        const double * const cp = &contact_history[cpOffset_];
        const double k1 = parti.calc_curvature_coefficient(1,cp); // gaussian radius
        const double k2 = partj.calc_curvature_coefficient(1,cp);

        const double bbScale = 10;
        double r1 = 1.0/(k1+SMALL_DEVISION);
        if (r1 > radius[i]*bbScale)
            r1 = radius[i]*bbScale;
        double r2 = 1.0/(k2+SMALL_DEVISION);
        if (r2 > radius[j]*bbScale)
            r2 = radius[j]*bbScale;

        double delta = std::fabs(contact_history[alpha1Offset_]) + std::fabs(contact_history[alpha2Offset_]);
        if (areaCorrectionFlag_)
            delta *= deltan_ratio; 
        const double dist = r1 + r2 - delta;

        contactArea = MathExtraLiggghts::contactAreaTwoSpheres(r1,r2,dist);
        // restrict contactArea to positive values; this may happen in very rare cases of extrem small contact radii (r1,r2);
        if (contactArea < 0)
            contactArea = 0.0;
    }
#endif
#ifdef CONVEX_ACTIVE_FLAG
    else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONVEX)
    {
        if (atom->shapetype_flag && iradiusOffset_ >= 0)
        {
            double *contact_history = &all_contact_hist[hist_offset];
            const double cap_radius = contact_history[iradiusOffset_];
            contactArea = cap_radius*cap_radius*M_PI;

            if (areaCorrectionFlag_)
            {
                // Rather simple check to not scale 'blocky' contacts
                // For a perfect cube-cube contact the contact area does not change
                if (contactArea/Amin < 0.5)
                    // First-order approximation: Aorig = A * deltan_ratio
                    // This correction is valid for a Hertz contact area
                    contactArea *= deltan_ratio;
            }
        }
        else
            error->all(FLERR,"Non-convex particle or not convexhull/manifold used with heattransfer mode 'convex'");
    }
#endif
    // limit contactArea by cross sectional area of the smaller bounding sphere
    if ( checkContactArea && (contactArea > Amin) )
    {
        
        contactArea = Amin;
    }

    const double tcoi = conductivity_[type[i]-1];
    const double tcoj = conductivity_[type[j]-1];
    double hc;
    if (tcoi < SMALL_FIX_HEAT_GRAN || tcoj < SMALL_FIX_HEAT_GRAN)
        hc = 0.;
    else if (areaShapeSquare_)
        hc = 4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);
    else
        hc = 4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea/M_PI);

    return hc;
}

template <int CONTACTAREA>
double PhysicsHeatConduction::proxyComputeHtPropertiesPW(int ip, double delta_n, double *contact_history, int atom_type_wall_, double Temp_wall, double &contactArea)
{
    const int itype = atom->type[ip];
    const double ri = atom->radius[ip];

    // calculate cross sectional area of particle bounding sphere
    // it is used as a upper limit and for scaling in case of convex simulation
    const double Amin = ri*ri*M_PI;

    const double deltan_ratio = getDeltanRatio(itype,atom_type_wall_,temp_[ip],Temp_wall);

    contactArea = 0.0;
    if(CONTACTAREA == CONDUCTION_CONTACT_AREA_OVERLAP)
    {
        
        if (areaCorrectionFlag_)
           delta_n *= deltan_ratio;

        const double r = ri - delta_n;

        contactArea = (ri*ri-r*r)*M_PI; //contact area sphere-wall 
    }
    else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONSTANT)
        contactArea = fixedContactArea_;
    else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_PROJECTION)
    {
        contactArea = M_PI*ri*ri;
    }
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_SUPERQUADRIC)
    {
        Superquadric parti = Superquadric(atom->x[ip], atom->quaternion[ip], atom->shape[ip], atom->blockiness[ip]);

        const double * const cp = &contact_history[cpOffset_];
        const double k1 = parti.calc_curvature_coefficient(1,cp); // gaussian radius

        const double bbScale = 10;
        const double maxRad = ri * bbScale;
        double r1 = 1.0/(k1+SMALL_DEVISION);
        if (r1 > maxRad)
            r1 = maxRad;

        if (areaCorrectionFlag_)
            delta_n *= deltan_ratio; 
        const double dist = r1 - delta_n;

        contactArea = (r1*r1-dist*dist)*M_PI; 
    }
#endif
#ifdef CONVEX_ACTIVE_FLAG
    else if (CONTACTAREA == CONDUCTION_CONTACT_AREA_CONVEX)
    {
        const double cap_radius = contact_history[iradiusOffset_];
        contactArea = cap_radius*cap_radius*M_PI;

        if (areaCorrectionFlag_)
        {
            // Rather simple check to not scale 'blocky' contacts
            // For a perfect cube-cube contact the contact area does not change
            if (contactArea/Amin < 0.5)
                // First-order approximation: Aorig = A * deltan_ratio
                // This correction is valid for a Hertz contact area
                contactArea *= deltan_ratio;
        }
    }
#endif
    const double tcop = conductivity_[itype-1]; //types start at 1, array at 0
    const double tcowall = conductivity_[atom_type_wall_-1];

    double hc;
    if ((fabs(tcop) < SMALL_FIX_HEAT_GRAN) || (fabs(tcowall) < SMALL_FIX_HEAT_GRAN))
        hc = 0.;
    else if (areaShapeSquare_)
        hc = 4.*tcop*tcowall/(tcop+tcowall)*sqrt(contactArea);
    else
        hc = 4.*tcop*tcowall/(tcop+tcowall)*sqrt(contactArea/M_PI);

    return hc;
}
