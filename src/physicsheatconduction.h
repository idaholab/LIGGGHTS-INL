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

#ifndef PHYSICSHEATTRANSFER_H
#define PHYSICSHEATTRANSFER_H

#include "pointers.h"

namespace LAMMPS_NS
{
// forward declaration
class FixPropertyGlobal;
class FixPropertyAtom;
class PairGran;

class PhysicsHeatConduction : protected Pointers
{
    // modes for conduction contact area calaculation
    enum contactArea { CONDUCTION_CONTACT_AREA_OVERLAP,
                       CONDUCTION_CONTACT_AREA_CONSTANT,
                       CONDUCTION_CONTACT_AREA_PROJECTION,
                       CONDUCTION_CONTACT_AREA_SUPERQUADRIC,
                       CONDUCTION_CONTACT_AREA_CONVEX };

public:
    PhysicsHeatConduction(LAMMPS *lmp);
    ~PhysicsHeatConduction();

    void parseArgs(int &iarg, int narg, char **arg);

    void init();

    void beginPass();
    double computeHtPropertiesPP(int i, int j, double r, double radsum, double* all_contact_hist, int hist_offset, double &contactArea);
    double computeHtPropertiesPW(int ip, double delta_n, double *contact_history, int atom_type_wall_, double Temp_wall, double &contactArea);

private:
    template<int CONTACTAREA> double proxyComputeHtPropertiesPP(int i, int j, double r, double radsum, double* all_contact_hist, int hist_offset, double &contactArea);
    double (PhysicsHeatConduction::*currentComputeHtPropertiesPP)(int i, int j, double r, double radsum, double* all_contact_hist, int hist_offset, double &contactArea);

    template<int CONTACTAREA> double proxyComputeHtPropertiesPW(int ip, double delta_n, double *contact_history, int atom_type_wall_, double Temp_wall, double &contactArea);
    double (PhysicsHeatConduction::*currentComputeHtPropertiesPW)(int ip, double delta_n, double *contact_history, int atom_type_wall_, double Temp_wall, double &contactArea);

    double getDeltanRatio(int itype, int jtype, double ival, double jval);
    void initAreaCorrection();

    void updatePtrs();

    PairGran *pair_gran_;

    FixPropertyAtom* fixTemp_;              // fix holding particle temperature
    double *temp_;                          // particle temperature
    FixPropertyGlobal* fixConductivity_;    // fix holding per-type conductivity
    double *conductivity_;                  // per-type conductivity

    bool areaCorrectionFlag_;               // if area correction is enabled
    contactArea areaCalculationMode_;
    bool areaShapeSquare_;                  // if the area is considered square or circular

    // area correction - softness correction
    FixPropertyGlobal* fixYOrig_;           // fix holding the original Young's modulus
    double ** Yeff_;                        // current Young's modulus
    double * nu_;                           // current Poisson's ratio
    double const* const* deltanRatio_;      // ratio between current and original overlap

    // Fixed contact area
    double fixedContactArea_;

    // contact history offsets
    int iradiusOffset_;                     // convex contact area
    int cpOffset_;                          // superquadric contact area
    int alpha1Offset_;
    int alpha2Offset_;

    bool initialized_;
};
}

#endif // PHYSICSHEATTRANSFER_H
