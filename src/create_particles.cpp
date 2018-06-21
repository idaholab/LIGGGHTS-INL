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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2017-     DCS Computing GmbH, Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

#include <cmath>
#include <stdlib.h>
#include <string.h>
#include "create_particles.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "lattice.h"
#include "region.h"
#include "random_park.h"
#include "error.h"
#include "force.h"

#include "fix_template_sphere.h"
#include "particleToInsert.h"

// include last to ensure correct macros
#include "domain_definitions.h"

using namespace LAMMPS_NS;

#define EPSILON  1.0e-6

enum{REGION,SINGLE};

/* ---------------------------------------------------------------------- */

CreateParticles::CreateParticles(LAMMPS *lmp) :
    Pointers(lmp)
{
    vectorZeroize3D(vel);
    vectorZeroize3D(omega);
    quatIdentity4D(quat);
}

/* ---------------------------------------------------------------------- */

void CreateParticles::command(int narg, char **arg)
{
    if (domain->box_exist == 0)
        error->all(FLERR,"create_particles command before simulation box is defined");
    if (modify->nfix_restart_peratom)
        error->all(FLERR,"Cannot create_particles after "
                         "reading restart file with per-atom info");

    // parse arguments

    if (narg < 2) error->all(FLERR,"Illegal create_particles command");
    int ifix = modify->find_fix(arg[0]);

    if(ifix < 0)
        error->all(FLERR,"invalid ID for fix particletemplate provided");

    if(strncmp(modify->fix[ifix]->style,"particletemplate/",16))
        error->all(FLERR,"fix is not of type particletemplate");

    pt = static_cast<FixTemplateSphere*>(modify->fix[ifix]);

    int iarg = 0;
    if (strcmp(arg[1],"region") == 0)
    {
        style = REGION;
        if (narg < 3)
            error->all(FLERR,"Illegal create_particles command");
        nregion = domain->find_region(arg[2]);
        if (nregion == -1)
            error->all(FLERR,"Create_particles region ID does not exist");
        domain->regions[nregion]->init();
        iarg = 3;
    }
    else if (strcmp(arg[1],"single") == 0)
    {
        style = SINGLE;
        if (narg < 5)
            error->all(FLERR,"Illegal create_particles command");
        xone[0] = force->numeric(FLERR,arg[2]);
        xone[1] = force->numeric(FLERR,arg[3]);
        xone[2] = force->numeric(FLERR,arg[4]);
        iarg = 5;
    } else error->all(FLERR,"Illegal create_particles command");

    // process optional keywords

    int scaleflag = 0;

    nbasis = domain->lattice->nbasis;

    while (iarg < narg) {
        if (strcmp(arg[iarg],"units") == 0) {
            if (iarg+2 > narg) error->all(FLERR,"Illegal create_particles command");
            if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
            else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
            else error->all(FLERR,"Illegal create_particles command");
            iarg += 2;
        } else if (strcmp(arg[iarg],"velocity") == 0) {
            if (iarg+4 > narg) error->all(FLERR,"Illegal create_particles command");
            vel[0] = force->numeric(FLERR,arg[iarg+1]);
            vel[1] = force->numeric(FLERR,arg[iarg+2]);
            vel[2] = force->numeric(FLERR,arg[iarg+3]);
            iarg += 4;
        } else if (strcmp(arg[iarg],"omega") == 0) {
            if (iarg+4 > narg) error->all(FLERR,"Illegal create_particles command");
            omega[0] = force->numeric(FLERR,arg[iarg+1]);
            omega[1] = force->numeric(FLERR,arg[iarg+2]);
            omega[2] = force->numeric(FLERR,arg[iarg+3]);
            iarg += 4;
        } else if (strcmp(arg[iarg],"orientation") == 0) {
            if (iarg+5 > narg) error->all(FLERR,"Illegal create_particles command");
            quat[0] = force->numeric(FLERR,arg[iarg+1]);
            quat[1] = force->numeric(FLERR,arg[iarg+2]);
            quat[2] = force->numeric(FLERR,arg[iarg+3]);
            quat[3] = force->numeric(FLERR,arg[iarg+4]);
            MathExtraLiggghts::quat_normalize(quat);
            iarg += 5;
        } else if (strcmp(arg[iarg],"rotate") == 0) {
            if (iarg+7 > narg) error->all(FLERR,"Illegal create_particles command: not enough arguments");
            if (strcmp(arg[iarg+1],"axis") != 0)
                error->all(FLERR, "Expected keyword 'axis' after keyword 'rotate'");
            double axis[3], quatRot[4], temp[4];
            axis[0] = force->numeric(FLERR,arg[iarg+2]);
            axis[1] = force->numeric(FLERR,arg[iarg+3]);
            axis[2] = force->numeric(FLERR,arg[iarg+4]);
            if (strcmp(arg[iarg+5],"angle") != 0)
                error->all(FLERR, "Expected keyword 'angle' after axis definition");
            double angle = M_PI/180*force->numeric(FLERR,arg[iarg+6]);
            if (MathExtraLiggghts::compDouble(vectorMag3DSquared(axis),0.0))
                error->all(FLERR, "invalid axis for rotation");
            MathExtra::norm3(axis);
            MathExtra::axisangle_to_quat(axis,angle,quatRot);
            MathExtra::quatquat(quatRot,quat,temp);
            vectorCopy4D(temp,quat);
            iarg += 7;
        } else error->all(FLERR,"Illegal create_particles command");
    }

    // error checks

    // setup scaling for SINGLE
    // could use domain->lattice->lattice2box() to do conversion of
    //   lattice to box, but not consistent with other uses of units=lattice
    // triclinic remapping occurs in add_single()

    if (style == REGION)
    {
        if (nbasis == 0)
            error->all(FLERR,"Cannot create particles with undefined lattice");
    }
    else if (scaleflag == 1)
    {
        xone[0] *= domain->lattice->xlattice;
        xone[1] *= domain->lattice->ylattice;
        xone[2] *= domain->lattice->zlattice;
    }

    // set bounds for my proc in sublo[3] & subhi[3]
    // if periodic and style = BOX or REGION, i.e. using lattice:
    //   should create exactly 1 atom when 2 images are both "on" the boundary
    //   either image may be slightly inside/outside true box due to round-off
    //   if I am lo proc, decrement lower bound by EPSILON
    //     this will insure lo image is created
    //   if I am hi proc, decrement upper bound by 2.0*EPSILON
    //     this will insure hi image is not created
    //   thus insertion box is EPSILON smaller than true box
    //     and is shifted away from true boundary
    //     which is where atoms are likely to be generated

    triclinic = domain->triclinic;

    double epsilon[3];
    if (triclinic) epsilon[0] = epsilon[1] = epsilon[2] = EPSILON;
    else {
        epsilon[0] = domain->prd[0] * EPSILON;
        epsilon[1] = domain->prd[1] * EPSILON;
        epsilon[2] = domain->prd[2] * EPSILON;
    }

    if (triclinic == 0) {
        sublo[0] = domain->sublo[0]; subhi[0] = domain->subhi[0];
        sublo[1] = domain->sublo[1]; subhi[1] = domain->subhi[1];
        sublo[2] = domain->sublo[2]; subhi[2] = domain->subhi[2];
    } else {
        sublo[0] = domain->sublo_lamda[0]; subhi[0] = domain->subhi_lamda[0];
        sublo[1] = domain->sublo_lamda[1]; subhi[1] = domain->subhi_lamda[1];
        sublo[2] = domain->sublo_lamda[2]; subhi[2] = domain->subhi_lamda[2];
    }

    // call setup to initialize 'just_created' fixes
    // in particular, this is important for fix multisphere, since it would
    // overwrite body informations in new multisphere particles
    modify->init_fixes();

    // allow fixes to e.g. update some pointers before set_arrays is called
    // set_arrays called in ParticleToInsert::insert()
    int nfix = modify->nfix;
    Fix **fix = modify->fix;
    for (int j = 0; j < nfix; j++)
        if (fix[j]->create_attribute) fix[j]->pre_set_arrays();

    // currently only single/region is supported
    bigint natoms_previous = atom->natoms;

    if (style == SINGLE)
        add_single();
    else
        add_lattice();

    //no need to call set_arrays here, since pti->insert() did it for us

    // clean up

    // new total # of atoms

    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (atom->natoms < 0 || atom->natoms > MAXBIGINT)
        error->all(FLERR,"Too many total atoms");

    // print status

    if (comm->me == 0) {
        if (screen)
            fprintf(screen,"Created " BIGINT_FORMAT " atoms\n",
                    atom->natoms-natoms_previous);
        if (logfile)
            fprintf(logfile,"Created " BIGINT_FORMAT " atoms\n",
                    atom->natoms-natoms_previous);
    }

    // reset simulation now that more atoms are defined
    // add tags for newly created atoms if possible
    // if global map exists, reset it
    // change these to MAXTAGINT when allow tagint = bigint

    if (atom->natoms > MAXSMALLINT) {
        if (comm->me == 0)
            error->warning(FLERR,"Total atom count exceeds ID limit, "
                                 "atoms will not have individual IDs");
        atom->tag_enable = 0;
    }
    if (atom->natoms <= MAXSMALLINT) atom->tag_extend();

    if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
    }
}

/* ----------------------------------------------------------------------
   add single atom with coords at xone if it's in my sub-box
   if triclinic, xone is in lamda coords
------------------------------------------------------------------------- */

void CreateParticles::add_single()
{

    // if triclinic, convert to lamda coords (0-1)

    double lamda[3],*coord;
    if (triclinic) {
        domain->x2lamda(xone,lamda);
        coord = lamda;
    } else coord = xone;

    // if atom is in my subbox, create it

    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
            coord[1] >= sublo[1] && coord[1] < subhi[1] &&
            coord[2] >= sublo[2] && coord[2] < subhi[2])
    {
        pt->randomize_single();
        ParticleToInsert *pti = pt->pti;
        pti->set_x_v_omega(xone,vel,omega,quat);
        pti->insert();
        //no need to call set_arrays here, since pti->insert() did it for us
    }
    pt->finalize_insertion();
}

/* ----------------------------------------------------------------------
   add many atoms by looping over lattice
------------------------------------------------------------------------- */

void CreateParticles::add_lattice()
{
    // convert 8 corners of my subdomain from box coords to lattice coords
    // for orthogonal, use corner pts of my subbox
    // for triclinic, use bounding box of my subbox
    // xyz min to max = bounding box around the domain corners in lattice space

    double bboxlo[3],bboxhi[3];

    if (triclinic == 0) {
        bboxlo[0] = domain->sublo[0]; bboxhi[0] = domain->subhi[0];
        bboxlo[1] = domain->sublo[1]; bboxhi[1] = domain->subhi[1];
        bboxlo[2] = domain->sublo[2]; bboxhi[2] = domain->subhi[2];
    } else domain->bbox(domain->sublo_lamda,domain->subhi_lamda,bboxlo,bboxhi);

    double xmin,ymin,zmin,xmax,ymax,zmax;
    xmin = ymin = zmin = BIG;
    xmax = ymax = zmax = -BIG;

    domain->lattice->bbox(1,bboxlo[0],bboxlo[1],bboxlo[2],
            xmin,ymin,zmin,xmax,ymax,zmax);
    domain->lattice->bbox(1,bboxhi[0],bboxlo[1],bboxlo[2],
            xmin,ymin,zmin,xmax,ymax,zmax);
    domain->lattice->bbox(1,bboxlo[0],bboxhi[1],bboxlo[2],
            xmin,ymin,zmin,xmax,ymax,zmax);
    domain->lattice->bbox(1,bboxhi[0],bboxhi[1],bboxlo[2],
            xmin,ymin,zmin,xmax,ymax,zmax);
    domain->lattice->bbox(1,bboxlo[0],bboxlo[1],bboxhi[2],
            xmin,ymin,zmin,xmax,ymax,zmax);
    domain->lattice->bbox(1,bboxhi[0],bboxlo[1],bboxhi[2],
            xmin,ymin,zmin,xmax,ymax,zmax);
    domain->lattice->bbox(1,bboxlo[0],bboxhi[1],bboxhi[2],
            xmin,ymin,zmin,xmax,ymax,zmax);
    domain->lattice->bbox(1,bboxhi[0],bboxhi[1],bboxhi[2],
            xmin,ymin,zmin,xmax,ymax,zmax);

    // ilo:ihi,jlo:jhi,klo:khi = loop bounds for lattice overlap of my subbox
    // overlap = any part of a unit cell (face,edge,pt) in common with my subbox
    // in lattice space, subbox is a tilted box
    // but bbox of subbox is aligned with lattice axes
    // so ilo:khi unit cells should completely tile bounding box
    // decrement lo, increment hi to avoid round-off issues in lattice->bbox(),
    //   which can lead to missing atoms in rare cases
    // extra decrement of lo if min < 0, since static_cast(-1.5) = -1

    int ilo,ihi,jlo,jhi,klo,khi;
    ilo = static_cast<int> (xmin) - 1;
    jlo = static_cast<int> (ymin) - 1;
    klo = static_cast<int> (zmin) - 1;
    ihi = static_cast<int> (xmax) + 1;
    jhi = static_cast<int> (ymax) + 1;
    khi = static_cast<int> (zmax) + 1;

    if (xmin < 0.0) ilo--;
    if (ymin < 0.0) jlo--;
    if (zmin < 0.0) klo--;

    // iterate on 3d periodic lattice of unit cells using loop bounds
    // iterate on nbasis atoms in each unit cell
    // convert lattice coords to box coords
    // add atom if it meets all criteria

    double **basis = domain->lattice->basis;
    double x[3],lamda[3];
    double *coord;

    int i,j,k,m;
    for (k = klo; k <= khi; k++)
        for (j = jlo; j <= jhi; j++)
            for (i = ilo; i <= ihi; i++)
                for (m = 0; m < nbasis; m++) {

                    x[0] = i + basis[m][0];
                    x[1] = j + basis[m][1];
                    x[2] = k + basis[m][2];

                    // convert from lattice coords to box coords

                    domain->lattice->lattice2box(x[0],x[1],x[2]);

                    // if a region was specified, test if atom is in it

                    if (style == REGION) 
                    {
                        if (!domain->regions[nregion]->match(x[0],x[1],x[2])) continue;
                    }

                    // test if atom is in my subbox

                    if (triclinic) {
                        domain->x2lamda(x,lamda);
                        coord = lamda;
                    } else coord = x;

                    if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
                        coord[1] >= sublo[1] && coord[1] < subhi[1] &&
                        coord[2] >= sublo[2] && coord[2] < subhi[2])
                    {
                        // add the atom to my list of atoms

                        pt->randomize_single();
                        ParticleToInsert *pti = pt->pti;
                        pti->set_x_v_omega(x,vel,omega,quat);
                        pti->insert();
                        //no need to call set_arrays here, since pti->insert() did it for us
                    }
                    pt->finalize_insertion();
                }
}
