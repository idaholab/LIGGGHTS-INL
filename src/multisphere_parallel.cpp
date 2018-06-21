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

#define DELTA 10000

#include "multisphere_parallel.h"
#include "atom.h"
#include "atom_vec.h"
#include "vector_liggghts.h"
#include "domain.h"
#include "memory.h"
#include "comm_brick.h"
#include "comm_tiled.h"

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

/* ----------------------------------------------------------------------
   constructor / destructor
------------------------------------------------------------------------- */

MultisphereParallel::MultisphereParallel(LAMMPS *lmp) :
  Multisphere(lmp),
  mycomm_(NULL)
{
    if (this->comm->style == 0)
        mycomm_ = new CommBrick(this->comm, this);
    else
        mycomm_ = new CommTiled(this->comm, this);
}

/* ----------------------------------------------------------------------
   exchange bodies with neighbor procs
------------------------------------------------------------------------- */

int MultisphereParallel::pack_exchange(const int dim)
{
    const double lo = domain->sublo[dim];
    const double hi = domain->subhi[dim];

    int i = 0;
    int nsend = 0;
    double x[3];

    double *buf_send = mycomm_->get_buf_send();
    CommTiled *const commTiled = (mycomm_->style == 1) ? static_cast<CommTiled*>(mycomm_) : NULL;
    while (i < nbody_)
    {
          MathExtraLiggghts::local_coosys_to_cartesian(x,xcm_to_xbound_(i),ex_space_(i),ey_space_(i),ez_space_(i));
          vectorAdd3D(xcm_(i),x,x);

          if (x[dim] < lo || x[dim] >= hi)
          {
              if (commTiled)
              {
                  const int proc = commTiled->evaluate_point_drop(dim, x);
                  if (proc != mycomm_->me)
                      buf_send[nsend++] = proc;
              }
              buf_send = mycomm_->check_grow_send(nsend, 1, BUFEXTRA);
              nsend += pack_exchange_rigid(i,&buf_send[nsend]);
              
              remove_body(i);
              
          }
          else
              i++;
    }

    return nsend;
}

void MultisphereParallel::unpack_exchange(const int nrecv, const int dim)
{
    const double lo = domain->sublo[dim];
    const double hi = domain->subhi[dim];

    int m = 0;

    double *const buf = mycomm_->get_buf_recv();

    while (m < nrecv)
    {
        bool do_check = true;
        // only for tiled
        if (mycomm_->style == 1)
        {
            const int proc = static_cast<int> (buf[m++]);
            if (proc != mycomm_->me)
                do_check = false;
        }
        double value = buf[m+dim+1];
        
        if (do_check && value >= lo && value < hi)
            m += unpack_exchange_rigid(&buf[m]);
        else
            m += static_cast<int> (buf[m]);
    }
}

void MultisphereParallel::post_exchange()
{
    calc_nbody_all();
    Multisphere::generate_map();
}

/* ----------------------------------------------------------------------
   restart functionality - read all required data from restart buffer
   executed on all processes
------------------------------------------------------------------------- */

void MultisphereParallel::restart_serial(double *list)
{
    bool dummy = false;
    int m = 0, nrecv_this;

    int nbody_all_old = static_cast<int> (list[m++]);
    
    nbody_ = nbody_all_ = 0;

    for(int i = 0; i < nbody_all_old; i++)
    {
        nrecv_this = static_cast<int>(list[m]);
        
        double *x_bnd = &(list[m+1]);
        
        if(domain->is_in_subdomain(x_bnd))
        {
            
            customValues_.addZeroElement();
            customValues_.deleteRestartElement(nbody_,dummy,dummy,dummy);
            customValues_.popElemFromBuffer(&(list[m+4]),OPERATION_RESTART,dummy,dummy,dummy);

            nbody_++;
        }
        m += nrecv_this;
    }

    // do initialization tasks

    MPI_Sum_Scalar(nbody_,nbody_all_,world);
    generate_map();
    reset_forces(true);

}

/* ----------------------------------------------------------------------
   new restart functionality - write only global data into restart buffer
   executed on all processes, but only proc 0 writes into writebuf
------------------------------------------------------------------------- */

void MultisphereParallel::write_restart_parallel(FILE *fp)
{
    double nba = static_cast<double>(n_body_all());
    
    // write data to file
    if(comm->me == 0)
    {
        // size with 1 extra value (nba)
        int size = (1) * sizeof(double);
        // write size
        fwrite(&size,sizeof(int),1,fp);

        // write extra value
        fwrite(&nba,sizeof(double),1,fp);
        
    }
}

/* ----------------------------------------------------------------------
   new restart functionality - read onyl global data from restart buffer
   executed on all processes
------------------------------------------------------------------------- */

void MultisphereParallel::restart_parallel(double *list)
{
    nbody_all_restart = static_cast<int> (list[0]);
    
}

/* ----------------------------------------------------------------------
   pack multisphere for irregular communication
------------------------------------------------------------------------- */

int MultisphereParallel::pack_irregular(const int i, double * const buf)
{
    const int nsend_this = pack_exchange_rigid(i,&buf[1]);
    buf[0] = static_cast<double>(nsend_this);
    remove_body(i);
    return nsend_this;
}

int MultisphereParallel::unpack_irregular(double * const buf)
{
    // number of values is first in buffer
    const int nrecv_this = static_cast<int>(buf[0]);
    unpack_exchange_rigid(&buf[1]);
    return nrecv_this;
}

int MultisphereParallel::size_restart() const
{
    bool dummy = false;
    return n_body() * (customValues_.elemBufSize(OPERATION_RESTART, NULL, dummy,dummy,dummy) + 4) + 1; // (elemSize + elemPos[3] + elemData[]) *n_body
}

int MultisphereParallel::pack_restart(double * const buf) const
{
    int sizeLocal = 0;
    const int nb = n_body();
    buf[sizeLocal++] = static_cast<double>(nb); // currently n_body_all is used!

    for (int i = 0; i < nb; i++)
    {
        double xbnd[3];
        x_bound(xbnd,i);
        const bool dummy = false;
        int sizeOne = customValues_.pushElemToBuffer(i,&(buf[sizeLocal+4]),OPERATION_RESTART,dummy,dummy,dummy);
        buf[sizeLocal] = static_cast<double>(sizeOne+4);
        buf[sizeLocal+1] = xbnd[0];
        buf[sizeLocal+2] = xbnd[1];
        buf[sizeLocal+3] = xbnd[2];

        sizeLocal += sizeOne+4;
        
    }

    return sizeLocal;
}

int MultisphereParallel::unpack_restart(double * const buf)
{
    
    int m = 0;
    const int nbody_local = static_cast<int>(buf[m++]);
    
    for(int i = 0; i < nbody_local; i++)
    {
        const int nrecv_this = static_cast<int>(buf[m]);
        // add multisphere regardless of its position
        // it may be communicated later accordingly
        
        customValues_.addZeroElement();
        const bool dummy = false;
        customValues_.deleteRestartElement(nbody_,dummy,dummy,dummy);
        customValues_.popElemFromBuffer(&(buf[m+4]),OPERATION_RESTART,dummy,dummy,dummy);

        nbody_++;

        m += nrecv_this;
    }
    
    return m;
}

/* ----------------------------------------------------------------------
 initialize after restart data were read
------------------------------------------------------------------------- */
void MultisphereParallel::finalize_restart()
{
    MPI_Sum_Scalar(nbody_,nbody_all_,world);
    Multisphere::generate_map();
    reset_forces(true);
    
}
