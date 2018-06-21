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

    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2018-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_PARALLEL_BASE_H
#define LMP_PARALLEL_BASE_H

#include <string>
#include <list>

namespace LAMMPS_NS
{

// buffer operation types (for push and pop)
enum{ OPERATION_COMM_EXCHANGE,
    OPERATION_COMM_BORDERS,
    OPERATION_COMM_FORWARD,
    OPERATION_COMM_REVERSE,
    OPERATION_RESTART,
    OPERATION_UNDEFINED};

class ParallelBase
{
public:
    ParallelBase() :
        is_ready(false)
    {}

    ~ParallelBase() {}

    // setup functions
    virtual double get_cutoff()
    { return -1.; }
    virtual double get_slab_offset(const int dim)
    { return -1.; }

    // exchange functions
    virtual void clear_exchange() = 0;
    virtual int pack_exchange(const int dim) = 0;
    virtual void unpack_exchange(const int nrecv, const int dim) = 0;
    virtual void post_exchange() = 0;

    // forward comm functions
    virtual void get_nrecv_scale(int operation, double (&nrecv_scale)[2])
    { nrecv_scale[0] = 1; nrecv_scale[1] = 1; }

    // pack functions for forward type comm
    virtual int pack_comm(const int operation, const int n, const int *const list, double *const buf, const int pbc_flag, const int *const pbc)
    { return -1; }
    virtual int unpack_comm(const int operation, const int n, const int first, double *const buf)
    { return -1; }

    // and reverse type comm
    virtual int pack_reverse(const int operation, const int n, const int first, double *const buf)
    { return -1; }
    virtual int unpack_reverse(const int operation, const int n, const int *const list, double *const buf)
    { return -1; }

    // borders functions
    virtual int get_ntotal() const
    { return -1; }
    virtual int get_nlocal() const
    { return -1; }
    virtual double get_pos(const int i, const int dim) const
    { return -1.; }
    virtual bool has_radius() const
    { return false; }
    virtual void add_nghost(const int new_ghosts)
    { }
    virtual double get_radius(const int i) const
    { return 0.; }
    virtual int get_type(const int i) const
    { return 0; }

    // irregular functions
    virtual bool do_migrate() const
    { return true; }
    virtual int pack_irregular(const int i, double *const buf)
    { return -1; }
    virtual int unpack_irregular(double *const buf)
    { return -1; }
    virtual void generate_map()
    { }

    // restart functions
    virtual int size_restart() const = 0;
    virtual int pack_restart(double *const buf) const = 0;
    virtual int unpack_restart(double *const buf) = 0;
    virtual void finalize_restart() = 0;

    // rcb functions
    virtual double get_rcb_weight(const int i)
    { return 1.; }
    virtual double get_rcb_default_weight()
    { return 1.; }
    virtual void shrink_box(double (&shrinkt)[6])
    { }

    // general functions
    virtual void prepare()
    { }

protected:
    bool is_ready;
};

}

#endif
