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

#ifndef LMP_MULTI_NODE_MESH_PARALLEL_H
#define LMP_MULTI_NODE_MESH_PARALLEL_H

#include "mpi_liggghts.h"
#include "multi_node_mesh.h"
#include "comm.h"
#include "comm_tiled.h"
#include "comm_brick.h"
#include "error.h"
#include "vector_liggghts.h"
#include "neighbor.h"
#include "math_extra_liggghts.h"
#include "container_base.h"
#include "domain_wedge.h"
#include "parallel_base.h"
#include <string>
#include <list>
#include <cmath>
#include <algorithm>

namespace LAMMPS_NS
{

  template<int NUM_NODES>
  class MultiNodeMeshParallel : public MultiNodeMesh<NUM_NODES>, public ParallelBase
  {
      public:

        void initialSetup();
        void deleteUnowned();
        void pbcExchangeBorders(int setupFlag);
        void clearReverse();

        void forwardComm(std::string property);
        void forwardComm(std::list<std::string> * properties = NULL);
        // low level exchange functions that are called via callback from comm
        void get_nrecv_scale(int operation, double (&nrecv_scale)[2]);

        void reverseComm(std::string property);
        void reverseComm(std::list<std::string> * properties = NULL);

        double get_cutoff();
        double get_slab_offset(const int dim)
        { return this->neighbor->cutneighmax*0.5 + this->isMoving() ? this->neighbor->skin*0.5 : 0.; }

        void clear_exchange();
        int pack_exchange(const int dim);
        void unpack_exchange(const int nrecv, const int dim);
        void post_exchange();
        void add_nghost(const int new_ghosts)
        { nGhost_ += new_ghosts; }

        int pack_comm(const int operation, const int n, const int *const list, double *const buf, const int pbc_flag, const int *const pbc);
        int unpack_comm(const int operation, const int n, const int first, double *const buf);

        int pack_reverse(const int operation, const int n, const int first, double *const buf);
        int unpack_reverse(const int operation, const int n, const int *const list, double *const buf);

        bool do_migrate() const
        { return doParallellization_; }
        int pack_irregular(const int i, double *const buf);
        int unpack_irregular(double *const buf);
        void generate_map()
        { this->generateMap(); }

        virtual double get_rcb_weight(const int i)
        { return get_rcb_default_weight(); }

        virtual double get_rcb_default_weight()
        { return this->isMoving() ? 3.0 : 1.0; }

        void shrink_box(double (&shrink)[6]);
        void prepare();

        void write_restart_serial(FILE *fp) const;
        void restart_serial(double *list);
        void write_restart_parallel(FILE *fp) const;
        void restart_parallel(double *list);

        int size_restart() const;
        int pack_restart(double *const buf) const;
        int unpack_restart(double *const buf);
        void finalize_restart();

        bool allNodesInsideSimulationBox();
        void useAsInsertionMesh(bool parallel);

        inline bool isInsertionMesh() const
        { return isInsertionMesh_; }

        inline int sizeLocal() const
        { return nLocal_; }

        inline int sizeGhost() const
        { return nGhost_; }

        inline int sizeGlobal() const
        { return nGlobal_; }

        inline int sizeGlobalOrig() const
        { return nGlobalOrig_; }

        bool doParallelization() const
        { return doParallellization_; }

        inline bool isParallel() const
        { return isParallel_; }

        virtual int id(int i) = 0;

        int get_ntotal() const
        { return nLocal_+nGhost_; }

        int get_nlocal() const
        { return nLocal_; }

        double get_pos(const int i, const int dim) const
        { return this->center_(i)[dim]; }

        bool has_radius() const
        { return true; }

        double get_radius(const int i) const
        { return this->rBound_(i); }

        void set_type(const int type)
        { type_ = type; }

        int get_type(const int i) const
        { return type_; }

      protected:

        MultiNodeMeshParallel(LAMMPS *lmp);
        virtual ~MultiNodeMeshParallel();

        virtual bool addElement(double **nodeToAdd);
        
        virtual bool addElement(double **nodeToAdd,int lineNumb) = 0;
        virtual void deleteElement(int n);

        virtual void buildNeighbours() = 0;

        virtual void refreshOwned(int setupFlag);
        virtual void refreshGhosts(int setupFlag);

        virtual void qualityCheck() = 0;

        virtual void preSetup() {}
        virtual void preInitialSetup() {}
        virtual void postInitialSetup() {}
        virtual void postBorders() {}

        virtual void clearMap() = 0;
        virtual void generateMap() = 0;
        virtual void clearGhostForward(bool scale,bool translate,bool rotate);

        // lo-level parallelization also used by derived classes

        virtual int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate) const;

        virtual int elemBufSize(int operation, std::list<std::string> * properties, bool scale,bool translate,bool rotate) const;
        virtual int pushElemToBuffer(int i, double *buf,int operation,bool scale,bool translate,bool rotate) const;
        virtual int popElemFromBuffer(double *buf,int operation,bool scale,bool translate,bool rotate);

        virtual int meshPropsBufSize(int operation,bool scale,bool translate,bool rotate) const = 0;
        virtual int pushMeshPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate) const = 0;
        virtual int popMeshPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate) = 0;

        // flags if mesh should be parallelized
        bool doParallellization_;

        std::list<std::string> *properties_;    // properties requested in communication

      private:

        // parallelization functions

        void setup();
        void pbc();

        void exchange();

        void borders();
        void clearGhosts();

        int sizeRestartMesh() const;
        int sizeRestartElement() const;

        // number of local and ghost elements
        int nLocal_;
        int nGhost_;
        int nGlobal_;
        int nGlobal_restart_; // temporary buffer - used only in case of restart

        // initial number of elements
        int nGlobalOrig_;

        // flags if mesh is parallelized
        bool isParallel_;

        // flag indicating usage as insertion mesh
        bool isInsertionMesh_;

        int type_;

        // *************************************
        // comm stuff - similar to Comm class
        // *************************************

        void grow_swap(int n);
        void allocate_swap(int n);
        void free_swap();
        void grow_send(int,int);
        void grow_recv(int);

        // communication buffers
        int maxsend_,maxrecv_;       // current size of send/recv buffer
        double *buf_send_, *buf_recv_;

        // comm
        int sendneed_[3][2];         // # of procs away I send elements to
        int maxneed_[3];             // max procs away any proc needs, per dim
        double half_atom_cut_;       // half atom neigh cut

        int size_exchange_;          // # of per-element datums in exchange
        int size_forward_;           // # of per-element datums in forward comm
        int size_reverse_;           // # of datums in reverse comm
        int size_border_;            // # of datums in forward border comm

        int maxforward_,maxreverse_; // max # of datums in forward/reverse comm

        // comm swaps
        
        int nswap_;                  // # of swaps to perform = sum of maxneed
        int maxswap_;                // # of swaps where mem was allocated
        int *sendnum_,*recvnum_;     // # of atoms to send/recv in each swap
        int *sendproc_,*recvproc_;   // proc to send/recv to/from at each swap
        int *size_forward_recv_;     // # of values to recv in each forward comm
        int *size_reverse_recv_;     // # to recv in each reverse comm

        double *slablo_,*slabhi_;
        int **sendwraplist_;         // whether an element needs to be wrapped or not
        int *maxsendlist_;           // max size of send list for each swap

        int *pbc_flag_;              // general flag for sending atoms thru PBC
        int **pbc_;                  // dimension flags for PBC adjustments

        Comm *mycomm_;
  };

  // *************************************
  #include "multi_node_mesh_parallel_I.h"
  #include "multi_node_mesh_parallel_buffer_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* MULTINODEMESH_H_ */
