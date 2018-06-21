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

#ifndef LMP_MULTI_NODE_MESH_PARALLEL_I_H
#define LMP_MULTI_NODE_MESH_PARALLEL_I_H

#define BIG_MNMP 1.0e20
#define BUFFACTOR_MNMP 1.5
#define BUFMIN_MNMP 2000
#define BUFEXTRA_MNMP 2000

  /* ----------------------------------------------------------------------
   consturctors
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  MultiNodeMeshParallel<NUM_NODES>::MultiNodeMeshParallel(LAMMPS *lmp)
  : MultiNodeMesh<NUM_NODES>(lmp),
    doParallellization_(true),
    properties_(NULL),
    nLocal_(0), nGhost_(0), nGlobal_(0), nGlobalOrig_(0),
    isParallel_(false),
    isInsertionMesh_(false),
    maxsend_(0), maxrecv_(0),
    buf_send_(0), buf_recv_(0),
    half_atom_cut_(0.),
    size_exchange_(0),
    size_forward_(0),
    size_reverse_(0),
    size_border_(0),
    maxforward_(0),maxreverse_(0),
    nswap_(0),
    maxswap_(6),
    sendnum_(0),recvnum_(0),
    sendproc_(0),recvproc_(0),
    size_forward_recv_(0),
    size_reverse_recv_(0),
    slablo_(0),slabhi_(0),
    sendwraplist_(0),
    maxsendlist_(0),
    mycomm_(NULL)
  {
      // initialize comm buffers & exchange memory
      
      maxsend_ = BUFMIN_MNMP;
      this->memory->create(buf_send_,maxsend_+BUFEXTRA_MNMP,"MultiNodeMeshParallel:buf_send");
      maxrecv_ = BUFMIN_MNMP;
      this->memory->create(buf_recv_,maxrecv_,"MultiNodeMeshParallel:buf_recv");

      allocate_swap(maxswap_);

      if (this->comm->style == 0)
          mycomm_ = new CommBrick(this->comm, this);
      else
          mycomm_ = new CommTiled(this->comm, this);
  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  MultiNodeMeshParallel<NUM_NODES>::~MultiNodeMeshParallel()
  {
      free_swap();

      this->memory->destroy(buf_send_);
      this->memory->destroy(buf_recv_);

      delete mycomm_;
  }

  /* ----------------------------------------------------------------------
   add and delete elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMeshParallel<NUM_NODES>::addElement(double **nodeToAdd)
  {
    
    if(MultiNodeMesh<NUM_NODES>::addElement(nodeToAdd))
    {
        nLocal_++;
        return true;
    }
    return false;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::deleteElement(int n)
  {
    
    if(n < nLocal_ && nGhost_ != 0)
        this->error->one(FLERR,"Illegal call to MultiNodeMeshParallel<NUM_NODES>::deleteElement");

    MultiNodeMesh<NUM_NODES>::deleteElement(n);

    if(n >= nLocal_)
        nGhost_--;
    else
        nLocal_--;
  }

  /* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::refreshOwned(int setupFlag)
  {
    MultiNodeMesh<NUM_NODES>::refreshOwned(setupFlag);
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::refreshGhosts(int setupFlag)
  {
    MultiNodeMesh<NUM_NODES>::refreshGhosts(setupFlag);
  }

  /* ----------------------------------------------------------------------
   completely clear ghosts - called in borders()
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearGhosts()
  {
      // delete ghost data from container classes

      while(nGhost_ > 0)
      {
          
          deleteElement(nLocal_);
      }
  }

  /* ----------------------------------------------------------------------
   clear ghost data that is communicated via forward comm - called in forw comm
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearGhostForward(bool scale,bool translate,bool rotate)
  {
      // delete ghost data from container classes
      // delete only data that is communicated afterwards

      for(int i = this->sizeLocal()+this->sizeGhost()-1; i >= this->sizeLocal(); i--)
      {
          // clear ghost data that belongs to this class
          // must match push/pop implementation for forward comm in this class
          if(translate || rotate || scale)
          {
            this->node_.del(i);
            this->center_.del(i);
          }
          if(scale)
            this->rBound_.del(i);
      }
  }

  /* ----------------------------------------------------------------------
   check if all elements are in domain
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMeshParallel<NUM_NODES>::allNodesInsideSimulationBox()
  {
    int flag = 0;
    for(int i=0;i<sizeLocal();i++)
      for(int j=0;j<NUM_NODES;j++)
      {
        
        if(!this->domain->is_in_domain(this->node_(i)[j]))
        {
            flag = 1;
            break;
        }
      }

    MPI_Max_Scalar(flag,this->world);
    if(flag) return false;
    else return true;
  }

  /* ----------------------------------------------------------------------
   set flag if used as insertion mesh
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::useAsInsertionMesh(bool parallelflag)
  {
    
    isInsertionMesh_ = true;
    
    if(!parallelflag)
    {
        if(isParallel())
            this->error->all(FLERR,"If a run command is between the fix mesh/surface and the "
                             "fix insert command, you have to use fix mesh/surface/planar for "
                             "the insertion mesh");
        doParallellization_ = false;
    }
  }

  /* ----------------------------------------------------------------------
   setup of communication
  ------------------------------------------------------------------------- */

   template<int NUM_NODES>
   void MultiNodeMeshParallel<NUM_NODES>::setup()
   {
       if(!doParallellization_) return;

       // get required size of communication per element
       const bool scale = this->isScaling();
       const bool translate = this->isTranslating();
       const bool rotate = this->isRotating();

       size_exchange_ = elemBufSize(OPERATION_COMM_EXCHANGE, NULL, scale, translate, rotate) + 1;
       size_border_   = elemBufSize(OPERATION_COMM_BORDERS,  NULL, scale, translate, rotate);
       size_forward_  = elemBufSize(OPERATION_COMM_FORWARD,  NULL, scale, translate, rotate);
       size_reverse_  = elemBufSize(OPERATION_COMM_REVERSE,  NULL, scale, translate, rotate);

       // maxforward = # of datums in largest forward communication
       // maxreverse = # of datums in largest reverse communication

       maxforward_ = MathExtraLiggghts::max(size_exchange_,size_border_,size_forward_);
       maxreverse_ = size_reverse_;

       mycomm_->set_sizes(size_border_, size_forward_, size_reverse_, maxforward_, maxreverse_);

       mycomm_->setup();
   }

   template<int NUM_NODES>
   double MultiNodeMeshParallel<NUM_NODES>::get_cutoff()
   {
       // ghost elements are for computing interaction with owned particles
       // so need to aquire ghost elements that overlap my subbox extened by
       // half neigh cutoff
       
       half_atom_cut_ = this->neighbor->cutneighmax / 2.;

       if(this->isMoving())
         half_atom_cut_+= this->neighbor->skin / 2.;

       // calculate maximum bounding radius of elements across all procs
       double rBound_max = 0.;
       for(int i = 0; i < sizeLocal(); i++)
           rBound_max = std::max(this->rBound_(i),rBound_max);
       MPI_Max_Scalar(rBound_max,this->world);

       // mesh element ghost cutoff is element bounding radius plus half atom neigh cut
       return rBound_max + half_atom_cut_;
   }

/* ----------------------------------------------------------------------
   realloc the buffers needed for communication and swaps
------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_swap(int n)
  {
      free_swap();
      allocate_swap(n);

      maxswap_ = n;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::allocate_swap(int n)
  {
      this->memory->create(sendnum_,n,"MultiNodeMeshParallel:sendnum_");
      this->memory->create(recvnum_,n,"MultiNodeMeshParallel:recvnum_");
      this->memory->create(sendproc_,n,"MultiNodeMeshParallel:sendproc_");
      this->memory->create(recvproc_,n,"MultiNodeMeshParallel:recvproc_");
      this->memory->create(size_forward_recv_,n,"MultiNodeMeshParallel:size");
      this->memory->create(size_reverse_recv_,n,"MultiNodeMeshParallel:size");
      this->memory->create(slablo_,n,"MultiNodeMeshParallel:slablo_");
      this->memory->create(slabhi_,n,"MultiNodeMeshParallel:slabhi_");
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::free_swap()
  {
      this->memory->destroy(sendnum_);
      this->memory->destroy(recvnum_);
      this->memory->destroy(sendproc_);
      this->memory->destroy(recvproc_);
      this->memory->destroy(size_forward_recv_);
      this->memory->destroy(size_reverse_recv_);
      this->memory->destroy(slablo_);
      this->memory->destroy(slabhi_);
  }

  /* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_send(int n, int flag)
  {
      maxsend_ = static_cast<int> (BUFFACTOR_MNMP * n);
      
      if (flag)
        this->memory->grow(buf_send_,(maxsend_+BUFEXTRA_MNMP),"MultiNodeMeshParallel:buf_send");
      else {
        this->memory->destroy(buf_send_);
        this->memory->create(buf_send_,maxsend_+BUFEXTRA_MNMP,"MultiNodeMeshParallel:buf_send");
      }
  }

  /* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_recv(int n)
  {
      maxrecv_ = static_cast<int> (BUFFACTOR_MNMP * n);
      this->memory->destroy(buf_recv_);
      this->memory->create(buf_recv_,maxrecv_,"MultiNodeMeshParallel:buf_recv");
  }

  /* ----------------------------------------------------------------------
   parallelization -
   initially, all processes have read the whole data
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::initialSetup()
  {
      this->prepare();

      // check for possible round-off isues
      double span = 0.;
      if (this->node_.size() > 0)
      {
          double max_scalar = this->node_.max_scalar();
          double min_scalar = this->node_.min_scalar();
          span = max_scalar - min_scalar;
      }
      MPI_Max_Scalar(span, this->world);
      if(span < 1e-4)
        this->error->all(FLERR,"Mesh error - root causes: (a) mesh empty or (b) dimensions too small - use different unit system");

      double comBefore[3];
      this->center_of_mass(comBefore);
      
      // delete all elements that do not belong to this processor
      
      if(sizeGlobal() != sizeGlobalOrig())
      {
        
        char errstr[1024];

        if(sizeGlobal() == 0)
        {
            sprintf(errstr,"Mesh (id %s): All %d mesh elements have been lost / left the domain. \n"
                           "Please use 'boundary m m m' or scale/translate/rotate the mesh or change its dynamics\n"
                           "FYI: center of mass of mesh including scale/tranlate/rotate is %f / %f / %f\n"
                           "     simulation box x from %f to %f y  from %f to %f z from %f to %f\n"
                           "     (gives indication about changes in scale/tranlate/rotate necessary to make simulation run)\n",
                       this->mesh_id_,sizeGlobalOrig()-sizeGlobal(),comBefore[0],comBefore[1],comBefore[2],
                       this->domain->boxlo[0],this->domain->boxhi[0],this->domain->boxlo[1],this->domain->boxhi[1],this->domain->boxlo[2],this->domain->boxhi[2]);
        }
        else
        {
            double comAfter[3];
            this->center_of_mass(comAfter);

            sprintf(errstr,"Mesh (id %s): %d mesh elements have been lost / left the domain. \n"
                           "Please use 'boundary m m m' or scale/translate/rotate the mesh or change its dynamics\n"
                           "FYI: center of mass of mesh including scale/tranlate/rotate before cutting out elements is %f / %f / %f\n"
                           "     simulation box x from %f to %f y  from %f to %f z from %f to %f\n"
                           "     center of mass of mesh after cutting out elements outside simulation box is is        %f / %f / %f\n"
                           "     (gives indication about changes in scale/tranlate/rotate necessary to make simulation run)\n",
                       this->mesh_id_,sizeGlobalOrig()-sizeGlobal(),comBefore[0],comBefore[1],comBefore[2],
                       this->domain->boxlo[0],this->domain->boxhi[0],this->domain->boxlo[1],this->domain->boxhi[1],this->domain->boxlo[2],this->domain->boxhi[2],
                       comAfter[0],comAfter[1],comAfter[2]);
        }
        this->error->all(FLERR,errstr);
      }

      // perform operations that should be done before initial setup
      
      preInitialSetup();

      // set-up mesh parallelism
      
      setup();

      // re-calculate properties for owned particles
      
      refreshOwned(1);

      // identify elements that are near borders
      // forward communicate them
      
      borders();

      // re-calculate properties for ghost particles
      
      refreshGhosts(1);

      // build mesh topology and neigh list
      
      buildNeighbours();

      // perform quality check on the mesh
      
      qualityCheck();

      if(doParallellization_) isParallel_ = true;

      postInitialSetup();

      // stuff that should be done before resuming simulation
      
      postBorders();
      
  }

  /* ----------------------------------------------------------------------
   parallelization - aggregates pbc, exchange and borders
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::pbcExchangeBorders(int setupFlag)
  {
      // need not do this during simulation for non-moving mesh and non-changing simulation box
      
      if(setupFlag) this->reset_stepLastReset();

      // perform operations that should be done before setting up parallellism and exchanging elements
      preSetup();

      if(!setupFlag && !this->isMoving() && !this->isDeforming() && !this->domain->box_change) return;

      // set-up mesh parallelism
      setup();

      // enforce pbc
      pbc();

      // communicate particles
      exchange();

      if(sizeGlobal() != sizeGlobalOrig())
      {
        
        //this->error->all(FLERR,"Mesh elements have been lost");
        char errstr[500];
        sprintf(errstr,"Mesh (id %s): Mesh elements have been lost / left the domain. Please use "
                       "'boundary m m m' or scale/translate/rotate the mesh or change its dynamics",
                       this->mesh_id_);
        this->error->all(FLERR,errstr);
      }

      // re-calculate properties for owned particles
      
      refreshOwned(setupFlag);

      // identify elements that are near borders
      // forward communicate them
      
      borders();

      // re-calculate properties for ghosts
      refreshGhosts(setupFlag);

      // stuff that should be done before resuming simulation
      postBorders();

  }

  /* ----------------------------------------------------------------------
   parallelization - clear data of reverse comm properties
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearReverse()
  {
      // nothing to do here
  }

  /* ----------------------------------------------------------------------
   delete all particles which are not owned on this proc
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::deleteUnowned()
  {
      nGlobalOrig_ = sizeLocal();

      int i = 0;

      if(doParallellization_)
      {

          while(i < nLocal_)
          {
              if(!this->domain->is_in_subdomain(this->center_(i)))
                  this->deleteElement(i);
              else i++;
          }

          // calculate nGlobal for the first time
          MPI_Sum_Scalar(nLocal_,nGlobal_,this->world);
          isParallel_ = true;
      }
      else
        nGlobal_ = nLocal_;

  }

  /* ----------------------------------------------------------------------
   enforce periodic boundary conditions
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::pbc()
  {
      if(!doParallellization_) return;

      double centerNew[3], delta[3];

      for(int i = 0; i < this->sizeLocal(); i++)
      {
          vectorCopy3D(this->center_(i),centerNew);
          this->domain->remap(centerNew);
          vectorSubtract3D(centerNew,this->center_(i),delta);

          // move element i incremental
          if(vectorMag3DSquared(delta) > 1e-9)
            this->moveElement(i,delta);
      }
  }

/* ----------------------------------------------------------------------
exchange elements with nearby processors
------------------------------------------------------------------------- */

template<int NUM_NODES>
void MultiNodeMeshParallel<NUM_NODES>::exchange()
{
    if(!doParallellization_)
        return;

    mycomm_->exchange();
}

/* ----------------------------------------------------------------------
clear before exchange
------------------------------------------------------------------------- */

template<int NUM_NODES>
void MultiNodeMeshParallel<NUM_NODES>::clear_exchange()
{
    // clear global->local map for owned and ghost atoms
    
    clearMap();

    // clear old ghosts
    
    clearGhosts();
}

/* ----------------------------------------------------------------------
work to be done after exchange
------------------------------------------------------------------------- */

template<int NUM_NODES>
void MultiNodeMeshParallel<NUM_NODES>::post_exchange()
{
    // re-calculate nGlobal as some element might have been lost
    MPI_Sum_Scalar(nLocal_,nGlobal_,this->world);
}

  /* ----------------------------------------------------------------------
   generate ghost elements, refresh global map
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::borders()
  {
      if(doParallellization_)
          mycomm_->borders();

      // build global-local map
      this->generateMap();
  }

  /* ----------------------------------------------------------------------
   communicate properties to ghost elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::forwardComm(std::string property)
  {
      std::list<std::string> properties (1, property);
      forwardComm(&properties);
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::forwardComm(std::list<std::string> * properties)
  {
      properties_ = properties;
      mycomm_->forward_comm();
      properties_ = NULL;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::get_nrecv_scale(int operation, double (&nrecv_scale)[2])
  {
      
      const bool scale = this->isScaling();
      const bool translate = this->isTranslating();
      const bool rotate = this->isRotating();

      int size_comm = 0;
      if (operation == OPERATION_COMM_FORWARD)
          size_comm = size_forward_;
      else if (operation == OPERATION_COMM_REVERSE)
          size_comm = size_reverse_;

      // nrecv_scale[0] is the size of the buffers listed in properties (if set)
      // if no properties_ are set we use either 1 or 0, the latter only in case size_comm is equal to 0
      nrecv_scale[0] = properties_ ?
                       elemBufSize(operation, properties_, scale, translate, rotate) :
                       std::min(size_comm, 1);
      // size_comm is the size of all comm (forward or reverse) buffers
      nrecv_scale[1] = properties_ ? size_comm : 1;
  }

  /* ----------------------------------------------------------------------
   reverse communication of properties on atoms every timestep
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::reverseComm(std::string property)
  {
      std::list<std::string> properties (1, property);
      reverseComm(&properties);
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::reverseComm(std::list<std::string> * properties)
  {
      properties_ = properties;
      mycomm_->reverse_comm();
      properties_ = NULL;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::shrink_box(double (&shrink)[6])
  {
      for (int i = 0; i < this->nLocal_; i++)
      {
          shrink[0] = MIN(shrink[0],this->center_(i)[0]);
          shrink[1] = MIN(shrink[1],this->center_(i)[1]);
          shrink[2] = MIN(shrink[2],this->center_(i)[2]);
          shrink[3] = MAX(shrink[3],this->center_(i)[0]);
          shrink[4] = MAX(shrink[4],this->center_(i)[1]);
          shrink[5] = MAX(shrink[5],this->center_(i)[2]);
      }
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::prepare()
  {
      if (is_ready)
          return;
      deleteUnowned();
      is_ready = true;
  }

#endif
