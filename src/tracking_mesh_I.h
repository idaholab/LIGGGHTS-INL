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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_TRACKING_MESH_I_H
#define LMP_TRACKING_MESH_I_H

  /* ----------------------------------------------------------------------
   constructor
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  TrackingMesh<NUM_NODES>::TrackingMesh(LAMMPS *lmp)
  : MultiNodeMeshParallel<NUM_NODES>(lmp),
    customValues_(*(new CustomValueTracker(lmp,this))),
    id_ (*this->prop().template addElementProperty< ScalarContainer<int> >("id","comm_exchange_borders"/*ID does never change*/,"frame_invariant","restart_yes")),
    lineNo_(this->prop().template addElementProperty< ScalarContainer<int> >("lineNo","comm_none"/*so deleting after setup does not interefere*/,"frame_invariant","restart_no")),
    mapTagMax_(0),
    verbose_(false)
  {
  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  TrackingMesh<NUM_NODES>::~TrackingMesh()
  {
     delete &customValues_;

     // deallocate map memory if exists
     if(!mapArray_.empty()) clearMap();
  }

  /* ----------------------------------------------------------------------
   add / delete element
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool TrackingMesh<NUM_NODES>::addElement(double **nodeToAdd,int lineNumb)
  {
    // this function is always called in serial mode
    
    if(MultiNodeMeshParallel<NUM_NODES>::addElement(nodeToAdd))
    {
        // tracking mesh add one element, grow memory if necessary
        customValues_.addZeroElement();

        // set ID for element
        // ID starts from 0
        id_(this->sizeLocal()-1) = this->sizeLocal()-1;
        if(lineNo_)
            (*lineNo_)(this->sizeLocal()-1) = lineNumb;

        return true;
    }
    return false;
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::deleteElement(int n)
  {
    MultiNodeMeshParallel<NUM_NODES>::deleteElement(n);

    // tracking mesh delete code
    customValues_.deleteElement(n);
  }

  /* ----------------------------------------------------------------------
   reset global properties to original value
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::setVerbose()
  {
    verbose_ = true;
  }

  /* ----------------------------------------------------------------------
   reset global properties to original value
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool TrackingMesh<NUM_NODES>::resetToOrig()
  {
    if(MultiNodeMesh<NUM_NODES>::resetToOrig())
    {
        customValues_.resetToOrig();
        return true;
    }
    return false;
  }

  /* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::postInitialSetup()
  {
    this->prop().removeElementProperty("lineNo");
    lineNo_ = 0;
  }

  /* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::refreshOwned(int setupFlag)
  {
    MultiNodeMeshParallel<NUM_NODES>::refreshOwned(setupFlag);
    if(setupFlag) customValues_.storeOrig();
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::refreshGhosts(int setupFlag)
  {
    MultiNodeMeshParallel<NUM_NODES>::refreshGhosts(setupFlag);
  }

  /* ----------------------------------------------------------------------
   clear and generate a global map for global-local lookup
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::clearMap()
  {
      // deallocate old memory
      mapArray_.clear();
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::generateMap()
  {
      int nall = this->sizeLocal() + this->sizeGhost();

      // deallocate old memory if exists
      if(!mapArray_.empty()) clearMap();

      // get max ID of all proc
      int idmax = id_.max(nall);
      MPI_Max_Scalar(idmax,mapTagMax_,this->world);

      // build map for owned and ghost particles
      for (int i = 0; i < nall; i++)
      {
          
          const int id = id_(i);
          mapArray_[id].push_back(i);
      }
  }

  /* ----------------------------------------------------------------------
   clear ghost data that is communicated via forward comm - called in forw comm
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::clearGhostForward(bool scale,bool translate,bool rotate)
  {
      MultiNodeMeshParallel<NUM_NODES>::clearGhostForward(scale,translate,rotate);

      // delete ghost data from container classes
      // delete only data that is communicated afterwards
      for(int i = this->sizeLocal()+this->sizeGhost()-1; i >= this->sizeLocal(); i--)
          customValues_.deleteForwardElement(i,scale,translate,rotate);

  }

  /* ----------------------------------------------------------------------
   clear reverse properties, i.e. reset all of them to 0
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::clearReverse()
  {
      MultiNodeMeshParallel<NUM_NODES>::clearReverse();
      customValues_.clearReverse(this->isScaling(),this->isTranslating(),this->isRotating());
  }

  /* ----------------------------------------------------------------------
   push / pop functions for a list of elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate) const
  {
    int buf_size = 0;
    buf_size += MultiNodeMeshParallel<NUM_NODES>::elemListBufSize(n,operation,scale,translate,rotate);
    buf_size += customValues_.elemListBufSize(n,operation,scale,translate,rotate);
    return buf_size;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::pack_comm(const int operation, const int n, const int *const list, double *const buf, const int pbc_flag, const int *const pbc)
  {
    int nsend = 0;
    
    const bool scale = this->isScaling();
    const bool translate = this->isTranslating();
    const bool rotate = this->isRotating();

    // data from mother class
    nsend += MultiNodeMeshParallel<NUM_NODES>::pack_comm(operation, n, list, buf, pbc_flag, pbc);
    // custom properties
    nsend += customValues_.pushElemListToBuffer(n, list, &buf[nsend], operation, this->properties_, pbc_flag, pbc, this->domain, scale, translate, rotate);
    return nsend;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::unpack_comm(const int operation, const int n, const int first, double *const buf)
  {
    int nrecv = 0;
    
    const bool scale = this->isScaling();
    const bool translate = this->isTranslating();
    const bool rotate = this->isRotating();

    // data to mother class
    nrecv += MultiNodeMeshParallel<NUM_NODES>::unpack_comm(operation, n, first, buf);
    // custom properties
    nrecv += customValues_.popElemListFromBuffer(first, n, &buf[nrecv], operation, this->properties_, scale, translate, rotate);
    return nrecv;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::pack_reverse(const int operation, const int n, const int first, double *const buf)
  {
    int nrecv = 0;
    
    const bool scale = this->isScaling();
    const bool translate = this->isTranslating();
    const bool rotate = this->isRotating();

    nrecv += MultiNodeMeshParallel<NUM_NODES>::pack_reverse(operation, n, first, buf);
    nrecv += customValues_.pushElemListToBufferReverse(first, n, &buf[nrecv], OPERATION_COMM_REVERSE, this->properties_, scale, translate, rotate);
    return nrecv;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::unpack_reverse(const int operation, const int n, const int *const list, double *const buf)
  {
    int nsend = 0;
    
    const bool scale = this->isScaling();
    const bool translate = this->isTranslating();
    const bool rotate = this->isRotating();

    nsend += MultiNodeMeshParallel<NUM_NODES>::unpack_reverse(operation, n, list, buf);
    nsend += customValues_.popElemListFromBufferReverse(n, list, &buf[nsend], OPERATION_COMM_REVERSE, this->properties_, scale, translate, rotate);
    return nsend;
  }

  /* ----------------------------------------------------------------------
   push / pop functions for a single element
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::elemBufSize(int operation, std::list<std::string> * properties, bool scale,bool translate,bool rotate) const
  {
    int buf_size = 0;
    buf_size += MultiNodeMeshParallel<NUM_NODES>::elemBufSize(operation, properties, scale,translate,rotate);
    
    buf_size += customValues_.elemBufSize(operation, properties, scale,translate,rotate);
    
    return buf_size;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate) const
  {
    int nsend = 0;
    nsend += MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer(n,&buf[nsend],operation,scale,translate,rotate);
    nsend += customValues_.pushElemToBuffer(n,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    nrecv += MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    
    nrecv += customValues_.popElemFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  /* ----------------------------------------------------------------------
   push / pop functions for mesh properties
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::meshPropsBufSize(int operation,bool scale,bool translate,bool rotate) const
  {
    int buf_size = 0;
    buf_size += customValues_.globalPropsBufSize(operation,scale,translate,rotate);
    
    return buf_size;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::pushMeshPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate) const
  {
    int nsend = 0;
    nsend += customValues_.pushGlobalPropsToBuffer(&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<int NUM_NODES>
  int TrackingMesh<NUM_NODES>::popMeshPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    nrecv += customValues_.popGlobalPropsFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  /* ----------------------------------------------------------------------
   move / rotate  / scale
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::move(const double * const vecTotal, const double * const vecIncremental)
  {
    
    MultiNodeMesh<NUM_NODES>::move(vecTotal, vecIncremental);
    customValues_.move(vecTotal,vecIncremental);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::move(const double * const vecIncremental)
  {
    
    MultiNodeMesh<NUM_NODES>::move(vecIncremental);
    customValues_.move(vecIncremental);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::moveElement(const int i, const double * const vecIncremental)
  {
    MultiNodeMesh<NUM_NODES>::moveElement(i,vecIncremental);
    customValues_.moveElement(i,vecIncremental);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::rotate(const double * const totalQ, const double * const dQ, const double * const origin)
  {
    double negorigin[3];
    bool trans = vectorMag3DSquared(origin) > 0.;
    vectorNegate3D(origin,negorigin);

    MultiNodeMesh<NUM_NODES>::rotate(totalQ,dQ,origin);

    if(trans) customValues_.move(negorigin);
    customValues_.rotate(totalQ,dQ);
    if(trans) customValues_.move(origin);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::rotate(const double * const dQ, const double * const origin)
  {
    double negorigin[3];
    bool trans = vectorMag3DSquared(origin) > 0.;
    vectorNegate3D(origin,negorigin);

    MultiNodeMesh<NUM_NODES>::rotate(dQ,origin);

    if(trans) customValues_.move(negorigin);
    customValues_.rotate(dQ);
    if(trans) customValues_.move(origin);
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::scale(double factor)
  {
    
    MultiNodeMesh<NUM_NODES>::scale(factor);
    customValues_.scale(factor);
  }

  /* ----------------------------------------------------------------------
   return container classes
  ------------------------------------------------------------------------- */
/*
  template <typename U>
  containerTmp(

  template<int NUM_NODES> template<typename U>
  U* TrackingMesh<NUM_NODES>::containerTmp(int len)
  {
      if(len ==1) return MultiVectorConainer
  }

  template<int NUM_NODES>
  void TrackingMesh<NUM_NODES>::move(double *vecTotal, double *vecIncremental)
  {
    
    MultiNodeMesh<NUM_NODES>::move(vecTotal, vecIncremental);
    customValues_.move(vecIncremental);
  }*/

#endif
