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
    Arno Mayrhofer (DCS Computing GmbH, Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_MULTI_NODE_MESH_PARALLEL_BUFFER_I_H
#define LMP_MULTI_NODE_MESH_PARALLEL_BUFFER_I_H

  /* ----------------------------------------------------------------------
   push / pop for irregular (single element)
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pack_irregular(const int i, double *const buf)
  {
      const int nsend_this = pushElemToBuffer(i, &buf[1], OPERATION_COMM_EXCHANGE, false, false, false)+1;
      buf[0] = static_cast<double>(nsend_this);
      this->deleteElement(i);
      return nsend_this;
  }

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::unpack_irregular(double *const buf)
  {
      // number of values is first in buffer
      const int nrecv_this = static_cast<int>(buf[0]);
      popElemFromBuffer(&buf[1], OPERATION_COMM_EXCHANGE, false, false, false);
      nLocal_++;
      return nrecv_this;
  }

  /* ----------------------------------------------------------------------
   push / pop for exchange
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pack_exchange(const int dim)
  {
      // scale translate rotate not needed here
      bool dummy = false;

      double checklo = this->domain->sublo[dim];
      double checkhi;
      if(this->domain->subhi[dim] == this->domain->boxhi[dim])
        checkhi = this->domain->boxhi[dim] + SMALL_DMBRDR;
      else
        checkhi = this->domain->subhi[dim];

      int nsend = 0, nsend_this = 0;
      int i = 0;
      double *buf_send = mycomm_->get_buf_send();
      CommTiled *const commTiled = (mycomm_->style == 1) ? static_cast<CommTiled*>(mycomm_) : NULL;

      while(i < nLocal_)
      {
          if(!(this->center_(i)[dim] >= checklo && this->center_(i)[dim] < checkhi))
          {
              if (commTiled)
              {
                  const int proc = commTiled->evaluate_point_drop(dim, this->center_(i));
                  if (proc != mycomm_->me)
                      buf_send[nsend++] = proc;
              }
              nsend_this = pushElemToBuffer(i,&(buf_send[nsend+1]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);
              buf_send[nsend] = static_cast<double>(nsend_this+1);
              nsend += (nsend_this+1);
              
              buf_send = mycomm_->check_grow_send(nsend, 1, BUFEXTRA_MNMP);
              this->deleteElement(i); 
          }
          else i++;
      }
      return nsend;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::unpack_exchange(const int nrecv, const int dim)
  {
      double center_elem[3];
      double checklo,checkhi;
      int m = 0, nrecv_this;

      // scale translate rotate not needed here
      bool dummy = false;

      checklo = this->domain->sublo[dim];
      if(this->domain->subhi[dim] == this->domain->boxhi[dim])
        checkhi = this->domain->boxhi[dim] + SMALL_DMBRDR;
      else
        checkhi = this->domain->subhi[dim];

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
          // number of values is first in buffer
          nrecv_this = static_cast<int>(buf[m]);
          if (do_check)
          {
              // center is next in buffer, test it
              vectorCopy3D(&(buf[m+1]),center_elem);

              if(center_elem[dim] >= checklo && center_elem[dim] < checkhi)
              {
                popElemFromBuffer(&(buf[m+1]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);
                nLocal_++;
                
              }
              
          }
          m += nrecv_this;
      }
  }

  /* ----------------------------------------------------------------------
   restart functionality - write all required data into restart buffer
   executed on all processes, but only proc 0 writes into writebuf
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::write_restart_serial(FILE *fp) const
  {
      int size_this;

      // # elements
      int nlocal = this->sizeLocal();
      int nglobal = sizeGlobal();

      // buffer sizes
      int sizeMesh, sizeElements, sizeElements_all;

      sizeMesh = sizeRestartMesh();
      sizeElements = nlocal * (sizeRestartElement() + 1); 

      double *bufMesh = NULL, *sendbufElems = NULL, *recvbufElems = NULL;
      bool dummy = false;

      // pack global data into buffer
      // do this only on proc 0
      if(this->comm->me == 0)
      {
          this->memory->create(bufMesh,sizeMesh,"MultiNodeMeshParallel::writeRestart:bufMesh");
          pushMeshPropsToBuffer(bufMesh, OPERATION_RESTART,dummy,dummy,dummy);
      }

      // allocate send buffer and pack element data
      // all local elements are in list
      this->memory->create(sendbufElems,sizeElements,"MultiNodeMeshParallel::writeRestart:sendbuf");
      sizeElements = 0;
      for(int i = 0; i < nlocal; i++)
      {
          size_this = pushElemToBuffer(i,&(sendbufElems[sizeElements+1]),OPERATION_RESTART,dummy,dummy,dummy);
          sendbufElems[sizeElements] = static_cast<double>(size_this+1);
          sizeElements += (size_this+1);
      }

      // gather the per-element data
      
      sizeElements_all = MPI_Gather0_Vector(sendbufElems,sizeElements,recvbufElems,this->world);

      // actually write data to restart file
      // do this only on proc 0
      if(this->comm->me == 0)
      {
        double nG = static_cast<double>(nglobal);

        // for error check
        double sE = static_cast<double>(sizeRestartElement());
        double sM = static_cast<double>(sizeRestartMesh());

        // size with 3 extra values
        int size = (sizeMesh+sizeElements_all+3) * sizeof(double);

        // write size
        fwrite(&size,sizeof(int),1,fp);

        // write 3 extra values
        fwrite(&nG,sizeof(double),1,fp);
        fwrite(&sE,sizeof(double),1,fp);
        fwrite(&sM,sizeof(double),1,fp);

        // write per-element and mesh data
        fwrite(recvbufElems,sizeof(double),sizeElements_all,fp);
        fwrite(bufMesh,sizeof(double),sizeMesh,fp);
      }

      // free mem

      if(bufMesh)
        this->memory->destroy(bufMesh);

      this->memory->destroy(sendbufElems);

      if(recvbufElems)
        delete []recvbufElems;
  }

  /* ----------------------------------------------------------------------
   restart functionality - read all required data from restart buffer
   executed on all processes
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::restart_serial(double *list)
  {
      int m, nglobal, nrecv_this, sE, sM;
      bool dummy = false;

      m = 0;

      nglobal = static_cast<int> (list[m++]);
      sE = static_cast<int> (list[m++]);
      sM = static_cast<int> (list[m++]);

      if(sE != sizeRestartElement() || sM != sizeRestartMesh())
          this->error->all(FLERR,"Incompatible mesh restart file - mesh has different properties in restarted simulation");

      for(int i = 0; i < nglobal; i++)
      {
          nrecv_this = static_cast<int>(list[m]);
          
          popElemFromBuffer(&(list[m+1]),OPERATION_RESTART,dummy,dummy,dummy);
          m += nrecv_this;
      }

      this->prop().deleteRestartGlobal(dummy,dummy,dummy);
      popMeshPropsFromBuffer(&list[m],OPERATION_RESTART,dummy,dummy,dummy);
  }

  /* ----------------------------------------------------------------------
   new restart implementation - write only global data into restart buffer
   executed on all processes, but only proc 0 writes into writebuf
   if the mesh is parallized only global data are saved
   for insertion meshes all data are writtin into the global file
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::write_restart_parallel(FILE *fp) const
  {
      bool parallel = doParallelization();

      double *bufMesh = NULL, *sendbufElems = NULL, *recvbufElems = NULL;
      bool dummy = false;
      int sizeElements_all = 0;

      // if not parallel gather per-element information from all procs
      if (!parallel)
      {
          const int nlocal = this->sizeLocal();
          int sizeElements = nlocal * (sizeRestartElement() + 1); 

          // allocate send buffer and pack element data
          // all local elements are in list
          this->memory->create(sendbufElems,sizeElements,"MultiNodeMeshParallel::writeRestart:sendbuf");
          sizeElements = 0;
          for(int i = 0; i < nlocal; i++)
          {
              const int size_this = pushElemToBuffer(i,&(sendbufElems[sizeElements+1]),OPERATION_RESTART,dummy,dummy,dummy);
              sendbufElems[sizeElements] = static_cast<double>(size_this+1);
              sizeElements += (size_this+1);
          }

          // gather the per-element data
          
          sizeElements_all = MPI_Gather0_Vector(sendbufElems,sizeElements,recvbufElems,this->world);
      }

      // pack global data
      // actually write data to restart file
      // do this only on proc 0
      if(this->comm->me == 0)
      {
        // pack global data into buffer
        // buffer sizes
        int sizeMesh = sizeRestartMesh();
        this->memory->create(bufMesh,sizeMesh,"MultiNodeMeshParallel::writeRestart:bufMesh");
        pushMeshPropsToBuffer(bufMesh, OPERATION_RESTART,dummy,dummy,dummy);

        // for error check
        double pF = static_cast<double>(parallel);
        double nG = static_cast<double>(sizeGlobal());
        double sM = static_cast<double>(sizeRestartMesh());
        
        // size with 3 extra values
        int size = 3+sizeMesh; // size, parallel flag, nGlobal, sizeMesh
        if (!parallel)
            size += sizeElements_all + 1; // #elems and size of all elements
        size *= sizeof(double);

        // write size - always first!
        fwrite(&size,sizeof(int),1,fp);

        // write 3 extra values
        fwrite(&pF,sizeof(double),1,fp);
        fwrite(&nG,sizeof(double),1,fp);
        fwrite(&sM,sizeof(double),1,fp);

        // write mesh data
        fwrite(bufMesh,sizeof(double),sizeMesh,fp);

        // only for non-parallel meshes write also per-element data
        if (!parallel)
        {
            double sE = static_cast<double>(sizeRestartElement());
            
            fwrite(&sE,sizeof(double),1,fp);
            fwrite(recvbufElems,sizeof(double),sizeElements_all,fp);
        }
      }

      // free mem

      if(bufMesh)
        this->memory->destroy(bufMesh);

      if (sendbufElems)
          this->memory->destroy(sendbufElems);

      if(recvbufElems)
        delete []recvbufElems;
  }

  /* ----------------------------------------------------------------------
   new restart implementation - read only global data from restart buffer
   executed on all processes
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::restart_parallel(double *list)
  {
      bool dummy = false;

      int m = 0;
      bool parallel = static_cast<bool>(list[m++]);
      int nglobal = static_cast<int> (list[m++]);
      int sM = static_cast<int> (list[m++]);
      
      nGlobal_restart_ = nglobal;
      nGlobalOrig_ = nglobal; // save the original nGlobal; May be overwritten

      if(sM != sizeRestartMesh())
          this->error->all(FLERR,"Incompatible mesh restart file - mesh has different properties in restarted simulation");

      this->prop().deleteRestartGlobal(dummy,dummy,dummy);
      popMeshPropsFromBuffer(&list[m],OPERATION_RESTART,dummy,dummy,dummy);

      if (!parallel)
      {
          const int sE = static_cast<int> (list[m++]);

          if(sE != sizeRestartElement())
              this->error->all(FLERR,"Incompatible mesh restart file - mesh has different properties in restarted simulation");

          for(int i = 0; i < nglobal; i++)
          {
              const int nrecv_this = static_cast<int>(list[m]);
              
              popElemFromBuffer(&(list[m+1]),OPERATION_RESTART,dummy,dummy,dummy);
              m += nrecv_this;
          }
      }
  }

  /* ----------------------------------------------------------------------
   size of local buffer for per-element data
   return 0 for non-parallized meshes
  ------------------------------------------------------------------------- */
  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::size_restart() const
  {
      if (doParallelization())
          // need 1 extra per-element-space for per-element buffer length
          // and 1 extra space for number of local elements
          return sizeLocal() * (sizeRestartElement() + 1) + 1;
      else
          return 0;
  }

  /* ----------------------------------------------------------------------
   pack only per-element data into buffer
   global data are saved during write_restart
  ------------------------------------------------------------------------- */
  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pack_restart(double *const buf) const
  {
      // # elements
      const int nlocal = this->sizeLocal();

      // pack element data
      // all local elements are in list
      int sizeElements = 0;
      buf[sizeElements++] = static_cast<double>(nlocal);
      for(int i = 0; i < nlocal; i++)
      {
          // No scale/rotate/translate info needed
          const bool dummy = false;
          int size_this = pushElemToBuffer(i,&(buf[sizeElements+1]),OPERATION_RESTART,dummy,dummy,dummy);
          buf[sizeElements] = static_cast<double>(size_this+1);
          sizeElements += (size_this+1);
      }
      return sizeElements;
  }

  /* ----------------------------------------------------------------------
   unpack only per-element data from buffer
  ------------------------------------------------------------------------- */
  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::unpack_restart(double *const buf)
  {
      int m = 0;
      const int nlocal = static_cast<int>(buf[m++]);
      
      // pack element data
      // all local elements are in list
      for(int i = 0; i < nlocal; i++)
      {
          // No scale/rotate/translate info needed
          const bool dummy = false;
          const int nrecv_this = static_cast<int>(buf[m]);
          popElemFromBuffer(&(buf[m+1]),OPERATION_RESTART,dummy,dummy,dummy);
          m += nrecv_this;
      }
      return m;
  }

  /* ----------------------------------------------------------------------
   initialize after restart data were read
  ------------------------------------------------------------------------- */
  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::finalize_restart()
  {
      // calculate nGlobal after restart
      MPI_Sum_Scalar(nLocal_,nGlobal_,this->world);
      nGlobalOrig_ = nGlobal_restart_; // restore nGlobalOrig from restart data - for error checks
  }

  /* ----------------------------------------------------------------------
   size of restart data
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::sizeRestartMesh() const
  {
      
      bool dummy = false;
      return meshPropsBufSize(OPERATION_RESTART,dummy,dummy,dummy);
  }

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::sizeRestartElement() const
  {
      
      bool dummy = false;
      return elemBufSize(OPERATION_RESTART, NULL, dummy,dummy,dummy);
  }

  /* ----------------------------------------------------------------------
   return required buffer size for a list of elements for borders(),forwardComm()
   must match push / pop implementation
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate) const
  {
      return n*elemBufSize(operation, NULL, scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   push a list of elements for borders(), forwardComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
   
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pack_comm(const int operation, const int n, const int *const list, double *const buf, const int pbc_flag, const int *const pbc)
  {
      
      int nsend = 0;
      
      const bool scale = this->isScaling();
      const bool translate = this->isTranslating();
      const bool rotate = this->isRotating();

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          
          // triangle center
          if (!this->properties_ || this->center_.matches_any_id(this->properties_))
              nsend += this->center_.pushElemListToBuffer(n, list, &(buf[nsend]), operation, pbc_flag, pbc, this->domain, scale, translate, rotate);
          // triangle vertices
          if (!this->properties_ || this->node_.matches_any_id(this->properties_))
              nsend += this->node_.pushElemListToBuffer(n, list, &(buf[nsend]), operation, pbc_flag, pbc, this->domain, scale, translate, rotate);
          // bounding radius
          if (!this->properties_ || this->rBound_.matches_any_id(this->properties_))
              nsend += this->rBound_.pushElemListToBuffer(n, list, &(buf[nsend]), operation, pbc_flag, pbc, this->domain, scale, translate, rotate);
          // original vertex positions
          if(this->node_orig_)
          {
              if (!this->properties_ || this->node_orig_->matches_any_id(this->properties_))
                  nsend += this->node_orig_->pushElemListToBuffer(n, list, &(buf[nsend]), operation, pbc_flag, pbc, this->domain, scale, translate, rotate);
          }
          return nsend;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          
          return nsend;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::pack_comm");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop a list of elements for borders(), forwardComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
   
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::unpack_comm(const int operation, const int n, const int first, double *const buf)
  {
      int nrecv = 0;

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          // triangle center
          if (!this->properties_ || this->center_.matches_any_id(this->properties_))
              nrecv += this->center_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          // triangle vertices
          if (!this->properties_ || this->node_.matches_any_id(this->properties_))
              nrecv += this->node_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          // bounding radius
          if (!this->properties_ || this->rBound_.matches_any_id(this->properties_))
              nrecv += this->rBound_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          // original vertex positions
          if(this->node_orig_)
          {
              if (!this->properties_ || this->node_orig_->matches_any_id(this->properties_))
                  nrecv += this->node_orig_->popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          }
          return nrecv;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          
          //    nrecv += MultiNodeMesh<NUM_NODES>::node_.popListFromBuffer(first,n,&(buf[nrecv]),operation);
          return nrecv;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::unpack_comm");
      return 0;
  }

  /* ----------------------------------------------------------------------
   push a list of elements for reverseComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
   
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pack_reverse(const int operation, const int n, const int first, double *const buf)
  {
      int nsend = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        
        return nsend;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::pack_reverse");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop a list of elements for reverseComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
   
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::unpack_reverse(const int operation, const int n, const int *const list, double *const buf)
  {
      int nrecv = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        
        return nrecv;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::unpack_reverse");
      return 0;
  }

  /* ----------------------------------------------------------------------
   return required buffer size for one element for exchange()
   must match push / pop implementation
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
   
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::elemBufSize(int operation, std::list<std::string> * properties, bool, bool, bool) const
  {
      int size_buf = 0;

      if(operation == OPERATION_RESTART)
      {
          if (!properties || MultiNodeMesh<NUM_NODES>::node_.matches_any_id(properties))
              size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          return size_buf;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          if (!properties || MultiNodeMesh<NUM_NODES>::center_.matches_any_id(properties))
              size_buf += MultiNodeMesh<NUM_NODES>::center_.elemBufSize();
          if (!properties || MultiNodeMesh<NUM_NODES>::node_.matches_any_id(properties))
              size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          if (!properties || MultiNodeMesh<NUM_NODES>::rBound_.matches_any_id(properties))
              size_buf += MultiNodeMesh<NUM_NODES>::rBound_.elemBufSize();
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
          {
              if (!properties || MultiNodeMesh<NUM_NODES>::node_orig_->matches_any_id(properties))
                  size_buf += MultiNodeMesh<NUM_NODES>::node_orig_->elemBufSize();
          }
          return size_buf;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          
          return size_buf;
      }

      if(operation == OPERATION_COMM_REVERSE)
      {
        
        return size_buf;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::elemBufSize");
      return 0;
  }

  /* ----------------------------------------------------------------------
   push one element for exchange()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer(int i, double *buf, int operation, bool, bool, bool) const
  {
      int nsend = 0;

      if(operation == OPERATION_RESTART)
      {
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);

          return nsend;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          
          nsend += MultiNodeMesh<NUM_NODES>::center_.pushElemToBuffer(i,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::rBound_.pushElemToBuffer(i,&(buf[nsend]),operation);
          if(this->node_orig_)
              nsend += this->node_orig_->pushElemToBuffer(i,&(buf[nsend]),operation);
          return nsend;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop one element for exchange, borders or restart
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer(double *buf, int operation, bool, bool, bool)
  {
      int nrecv = 0;

      if(operation == OPERATION_RESTART)
      {
          MultiVectorContainer<double,NUM_NODES,3> nodeTmp("nodeTmp");
          
          nrecv += nodeTmp.popElemFromBuffer(&(buf[nrecv]),operation);
          this->addElement(nodeTmp.begin()[0],-1);

          this->prop().deleteRestartElement(nLocal_-1,false,false,false);

          return nrecv;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          nrecv += MultiNodeMesh<NUM_NODES>::center_.popElemFromBuffer(&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popElemFromBuffer(&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::rBound_.popElemFromBuffer(&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popElemFromBuffer(&(buf[nrecv]),operation);
          return nrecv;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

#endif
