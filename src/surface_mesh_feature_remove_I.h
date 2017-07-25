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

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_SURFACE_MESH_FEATURE_REMOVE_I_H
#define LMP_SURFACE_MESH_FEATURE_REMOVE_I_H

/* ----------------------------------------------------------------------
   handle exclusion of small features
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::handleExclusion(int *idListVisited)
{
    
    int nall = this->sizeLocal()+this->sizeGhost();

    if(MultiNodeMesh<NUM_NODES>::minFeatureLength() > 0. && MultiNodeMesh<NUM_NODES>::elementExclusionList())
    {
        for(int i = 0; i < nall; i++)
        {
            int nIdListVisited = 0;
            if(!checkFeatureRecursive(i,nIdListVisited,idListVisited,10))
            {
                //fprintf(this->screen,"excluded i %d\n",i);
                fprintf(MultiNodeMesh<NUM_NODES>::elementExclusionList(),"%d\n",TrackingMesh<NUM_NODES>::lineNo(i));
            }
        }
    }
}

/* ----------------------------------------------------------------------
   handle exclusion of small features - recursive algorithm
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::checkFeatureRecursive(int iSrf,
        int &nIdListVisited,int *idListVisited,int i_rec)
{
    // check if I have been here already
    for(int i = 0; i < nIdListVisited; i++)
        if(idListVisited[i] == TrackingMesh<NUM_NODES>::id(iSrf))
            return false;

    // add to visited list
    idListVisited[nIdListVisited++] = TrackingMesh<NUM_NODES>::id(iSrf);

    if(MultiNodeMesh<NUM_NODES>::rBound_(iSrf) > MultiNodeMesh<NUM_NODES>::minFeatureLength())
        return true;

    if(i_rec < 0)
        return false;

    for(int iN = 0; iN < nNeighs_(iSrf); iN++)
    {
        
        int idNeigh = neighFaces_(iSrf)[iN];
        if(idNeigh < 0) return false;
        const int nTri_j = this->map_size(idNeigh);
        for (int j = 0; j < nTri_j; j++)
        {
            int iNeigh = this->map(idNeigh, j);
            if(iNeigh >= 0 && checkFeatureRecursive(iNeigh,nIdListVisited,idListVisited,i_rec-1))
                return true;
        }
    }

    return false;
}

#endif
