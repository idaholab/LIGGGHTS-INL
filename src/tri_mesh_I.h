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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_TRI_MESH_I_H
#define LMP_TRI_MESH_I_H

#ifndef SMALL_TRIMESH
#define SMALL_TRIMESH (1.e-10)  
#endif
#define LARGE_TRIMESH 1000000

/* ---------------------------------------------------------------------- */

inline double TriMesh::resolveTriSphereContact(const int iPart, const int nTri, const double rSphere, const double *const cSphere, double *delta, int &barysign)
{
    // this is the overlap algorithm, neighbor list build is
    // coded in resolveTriSphereNeighbuild

    double bary[3];
    return resolveTriSphereContactBary(iPart,nTri,rSphere,cSphere,delta,bary,barysign);
}

/* ---------------------------------------------------------------------- */

inline double TriMesh::resolveTriSphereContactBary(const int iPart, const int nTri, const double rSphere,
                                                 const double *const cSphere, double *delta, double *bary,
                                                 int &barySign,bool skip_inactive)
{
    double **n = node_(nTri);
    int obtuseAngleIndex = SurfaceMeshBase::obtuseAngleIndex(nTri);

    bary[0] = bary[1] = bary[2] = 0.;

    double node0ToSphereCenter[3];
    //double *surfNorm = SurfaceMeshBase::surfaceNorm(nTri);
    vectorSubtract3D(cSphere,n[0],node0ToSphereCenter);

    MathExtraLiggghts::calcBaryTriCoords(node0ToSphereCenter,edgeVec(nTri),edgeLen(nTri),bary);

    double invlen = 1./(2.*rBound_(nTri));
    barySign = (bary[0] > -precision_trimesh()*invlen) + 2*(bary[1] > -precision_trimesh()*invlen) + 4*(bary[2] > -precision_trimesh()*invlen);

    if (weightedWallFormulation_)
    {
        return calculateWeightedProperties(cSphere, n, SurfaceMeshBase::surfaceNorm(nTri), bary, node0ToSphereCenter, rSphere, delta);
    }
    else
    {
        
        double d(0.);

        switch(barySign)
        {
        case 1: 
            d = resolveCornerContactBary(nTri,0,obtuseAngleIndex == 0,cSphere,delta,bary,skip_inactive);
            break;
        case 2: 
            d = resolveCornerContactBary(nTri,1,obtuseAngleIndex == 1,cSphere,delta,bary,skip_inactive);
            break;
        case 3: 
            d = resolveEdgeContactBary(nTri,0,cSphere,delta,bary,skip_inactive);
            break;
        case 4: 
            d = resolveCornerContactBary(nTri,2,obtuseAngleIndex == 2,cSphere,delta,bary,skip_inactive);
            break;
        case 5: 
            d = resolveEdgeContactBary(nTri,2,cSphere,delta,bary,skip_inactive);
            break;
        case 6: 
            d = resolveEdgeContactBary(nTri,1,cSphere,delta,bary,skip_inactive);
            break;
        case 7: // face contact - all three barycentric coordinates are > 0
            d = resolveFaceContactBary(nTri,cSphere,node0ToSphereCenter,delta);
            break;
        default:
            
            this->error->one(FLERR,"Internal error");
            d = 1.; // doesn't exist, just to satisfy the compiler
            break;
        }

        // return distance - radius of the particle
        return d - rSphere;
    }
}

/* ---------------------------------------------------------------------- */

inline double TriMesh::resolveEdgeContactBary(const int iTri, const int iEdge, const double *const p,
                                            double *const delta, double *const bary, const bool skip_inactive)
{
    int ip = (iEdge+1)%3, ipp = (iEdge+2)%3;
    double nodeToP[3], d(1.);
    double **n = node_(iTri);

    vectorSubtract3D(p,n[iEdge],nodeToP);

    double distFromNode =  vectorDot3D(nodeToP,edgeVec(iTri)[iEdge]);

    if (distFromNode < -SMALL_TRIMESH)
    {
        
        if(skip_inactive && !cornerActive(iTri)[iEdge])
            return LARGE_TRIMESH;
        d = calcDist(p,n[iEdge],delta);
        bary[iEdge] = 1.; bary[ip] = 0.; bary[ipp] = 0.;
    }
    else if(distFromNode > edgeLen(iTri)[iEdge] + SMALL_TRIMESH)
    {
        
        if(skip_inactive && !cornerActive(iTri)[ip])
            return LARGE_TRIMESH;
        d = calcDist(p,n[ip],delta);
        bary[iEdge] = 0.; bary[ip] = 1.; bary[ipp] = 0.;
    }
    else
    {
        
        double closestPoint[3];

        if(skip_inactive && !edgeActive(iTri)[iEdge])
            return LARGE_TRIMESH;

        vectorAddMultiple3D(n[iEdge],distFromNode,edgeVec(iTri)[iEdge],closestPoint);

        d = calcDist(p,closestPoint,delta);

        bary[ipp] = 0.;
        bary[iEdge] = 1. - distFromNode/edgeLen(iTri)[iEdge];
        bary[ip] = 1. - bary[iEdge];
    }

    return d;
}

/* ---------------------------------------------------------------------- */

inline double TriMesh::resolveCornerContactBary(const int iTri, const int iNode, const bool obtuse,
                                                const double *const p, double *const delta,
                                                double *const bary, const bool skip_inactive)
{
    int ip = (iNode+1)%3, ipp = (iNode+2)%3;
    //double d(1.);
    double *n = node_(iTri)[iNode];

    if(obtuse)
    {
        
        double **edge = edgeVec(iTri);
        double nodeToP[3], closestPoint[3];

        vectorSubtract3D(p,n,nodeToP);

        double distFromNode = vectorDot3D(nodeToP,edge[ipp]);
        if(distFromNode < SMALL_TRIMESH)
        {
            if(distFromNode > -edgeLen(iTri)[ipp])
            {
                
                if(skip_inactive && !edgeActive(iTri)[ipp])
                    return LARGE_TRIMESH;

                vectorAddMultiple3D(n,distFromNode,edge[ipp],closestPoint);

                bary[ip] = 0.;
                bary[iNode] = 1. + distFromNode/edgeLen(iTri)[ipp];
                bary[ipp] = 1. - bary[iNode];

                return calcDist(p,closestPoint,delta);
            }
            else
            {
                
                if(skip_inactive && !cornerActive(iTri)[ipp])
                    return LARGE_TRIMESH;

                bary[ipp] = 1.; bary[iNode] = bary[ip] = 0.;
                return calcDist(p,node_(iTri)[ipp],delta);
            }
        }

        distFromNode = vectorDot3D(nodeToP,edge[iNode]);
        if(distFromNode > -SMALL_TRIMESH)
        {
            if(distFromNode < edgeLen(iTri)[iNode])
            {
                
                if(skip_inactive && !edgeActive(iTri)[iNode])
                    return LARGE_TRIMESH;

                vectorAddMultiple3D(n,distFromNode,edge[iNode],closestPoint);

                bary[ipp] = 0.;
                bary[iNode] = 1. - distFromNode/edgeLen(iTri)[iNode];
                bary[ip] = 1. - bary[iNode];

                return calcDist(p,closestPoint,delta);
            }
            else
            {
                
                if(skip_inactive && !cornerActive(iTri)[ip])
                    return LARGE_TRIMESH;

                bary[ip] = 1.; bary[iNode] = bary[ipp] = 0.;
                return calcDist(p,node_(iTri)[ip],delta);
            }
        }
    }

    if(skip_inactive && !cornerActive(iTri)[iNode])
        return LARGE_TRIMESH;

    bary[iNode] = 1.; bary[ip] = bary[ipp] = 0.;
    return calcDist(p,node_(iTri)[iNode],delta);
}

/* ---------------------------------------------------------------------- */

inline double TriMesh::resolveFaceContactBary(const int iTri, const double *const p, const double *const node0ToSphereCenter, double *const delta)
{
    double *surfNorm = SurfaceMeshBase::surfaceNorm(iTri);

    double dNorm = vectorDot3D(surfNorm,node0ToSphereCenter);

    double csPlane[3], tmp[3];
    vectorScalarMult3D(surfNorm,dNorm,tmp);
    vectorSubtract3D(p,tmp,csPlane);

    return calcDist(p,csPlane,delta);
}

/* ---------------------------------------------------------------------- */

inline bool TriMesh::resolveTriSphereNeighbuild(int nTri, double rSphere,
  double *cSphere, double treshold)
{
    
    double maxDist = rSphere + treshold;

    double dNorm = fabs( calcDistToPlane(cSphere,
        SurfaceMeshBase::center_(nTri),SurfaceMeshBase::surfaceNorm(nTri)) );
    if (dNorm > maxDist)
        return false;

    double **node = MultiNodeMesh<3>::node_(nTri);
    double **edgeNorm = SurfaceMeshBase::edgeNorm(nTri);

    // d_para^2 + d_norm^2 > maxDist^2 --> return false
    double dParaMax = maxDist*maxDist;// - dNorm*dNorm;

    for(int i=0;i<3;i++)
    {
        double d = calcDistToPlane(cSphere,node[i],edgeNorm[i]);
        if(d>0 && d*d > dParaMax)
            return false;
    }

    /*
    for(int i=0;i<3;i++)
      if(calcDist(cSphere,node_[i]) > maxDist)
        return false;
    */
    return true;
}

/* ---------------------------------------------------------------------- */

inline double TriMesh::calcDist(const double *const cs, const double *const closestPoint, double *const delta)
{
    vectorSubtract3D(closestPoint,cs,delta);
    return pointDistance(cs,closestPoint);
}

inline double TriMesh::calcDistToPlane(const double *const p, const double *const pPlane, const double *const nPlane)
{
    double v[3];
    vectorSubtract3D(p,pPlane,v);
    // normal distance of sphere center_ to plane
    return vectorDot3D(nPlane,v);
}

/* ----------------------------------------------------------------------
 calculate area of a triangle
------------------------------------------------------------------------- */

inline double TriMesh::calcArea(const int n)
{
    double *vecTmp3 = new double[3];

    vectorCross3D(SurfaceMeshBase::edgeVec(n)[0],
                  SurfaceMeshBase::edgeVec(n)[1],vecTmp3);

    // edgevecs are normalized so have to multiply with their lengths
    double area = 0.5*vectorMag3D(vecTmp3) * edgeLen(n)[0] * edgeLen(n)[1];
    delete[] vecTmp3;
    return area;
}

/* ----------------------------------------------------------------------
 check if point in triangle, within round-off
 from http://www.blackpawn.com/texts/pointinpoly/default.html
------------------------------------------------------------------------- */

inline bool TriMesh::isInElement(double *pos,int i)
{
    double v0[3],v1[3],v2[3];
    double dot00,dot01,dot02,dot11,dot12,invDenom,u,v;
    double ***node = node_.begin();

    vectorSubtract3D(node[i][2], node[i][0], v0);
    vectorSubtract3D(node[i][1], node[i][0], v1);
    vectorSubtract3D(pos,        node[i][0], v2);

    dot00 = vectorDot3D(v0, v0);
    dot01 = vectorDot3D(v0, v1);
    dot02 = vectorDot3D(v0, v2);
    dot11 = vectorDot3D(v1, v1);
    dot12 = vectorDot3D(v1, v2);

    invDenom = 1. / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    if((u > -SMALL_TRIMESH) && (v > -SMALL_TRIMESH) && (u + v < 1.+SMALL_TRIMESH))
        return true;
    else
        return false;
}

/* ----------------------------------------------------------------------
 generates a random point on the surface that lies within my subbox
------------------------------------------------------------------------- */

inline int TriMesh::generateRandomSubbox(double *pos)
{
    int index;
    do
    {
        index = generateRandomOwnedGhost(pos);
    }
    while (!domain->is_in_subdomain(pos));

    return index;
}

/* ----------------------------------------------------------------------
 generates a random point on the surface that lies on an owned or ghost element
------------------------------------------------------------------------- */

inline int TriMesh::generateRandomOwnedGhost(double *pos)
{
    double u,v, bary_0,bary_1,bary_2;
    double ***node = node_.begin();
    int nTri = sizeLocal() + sizeGhost();

    // step 1 - choose triangle
    int chosen = randomOwnedGhostElement();
    
    if(chosen >= nTri || chosen < 0)
    {
        
        error->one(FLERR,"TriMesh::generate_random error");
        return -1;
    }

    // step 2 - random bary coords
    
    do
    {
        u = random_->uniform();
        v = random_->uniform();
    }
    while( u+v > 1);

    bary_0 = 1. - u - v;
    bary_1 = v;
    bary_2 = u;

    pos[0] = bary_0 * node[chosen][0][0] + bary_1 * node[chosen][1][0] + bary_2 * node[chosen][2][0];
    pos[1] = bary_0 * node[chosen][0][1] + bary_1 * node[chosen][1][1] + bary_2 * node[chosen][2][1];
    pos[2] = bary_0 * node[chosen][0][2] + bary_1 * node[chosen][1][2] + bary_2 * node[chosen][2][2];

    return chosen;
}

/* ----------------------------------------------------------------------
 not implemented in this class
------------------------------------------------------------------------- */

inline int TriMesh::generateRandomOwnedGhostWithin(double *pos,double delta)
{
    UNUSED(pos);
    UNUSED(delta);
    error->one(FLERR,"internal error");
    return 0;
}

inline double TriMesh::calculateWeightedProperties(const double *const x, const double *const *const tri_nodes, const double *const tri_normal, const double *const bary, const double *const x_to_node0, const double radius, double *const delta)
{
    
    const double delta_exponent = 1.0;

    double delta_face, delta_edge[3], delta_vertex[3];
    int contact_edge[3];
    bool contact_vertex[3];
    double normal_face[3], normal_edge[3][3], normal_vertex[3][3], vertex_angle[3];
    const bool contact_face = weightedFaceContact(tri_normal, bary, x_to_node0, radius, delta_face, normal_face);
    weightedEdgeContact(x, tri_nodes, radius, contact_edge, delta_edge, normal_edge);
    weightedVertexContact(x, tri_nodes, radius, contact_vertex, delta_vertex, normal_vertex, vertex_angle);
    
    double prev_weight[3] = {0., 0., 0.};
    delta[0] = 0.;
    delta[1] = 0.;
    delta[2] = 0.;
    for (int j = 0; j < 3; j++)
    {
        if (contact_edge[j] != 0)
        {
            const double weight = contact_edge[j] == 1 ? -0.5 : 0.5;
            vectorAddMultiple3D(delta, weight*pow(delta_edge[j], delta_exponent), normal_edge[j], delta);
            if (contact_edge[j] == 1)
            {
                prev_weight[j] -= 0.5;
                prev_weight[(j+1)%3] -= 0.5;
            }
            else
            {
                prev_weight[j] += 0.5;
                prev_weight[(j+1)%3] += 0.5;
            }
        }
    }
    if (contact_face)
    {
        for (int j = 0; j < 3; j++)
            prev_weight[j] += 1.0;
        
        vectorAddMultiple3D(delta, pow(delta_face, delta_exponent), normal_face, delta);
    }
    for (int j = 0; j < 3; j++)
    {
        if (contact_vertex[j])
        {
            const double target_weight = vertex_angle[j]/M_PI*0.5;
            // adjust weights for small points
            int sum_contact_edge_outside = 0;
            if (contact_vertex[(j+1)%3] && contact_edge[j] == 0 && delta_vertex[j] < delta_vertex[(j+1)%3])
            {
                if (bary[(j+2)%3] > 0)
                {
                    prev_weight[j] -= 0.5;
                    sum_contact_edge_outside += 1;
                }
                else
                {
                    prev_weight[j] += 0.5;
                    sum_contact_edge_outside -= 1;
                }
            }
            else if (contact_edge[j] != 0)
                sum_contact_edge_outside += contact_edge[j];
            if (contact_vertex[(j+2)%3] && contact_edge[(j+2)%3] == 0 && delta_vertex[j] < delta_vertex[(j+2)%3])
            {
                if (bary[(j+1)%3] > 0)
                {
                    prev_weight[j] -= 0.5;
                    sum_contact_edge_outside += 1;
                }
                else
                {
                    prev_weight[j] += 0.5;
                    sum_contact_edge_outside -= 1;
                }
            }
            else if (contact_edge[(j+2)%3] != 0)
                sum_contact_edge_outside += contact_edge[(j+2)%3];
            // sum = 1 => one inside contact + one with no contact
            // sum = 2 => both inside
            if (sum_contact_edge_outside >= 1 && !contact_face)
                prev_weight[j] += 1.0;
            double weight = target_weight - prev_weight[j];
            vectorAddMultiple3D(delta, weight*pow(delta_vertex[j], delta_exponent), normal_vertex[j], delta);
        }
    }
    const double lenDelta = vectorMag3D(delta);
    if(lenDelta < 1e-300)
        return 1.0;
    const double deltan = pow(lenDelta, 1.0/delta_exponent);
    // scale delta to be the normal * (radius - deltan).
    vectorScalarMult3D(delta, -(radius - deltan)/lenDelta);
    
    // the negative sign is important here
    return -deltan;
}

inline double TriMesh::weightedFaceContact(const double *const tri_normal, const double *const bary, const double *const x_to_node0, const double radius, double &delta, double *const normal)
{
    delta = 0.;
    vectorZeroize3D(normal);
    const double dist = vectorDot3D(tri_normal, x_to_node0);
    if (fabs(dist) > radius)
        return false;
    // get and test barycentric coords
    if (bary[0] < 0.0 || bary[1] < 0.0 || bary[2] < 0.0)
        return false;
    delta = radius - fabs(dist);
    vectorCopy3D(tri_normal, normal);
    if (dist < 0.)
        vectorScalarMult3D(normal, -1.0);
    return true;
}

inline void TriMesh::weightedEdgeContact(const double *const x, const double *const *const tri_nodes, const double radius, int *const contact, double *const delta, double normal[3][3])
{
    vectorZeroize3D(delta);
    for (int i = 0; i < 3; i++)
    {
        vectorZeroize3D(normal[i]);
        double edge[3], ax[3], bx[3];
        vectorSubtract3D(tri_nodes[(i+1)%3], tri_nodes[i], edge);
        vectorNormalize3D(edge);
        vectorSubtract3D(tri_nodes[i], x, ax);
        vectorSubtract3D(tri_nodes[(i+1)%3], x, bx);
        const double dot_a = vectorDot3D(ax, edge);
        const double dot_b = vectorDot3D(bx, edge);
        if (dot_a <= 0. && dot_b >= 0.)
        {
            const double dist = ::sqrt(vectorMag3DSquared(ax) - dot_a*dot_a);
            if (dist < radius)
            {
                delta[i] = radius - dist;
                vectorAddMultiple3D(ax, -dot_a, edge, normal[i]);
                vectorScalarMult3D(normal[i], -1.0/vectorMag3D(normal[i]));
                // check if we are on the side of the triangle face or outside
                double edge2[3];
                vectorSubtract3D(tri_nodes[(i+2)%3], tri_nodes[i], edge2);
                double scale = vectorDot3D(edge, edge2)/vectorMag3D(edge);
                vectorAddMultiple3D(edge2, -scale, edge, edge2); // this is the edge normal
                if (vectorDot3D(edge2, normal[i]) > 0)
                    contact[i] = 1;
                else
                    contact[i] = -1;
            }
            else
                contact[i] = 0;
        }
        else
            contact[i] = 0;
    }
}

inline void TriMesh::weightedVertexContact(const double *const x, const double *const *const tri_nodes, const double radius, bool *const contact, double *const delta, double normal[3][3], double *const vertex_angle)
{
    vectorZeroize3D(delta);
    for (int i = 0; i < 3; i++)
    {
        vectorZeroize3D(normal[i]);
        const double d = pointDistance(x, tri_nodes[i]);
        if (d < radius)
        {
            delta[i] = radius - d;
            contact[i] = true;
            vectorSubtract3D(x, tri_nodes[i], normal[i]);
            vectorNormalize3D(normal[i]);
            double e1[3], e2[3];
            vectorSubtract3D(tri_nodes[(i+1)%3], tri_nodes[i], e1);
            vectorSubtract3D(tri_nodes[(i+2)%3], tri_nodes[i], e2);
            const double cos_angle = vectorDot3D(e1, e2)/(vectorMag3D(e1)*vectorMag3D(e2));
            vertex_angle[i] = acos(cos_angle);
        }
        else
            contact[i] = false;
    }
}

#endif
