
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "paramset.h"

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    z = new float[nx*ny];
    point = new Point[nx*ny];
    normal = new Normal[nx * ny];
    memcpy(z, zs, nx*ny*sizeof(float));

    nVoxels[0] = x-1;
    nVoxels[1] = y-1;
    bounds = ObjectBound();
    for (int axis = 0; axis < 2; ++axis) {
        width[axis] = 1.f / nVoxels[axis];
        invWidth[axis] = (width[axis] == 0.f) ? 0.f : 1.f / width[axis];
    }

    int index = 0;
    for(int _y = 0; _y < y; _y++){
        for(int _x = 0; _x < x; _x++){
            point[index] = Point(_x * width[0], _y * width[1], z[index]);
            normal[index] = Normal(0, 0, 0);
            index++;
        }
    }

    for (int _y = 0; _y < nVoxels[1]; _y++) {
        for (int _x = 0; _x < nVoxels[0]; _x++) {
            int index1 = _x + _y*nx;
            int index2 = _x+1 + _y*nx;
            int index3 = _x+1 + (_y+1)*nx;
            int index4 = _x + (_y+1)*nx;

            const Point &p1 = point[index1];
            const Point &p2 = point[index2];
            const Point &p3 = point[index3];
            const Point &p4 = point[index4];
            
            Normal &normal1 = normal[index1];
            Normal &normal2 = normal[index2];
            Normal &normal3 = normal[index3];
            Normal &normal4 = normal[index4];
            Normal normal;

            normal = Normal(Normalize(Cross(p2 - p1, p3 - p1)));
            normal1 += normal;
            normal2 += normal;
            normal3 += normal;

            normal = Normal(Normalize(Cross(p3 - p1, p4 - p1)));
            normal1 += normal;
            normal3 += normal;
            normal4 += normal;
        }
    }

    for (int i = 0; i < nx*ny; i++) {
        normal[i] = Normalize(normal[i]);
    }
}


Heightfield2::~Heightfield2() {
    delete[] z;
    delete[] point;
}


BBox Heightfield2::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}

bool Heightfield2::TriangleIntersect(const Ray &ray, float *tHit, float *rayEpsilon, DifferentialGeometry *dg, int index[]) const{
    PBRT_RAY_TRIANGLE_INTERSECTION_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));

    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1o = point[index[0]];
    const Point &p2o = point[index[1]];
    const Point &p3o = point[index[2]];
    const Point &p1 = (*ObjectToWorld)(p1o);
    const Point &p2 = (*ObjectToWorld)(p2o);
    const Point &p3 = (*ObjectToWorld)(p3o);
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;

    // Compute triangle partial derivatives
    Vector dpdu, dpdv;

    // Compute deltas for triangle partial derivatives
    float du1 = p1o.x - p3o.x;
    float du2 = p2o.x - p3o.x;
    float dv1 = p1o.y - p3o.y;
    float dv2 = p2o.y - p3o.y;

    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

    // Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*p1o.x + b1*p2o.x + b2*p3o.x;
    float tv = b0*p1o.y + b1*p2o.y + b2*p3o.y;

    // Fill in _DifferentialGeometry_ from triangle hit
    *dg = DifferentialGeometry(ray(t), dpdu, dpdv,
                               Normal(0,0,0), Normal(0,0,0),
                               tu, tv, this);
    *tHit = t;
    *rayEpsilon = 1e-3f * *tHit;
    ray.maxt = t;
    PBRT_RAY_TRIANGLE_INTERSECTION_HIT(const_cast<Ray *>(&ray), t);
    return true;
}

bool Heightfield2::TriangleIntersectP(const Ray &ray, int index[]) const{
    PBRT_RAY_TRIANGLE_INTERSECTIONP_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));
    // Get triangle vertices in _p1_, _p2_, and _p3_
    const Point &p1o = point[index[0]];
    const Point &p2o = point[index[1]];
    const Point &p3o = point[index[2]];
    const Point &p1 = (*ObjectToWorld)(p1o);
    const Point &p2 = (*ObjectToWorld)(p2o);
    const Point &p3 = (*ObjectToWorld)(p3o);
    Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;

    // Compute first barycentric coordinate
    Vector d = ray.o - p1;
    float b1 = Dot(d, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;

    // Compute second barycentric coordinate
    Vector s2 = Cross(d, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
    
    PBRT_RAY_TRIANGLE_INTERSECTIONP_HIT(const_cast<Ray *>(&ray), t);

    ray.maxt = t;
    return true;
}

bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const{
    Ray ray;
    (*WorldToObject)(r, &ray);

    float rayT;
    if(bounds.Inside(ray(ray.mint)))
        rayT = ray.mint;
    else if(!bounds.IntersectP(ray, &rayT))
        return false;
    Point gridIntersect = ray(rayT);
    
    // Set up 2D DDA for ray
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = width[axis] / ray.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = -width[axis] / ray.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk ray through voxel grid
    bool hitSomething = false;
    int index[3];
    for (;;) {

        // Check for intersection in current voxel and advance to next
        index[0] = Pos[0] + Pos[1]*nx;
        index[1] = Pos[0]+1 + Pos[1]*nx;
        index[2] = Pos[0]+1 + (Pos[1]+1)*nx;
        hitSomething |= TriangleIntersect(r, tHit, rayEpsilon, dg, index);
        
        index[1] = Pos[0]+1 + (Pos[1]+1)*nx;
        index[2] = Pos[0] + (Pos[1]+1)*nx;
        hitSomething |= TriangleIntersect(r, tHit, rayEpsilon, dg, index);
        if(hitSomething) return hitSomething;

        // Advance to next voxel
        // Find _stepAxis_ for stepping to next voxel
        int stepAxis = (NextCrossingT[0] < NextCrossingT[1])? 0 : 1;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
    return hitSomething;
}

bool Heightfield2::IntersectP(const Ray &r) const{
    Ray ray;
    (*WorldToObject)(r, &ray);

    float rayT;
    if(bounds.Inside(ray(ray.mint)))
        rayT = ray.mint;
    else if(!bounds.IntersectP(ray, &rayT))
        return false;
    Point gridIntersect = ray(rayT);
    
    // Set up 2D DDA for ray
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = posToVoxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = width[axis] / ray.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = -width[axis] / ray.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }

    // Walk ray through voxel grid
    bool hitSomething = false;
    int index[3];
    for (;;) {

        // Check for intersection in current voxel and advance to next
        index[0] = Pos[0] + Pos[1]*nx;
        index[1] = Pos[0]+1 + Pos[1]*nx;
        index[2] = Pos[0]+1 + (Pos[1]+1)*nx;
        hitSomething |= TriangleIntersectP(r, index);
        if(hitSomething) return hitSomething;

        index[1] = Pos[0]+1 + (Pos[1]+1)*nx;
        index[2] = Pos[0] + (Pos[1]+1)*nx;
        hitSomething |= TriangleIntersectP(r, index);
        if(hitSomething) return hitSomething;

        // Advance to next voxel
        // Find _stepAxis_ for stepping to next voxel
        int stepAxis = (NextCrossingT[0] < NextCrossingT[1])? 0 : 1;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
    return hitSomething;
}


void Heightfield2::GetShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry *dgShading) const {
    // Initialize _Triangle_ shading geometry with _n_ and _s_

    // *dgShading = dg;
    // return;

    Point p = Point(dg.u, dg.v, 0);
    int x = posToVoxel(p, 0);
    int y = posToVoxel(p, 1);
    int index1 = x + y*nx;
    int index2, index3;

    const Point &p1o = point[index1];
    if((dg.u-p1o.x)*width[1] > (dg.v-p1o.y)*width[0]){
        index2 = x+1 + y*nx;
        index3 = x+1 + (y+1)*nx;
    }else{
        index2 = x+1 + (y+1)*nx;
        index3 = x + (y+1)*nx;
    }
    const Point &p2o = point[index2];
    const Point &p3o = point[index3];
    const Normal &normal0 = normal[index1];
    const Normal &normal1 = normal[index2];
    const Normal &normal2 = normal[index3];

    // Compute barycentric coordinates for point
    float b[3];

    // Initialize _A_ and _C_ matrices for barycentrics
    float A[2][2] =
        { { p2o.x - p1o.x, p3o.x - p1o.x },
          { p2o.y - p1o.y, p3o.y - p1o.y } };
    float C[2] = { dg.u - p1o.x, dg.v - p1o.y };
    if (!SolveLinearSystem2x2(A, C, &b[1], &b[2])) {
        // Handle degenerate parametric mapping
        b[0] = b[1] = b[2] = 1.f/3.f;
    }
    else
        b[0] = 1.f - b[1] - b[2];

    // Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
    Normal ns;
    Vector ss, ts;
    ns = Normalize(obj2world(b[0] * normal0 + b[1] * normal1 + b[2] * normal2));
    ss = Normalize(dg.dpdu);
    
    ts = Cross(ss, ns);
    if (ts.LengthSquared() > 0.f) {
        ts = Normalize(ts);
        ss = Cross(ts, ns);
    }
    else
        CoordinateSystem((Vector)ns, &ss, &ts);
    Normal dndu, dndv;

    // Compute $\dndu$ and $\dndv$ for triangle shading geometry
    // Compute deltas for triangle partial derivatives of normal
    float du1 = p1o.x - p3o.x;
    float du2 = p2o.x - p3o.x;
    float dv1 = p1o.y - p3o.y;
    float dv2 = p2o.y - p3o.y;
    Normal dn1 = normal0 - normal2;
    Normal dn2 = normal1 - normal2;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f)
        dndu = dndv = Normal(0,0,0);
    else {
        float invdet = 1.f / determinant;
        dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
        dndv = (-du2 * dn1 + du1 * dn2) * invdet;
    }

    *dgShading = DifferentialGeometry(dg.p, ss, ts, obj2world(dndu), obj2world(dndv), dg.u, dg.v, dg.shape);
    dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
    dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
    dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
}

Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}


