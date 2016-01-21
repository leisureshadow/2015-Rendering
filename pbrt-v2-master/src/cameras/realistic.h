
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

using namespace std;

struct Len{
	float radius, axpos, N, aperture, center;

	Len(float radius, float axpos, float N, float aperture, float center): radius(radius), axpos(axpos), N(N), aperture(aperture), center(center) {
		if(N == 0.f)
			this->N = 1.f;
	}

	bool Intersect(const Ray *ray, Point &point) const;
};

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *ray) const;
  
private:
	Transform RasterToCamera;
	float hither, yon, filmdistance, weightCoeff;
	vector<Len> lens;
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H