
#include "stdafx.h"
#include "core/sampler.h"
#include "cameras/realistic.h"
#include "core/montecarlo.h"
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                 float hither, float yon, 
                 float sopen, float sclose, 
                 float filmdistance, float aperture_diameter, string specfile, 
                 float filmdiag, Film *f)
    : Camera(cam2world, sopen, sclose, f) {


    //read spec file
    float z = 0;
    ifstream inputFile(specfile);
    if (!inputFile.is_open()) {
        Severe("Cannot open lens spec file!\n");
    }
    char line[100];
    while(inputFile.getline(line, 100)){
        if(line[0] == '#')
            continue;
        istringstream iss(line);
        float radius, axpos, N, aperture, center;
        iss >> radius >> axpos >> N >> aperture;
        center = z - radius;
        z -= axpos;
        lens.push_back(Len(radius, axpos, N, aperture, center));
    }

    this->hither = hither;
    this->yon = yon;
    this->filmdistance = filmdistance;

    //calculate Transform RasterToCamera
    float rasterFilmRatio = filmdiag / sqrt(f->xResolution*f->xResolution + f->yResolution*f->yResolution);
    RasterToCamera = Scale(-rasterFilmRatio, rasterFilmRatio, 1.f) * Translate(Vector(-f->xResolution/2.f, -f->yResolution/2.f, z-filmdistance));

    float lastLenAperture = lens.back().aperture / 2.f;
    weightCoeff = lastLenAperture*lastLenAperture*M_PI / powf(fabs(filmdistance), 2.f);

}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
  
  // use sample->imageX and sample->imageY to get raster-space coordinates
  // of the sample point on the film.
  // use sample->lensU and sample->lensV to get a sample position on the lens

    Point Pras(sample.imageX, sample.imageY, 0);
    Point Pcamera;
    RasterToCamera(Pras, &Pcamera);
    float x, y, z;
    UniformSampleDisk(sample.lensU, sample.lensV, &x, &y);
    float lastLenAperture = lens.back().aperture / 2.f;
    x *= lastLenAperture;
    y *= lastLenAperture;
    float lastLenRadius = lens.back().radius;
    float zdiff = sqrt(lastLenRadius*lastLenRadius - x*x - y*y);
    z = lens.back().center + ((lastLenRadius >= 0.f) ? zdiff : -zdiff);
    *ray = Ray(Pcamera, Normalize(Point(x, y, z) - Pcamera), hither, yon);
    Point hit;

    for(int i = lens.size()-1 ; i >= 0; i--) {

        if (lens[i].radius == 0.f) {
            float hitt = (lens[i].center - ray->o.z) / ray->d.z;
            hit = (*ray)(hitt);
        }
        else if (i == lens.size() - 1) {
            hit = Point(x, y, z);
        }
        else {
            if (!lens[i].Intersect(ray, hit))
                return 0.f;
        }
        if (sqrt(hit.x * hit.x + hit.y * hit.y) > (lens[i].aperture / 2.f))
            return 0.f;

        // Stop do not need to calculate Snell's Law
        if (lens[i].radius == 0.f)
            continue;

        // snell by Heckbert's method
        float eta = lens[i].N / ((i == 0)? 1 : lens[i-1].N );
        Point center = Point(0.f, 0.f, lens[i].center);
        Vector N = lens[i].radius > 0.f ? Normalize(center - hit) : Normalize(hit - center);
        float c1 = Dot(-ray->d, N);
        float c2 = 1 - eta*eta*(1 - c1*c1);
        if(c2 <= 0.f)
            return 0.f;
        c2 = sqrt(c2);
        Vector T = eta*ray->d + (eta*c1 - c2)*N;
 
        ray->o = hit;
        ray->d = Normalize(T);
    }

    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    CameraToWorld(*ray, ray);

    float coscine = Normalize(ray->o - Pcamera).z;

    return weightCoeff*powf(coscine, 4.f);  
}

bool Len::Intersect(const Ray *ray, Point &point) const {

    // Compute quadratic sphere coefficients
    Vector D = ray->o - Point(0.f, 0.f, center);
    float A = ray->d.x*ray->d.x + ray->d.y*ray->d.y + ray->d.z*ray->d.z;
    float B = 2.f * Dot(ray->d, D);
    float C = Dot(D, D) - radius*radius;

    // Solve quadratic equation for _t_ values
    float t0, t1;
    if (!Quadratic(A, B, C, &t0, &t1))
        return false;

    float phit = (radius >= 0.f)? t1 : t0;
    point = (*ray)(phit);

    return true;

}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
    // Extract common camera parameters from \use{ParamSet}
    float hither = params.FindOneFloat("hither", -1);
    float yon = params.FindOneFloat("yon", -1);
    float shutteropen = params.FindOneFloat("shutteropen", -1);
    float shutterclose = params.FindOneFloat("shutterclose", -1);

    // Realistic camera-specific parameters
    string specfile = params.FindOneString("specfile", "");
    float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
    float fstop = params.FindOneFloat("aperture_diameter", 1.0);    
    float filmdiag = params.FindOneFloat("filmdiag", 35.0);

    Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
        shutterclose != -1 && filmdistance!= -1);
    if (specfile == "") {
        Severe( "No lens spec file supplied!\n" );
    }
    return new RealisticCamera(cam2world, hither, yon,
                   shutteropen, shutterclose, filmdistance, fstop, 
                   specfile, filmdiag, film);
}
