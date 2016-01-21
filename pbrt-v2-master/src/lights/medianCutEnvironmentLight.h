#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_LIGHTS_MEDIANCUTENVIRONMENTLIGHT_H
#define PBRT_LIGHTS_MEDIANCUTENVIRONMENTLIGHT_H

// lights/medianCutEnvironmentLight.h*
#include "pbrt.h"
#include "light.h"
#include "texture.h"
#include "shape.h"
#include "scene.h"
#include "mipmap.h"

struct SimpleLight {
    SimpleLight(const Transform &light2world, float x, float y, RGBSpectrum _spectrum) {
      float theta = y * M_PI, phi = x * 2.f * M_PI;
      float costheta = cosf(theta), sintheta = sinf(theta);
      float sinphi = sinf(phi), cosphi = cosf(phi);
      wi = light2world(Vector(sintheta * cosphi, sintheta * sinphi, costheta));
      spectrum = _spectrum;
    }
    Vector wi;
    RGBSpectrum spectrum;
};

// MedianCutEnvironmentLight Declarations
class MedianCutEnvironmentLight : public Light {
public:
    // MedianCutEnvironmentLight Public Methods
    MedianCutEnvironmentLight(const Transform &light2world, const Spectrum &power, int ns,
                      const string &texmap);
    ~MedianCutEnvironmentLight();
    Spectrum Power(const Scene *) const;
    bool IsDeltaLight() const { return true; }
    Spectrum Le(const RayDifferential &r) const;
    Spectrum Sample_L(const Point &p, float pEpsilon, const LightSample &ls,
                      float time, Vector *wi, float *pdf, VisibilityTester *visibility) const;
    Spectrum Sample_L(const Scene *scene, const LightSample &ls, float u1, float u2,
                      float time, Ray *ray, Normal *Ns, float *pdf) const;
    float Pdf(const Point &, const Vector &) const;
    void SHProject(const Point &p, float pEpsilon, int lmax, const Scene *scene,
                   bool computeLightVis, float time, RNG &rng, Spectrum *coeffs) const;
private:
    // MedianCutEnvironmentLight Private Data
    MIPMap<RGBSpectrum> *radianceMap;
    Distribution2D *distribution;
    vector<SimpleLight> lights;
    float pdf;
};


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
                                       const ParamSet &paramSet);

#endif // PBRT_LIGHTS_MEDIANCUTENVIRONMENTLIGHT_H
