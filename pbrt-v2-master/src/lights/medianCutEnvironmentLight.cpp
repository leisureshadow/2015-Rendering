#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"
#include <queue>


// MedianCutEnvironmentLight Utility Classes
struct InfiniteAreaCube {
    // InfiniteAreaCube Public Methods
    InfiniteAreaCube(const MedianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
    : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};

struct Region {
    Region(int _x1, int _y1, int _x2, int _y2, int _n, float _energy) : x1(_x1), y1(_y1), x2(_x2), y2(_y2), n(_n), energy(_energy) {
        width = _x2 - _x1 + 1;
        height = _y2 - _y1 + 1;

        // printf("n=%3d width=%4d height=%4d left-up=(%4d, %4d) right-down=(%4d, %4d) %10lf\n", n, width, height, x1, y1, x2, y2, energy);
    }
    int n;
    int width, height;
    int x1, x2, y1, y2;
    float energy;
};

static float getAreaEnergy(float *sumAreaTable, int x1, int y1, int x2, int y2, int width){
    float energy = sumAreaTable[x2+y2*width];
    if( y1 > 0 )
        energy -= sumAreaTable[x2+(y1-1)*width];
    if( x1 > 0 )
        energy -= sumAreaTable[x1-1+y2*width];
    if( x1 > 0 && y1 > 0 )
        energy += sumAreaTable[x1-1+(y1-1)*width];
    return energy;
}


// MedianCutEnvironmentLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete distribution;
    delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
                                     const Spectrum &L, int ns, const string &texmap)
: Light(light2world, ns) {

    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);

    // Initialize sampling PDFs for infinite area light
    
    // Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    float *sumAreaTable = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta;

            if (u > 0 && v > 0)
                sumAreaTable[u+v*width] = img[u+v*width] + sumAreaTable[u-1+v*width] + sumAreaTable[u+(v-1)*width] - sumAreaTable[u-1+(v-1)*width];
            else if (u > 0 && v == 0)
                sumAreaTable[u+v*width] = img[u+v*width] + sumAreaTable[u-1+v*width];
            else if (u == 0 && v > 0)
                sumAreaTable[u+v*width] = img[u+v*width] + sumAreaTable[u+(v-1)*width];
            else
                sumAreaTable[u+v*width] = img[u+v*width];
        }
    }

    std::queue<Region> regionList;
    regionList.push(Region(0, 0, width-1, height-1, ns, sumAreaTable[width-1+(height-1)*width]));
    float leftEnergy, rightEnergy;
    float solidAngle = 2.f * M_PI * M_PI / ((float)width * (float)height);

    while(!regionList.empty()){
        Region& region = regionList.front();
        regionList.pop();

        if(region.n!=1){
            if(region.width >= region.height){
                if(region.width != 1){
                    int x;
                    for(x = region.x1; x < region.x2; x++){
                        leftEnergy = getAreaEnergy(sumAreaTable, region.x1, region.y1, x, region.y2, width);
                        rightEnergy = getAreaEnergy(sumAreaTable, x+1, region.y1, region.x2, region.y2, width);
                        if(leftEnergy > rightEnergy){
                            break;
                        }
                    }
                    if(x==region.x2) x = region.x2-1;
                    // printf("Region to split: n=%3d width=%4d height=%4d left-up=(%4d, %4d) right-down=(%4d, %4d) %10lf\n", region.n, region.width, region.height, region.x1, region.y1, region.x2, region.y2, region.energy);
                    regionList.push(Region(region.x1, region.y1, x, region.y2, region.n/2, leftEnergy));
                    regionList.push(Region(x+1, region.y1, region.x2, region.y2, region.n/2, rightEnergy));
                }
                else{
                    // printf("Region to split: n=%3d width=%4d height=%4d left-up=(%4d, %4d) right-down=(%4d, %4d) %10lf\n", region.n, region.width, region.height, region.x1, region.y1, region.x2, region.y2, region.energy);
                    regionList.push(Region(region.x1, region.y1, region.x2, region.y2, region.n/2, region.energy));
                    regionList.push(Region(region.x1, region.y1, region.x2, region.y2, region.n/2, region.energy));
                }
            }
            else{
                int y;
                for(y = region.y1; y < region.y2; y++){
                    leftEnergy = getAreaEnergy(sumAreaTable, region.x1, region.y1, region.x2, y, width);
                    rightEnergy = getAreaEnergy(sumAreaTable, region.x1, y+1, region.x2, region.y2, width);
                    if(leftEnergy > rightEnergy){
                        break;
                    }
                }
                if(y==region.y2) y = region.y2-1;
                // printf("Region to split: n=%3d width=%4d height=%4d left-up=(%4d, %4d) right-down=(%4d, %4d) %10lf\n", region.n, region.width, region.height, region.x1, region.y1, region.x2, region.y2, region.energy);
                regionList.push(Region(region.x1, region.y1, region.x2, y, region.n/2, leftEnergy));
                regionList.push(Region(region.x1, y+1, region.x2, region.y2, region.n/2, rightEnergy));
            }
        }
        else{
            RGBSpectrum spectrum;
            float x=0, y=0;
            if(region.energy != 0.f){
                for(int v=region.y1; v<=region.y2; v++){
                    float theta = M_PI * (float)(v+.5f) / (float)height;
                    float sintheta = sinf(theta);
                    for(int u=region.x1; u<=region.x2; u++){
                        x += u*img[u+v*width];
                        y += v*img[u+v*width];
                        spectrum += texels[u+v*width] * sintheta;
                    }
                }
                x /= region.energy;
                y /= region.energy;
                spectrum *= solidAngle;
            }
            else{
                x = region.x1;
                y = region.y1;
            }
            lights.push_back(SimpleLight(light2world, x/width, y/height, spectrum));
        }
    }

    pdf = 1.f/lights.size();

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);
    delete[] texels;
    delete[] sumAreaTable;
    delete[] img;
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
    Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
                                  int lmax, const Scene *scene, bool computeLightVis,
                                  float time, RNG &rng, Spectrum *coeffs) const {
    // Project _MedianCutEnvironmentLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _MedianCutEnvironmentLight_ to SH from lat-long representation
        
        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _MedianCutEnvironmentLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] * (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }
        
        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _MedianCutEnvironmentLight_ to SH from cube map sampling
        SHProjectCube(InfiniteAreaCube(this, scene, time, computeLightVis, pEpsilon), p, 200, lmax, coeffs);
    }
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
                                       const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
                                     const LightSample &ls, float time, Vector *wi, float *pdf,
                                     VisibilityTester *visibility) const {
    SimpleLight light = lights.at(Floor2Int(ls.uComponent * lights.size()));
    *wi = light.wi;
    visibility->SetRay(p, pEpsilon, *wi, time);
    Spectrum Ls = Spectrum(light.spectrum, SPECTRUM_ILLUMINANT);
    *pdf = this->pdf;
    return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
    (2.f * M_PI * M_PI * sintheta);
    return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
                                     const LightSample &ls, float u1, float u2, float time,
                                     Ray *ray, Normal *Ns, float *pdf) const {
    // Compute direction for infinite light sample ray
    
    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) 
        return Spectrum(0.f);
    
    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;
    
    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);
    
    // Compute _MedianCutEnvironmentLight_ ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
    return Ls;
}



