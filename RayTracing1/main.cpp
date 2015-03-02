//
//  main.cpp
//  RayTracing
//
//  Created by Walker Kennedy on 3/20/14.
//  Copyright (c) 2014 ait.hu.bud.kennedy. All rights reserved.
//

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>


#include <algorithm>
#include <float.h>
#include "float2.h"
#include "float3.h"
#include "LightSource.h"
#include <vector>
#include <map>
#include "perlin.h"

// Skeletal Material class. Feel free to add methods e.g. for illumination computation (shading).
class Material
{
public:
	bool reflective;
	bool refractive;
	bool textured;
    bool glowing;
	float3 f0;			// Fresnel coefficient
	float n;			// index of refraction
	float3 kd;			// diffuse reflection coefficient
	float3 ks;			// specular reflection coefficient
	float shininess;	// specular exponent
    float3 r;           // emitted radiance
	Material()
	{
		reflective = false;
		refractive = false;
		textured = false;
        glowing = false;
		f0 = float3(0.93, 0.85, 0.4);
		n = 1;
		kd = float3(0.5, 0.5, 0.5) + kd * 0.5;
		ks = float3(1, 1, 1);
		shininess = 15;
        r = float3(0, 0, 0);
        
	}
    
    // allows us to implement different kd functions for things like texturing
    virtual float3 getKD(float2 uv)
    {
        return kd;
    }
    
    // used for 3D procedural texturing
    virtual float getKDNoise(float3 pos)
    {
        return 1;
    }
    
};


// Checkerboard material
class Checkerboard : public Material
{
    float3 getKD(float2 uv) // get KD returns 1 of 2 colors
    {
        float3 color1(0,0,1);
        float3 color2(1,0,0);
        
        bool x = (fmodf(uv.x,0.1) < 0.05); // allows for horizontal stripes
        bool y = (fmodf(uv.y, 0.1) < 0.05); // allows for vertial stripes
        if (x != y)
        {
            return color1;
        } else
        {
            return color2;
        }
    }
};


// Marble material
class Marble : public Material
{
    float getKDNoise(float3 pos)
    {
        
        Perlin perlin = *new Perlin();
        float noise = perlin.marble(pos);// use perlin marble function
        return noise; // return noise
    }
};


// Skeletal Camera class. Feel free to add custom initialization, set aspect ratio to fit viewport dimensions, or animation.
class Camera
{
	float3 eye;
    
	float3 ahead;
	float3 right;
	float3 up;
    
	float2 tanFovHalf;
public:
	float3 getEye()
	{
		return eye;
	}
	Camera()
	{
		eye = float3(0, 0, 3);
		ahead = float3(0, 0, -1);
		right = float3(1, 0, 0);
		up = float3(0, 1, 0);
		tanFovHalf = float2(1, 1);
	}
    
	float3 rayDirFromNdc(const float2 ndc)
	{
		return (ahead + right * ndc.x * tanFovHalf.x + up * ndc.y * tanFovHalf.y).normalize();
	}
};

// Ray structure.
class Ray
{
public:
    float3 origin;
    float3 dir;
    float3 normal;
    Ray(float3 o, float3 d)
    {
        origin = o;
        dir = d;
        normal = dir.normalize();
    }
};

// Hit record structure. Contains all data that describe a ray-object intersection point.
class Hit
{
public:
	Hit()
	{
		t = -1;	// means no intersection found
	}
	float t;
    float tExit;
	float3 position;
	float3 normal;
	float2 uv;
	Material* material;
};

// Object abstract base class.
class Intersectable
{
protected:
	Material* material;
public:
	Intersectable(Material* material):material(material) {}
    virtual Hit intersect(const Ray& ray)=0;
};


// Egg/ Ellipsoid class
class Egg : public Intersectable
{
    float3 center;
	float3 coeff;
public:
    Egg(const float3& center, const float3& coeff, Material* material):
    Intersectable(material),
    center(center),
    coeff(coeff)
    {
    }
    Hit intersect(const Ray& ray)
    {
        float3 diff = ray.origin - center;
        // quadratic equation coefficients- recalculated for ellipsoid
        double a = (ray.dir/coeff).dot(ray.dir/coeff);
        double b = (diff/coeff).dot(ray.dir/coeff) * 2.0;
        double c = (diff/coeff).dot(diff/coeff) - 1;
        
        double discr = b * b - 4.0 * a * c;
        if ( discr < 0 )			// line of ray does not intersect sphere
            return Hit();			// empty hit record with negative t
        double sqrt_discr = sqrt( discr );
        double t1 = (-b + sqrt_discr)/2.0/a;
        double t2 = (-b - sqrt_discr)/2.0/a;
        
        float t = (t1<t2)?t1:t2;	// take minimum
        if(t < 0)					// minimum is behind us
            t = (t1<t2)?t2:t1;		// take maximum
        if (t < 0)					// sphere not intersected
            return Hit();			// empty hit record with negative t
        
        Hit h;
        h.t = t;
        h.material = material;
        h.position = ray.origin + ray.dir * t;
        h.normal = (h.position - center)/(coeff*coeff);
        h.normal.normalize();
        h.uv = float2( h.normal.dot(float3(0, 1, 0)) * 0.5 + 0.5,  atan2(h.normal.z, h.normal.x) / (2 * M_PI) + 0.5 );
		return h;
        
    }
};

// Plane class
class Plane : public Intersectable
{
    float3 normal;
	float3 x0;
public:
    Plane(const float3& normal, const float3& x0, Material* material):
    Intersectable(material),
    normal(normal),
    x0(x0)
    {
    }
    Hit intersect(const Ray& ray)
    {
		float t = (x0.dot(normal) - ray.origin.dot(normal)) / ray.dir.dot(normal);
		if(t<0)
			return Hit();
        
		Hit h;
		h.t = t;
		h.tExit = FLT_MAX;
		h.material = material;
		h.position = ray.origin + ray.dir * h.t;
		h.normal = normal;
		h.normal.normalize();
		float3 tangent = normal.cross(float3(0, 1, 0));
		if(tangent.norm() < 0.1)
			tangent = normal.cross(float3(1, 0, 0));
		float3 binormal = normal.cross(tangent);
		h.uv = float2( (h.position-x0).dot(tangent), (h.position-x0).dot(binormal));
        
		return h;
        
    }
};


// Cube using Smit's algorithm
/* Essentially, this class calculate a cube by analyzing
 * the intersection of a ray and the possible multiple
 * hits as it hits the planes of the cube
 * We must calculate that these hits are valid
 * and then return the one closest to the ray origin
*/

class Cube : public Intersectable
{
    float3 maxCorner;
    float3 minCorner;
    
public:
    Cube(float3 min, float3 max, Material* material):
    Intersectable(material),
    maxCorner(max),
    minCorner(min)
    {
    }
    
    
    Hit intersect(const Ray& ray)
    {
        float3 normal;
        float tmin, tmax, tymin, tymax, tzmin, tzmax;
        Hit h;
        h.normal = float3(1,0,0);
        
        // First we must find tmin and tmaz
        if (ray.dir.x >= 0)
        {
            tmin = (minCorner.x - ray.origin.x) / ray.dir.x;
            tmax = (maxCorner.x - ray.origin.x) / ray.dir.x;
        } else {
            normal = float3(1,0,0);
            tmin = (maxCorner.x - ray.origin.x) / ray.dir.x;
            tmax = (minCorner.x - ray.origin.x) / ray.dir.x;
        }
        
        if (ray.dir.y >= 0)
        {
            tymin = (minCorner.y - ray.origin.y) / ray.dir.y;
            tymax = (maxCorner.y - ray.origin.y) / ray.dir.y;
        } else {
            tymin = (maxCorner.y - ray.origin.y) / ray.dir.y;
            tymax = (minCorner.y - ray.origin.y) / ray.dir.y;
        }
        
        // invalid- return blank
        if ((tmin > tymax) || (tymin > tmax))
        {
            return Hit();
        }
        
        if (tmin < tymin)
        {
            tmin = tymin;
            h.normal = float3(0,1,0);
            
        }
        
        if (tmax > tymax)
        {
            tmax = tymax;
        }
        
        if (ray.dir.z >= 0)
        {
            tzmin = (minCorner.z - ray.origin.z) / ray.dir.z;
            tzmax = (maxCorner.z - ray.origin.z) / ray.dir.z;
        } else {
            tzmin = (maxCorner.z - ray.origin.z) / ray.dir.z;
            tzmax = (minCorner.z - ray.origin.z) / ray.dir.z;
        }
		
        // invalid- return blank
        if ((tmin > tzmax) || (tmax < tzmin))
        {
            return Hit();
        }
        
        if (tmin < tzmin)
        {
            tmin = tzmin;
            h.normal = float3(0,0,1); // set new normal
        }
        
        if (tmax > tzmax)
        {
            tmax = tzmax;
            normal = float3(0,0,-1); // set new normal
        }
        
        
        
        h.t = tmin;
        h.tExit = tmax;
		h.material = material;
		h.position = ray.origin + ray.dir * h.t;
		h.normal.normalize();
        
        return h;
        
    }
    
    
    
};


// Object realization.
class Sphere : public Intersectable
{
	float3 center;
	float radius;
public:
    Sphere(const float3& center, float radius, Material* material):
    Intersectable(material),
    center(center),
    radius(radius)
    {
    }
    Hit intersect(const Ray& ray)
    {
        float3 diff = ray.origin - center;
		// quadratic equation coefficients
        double a = ray.dir.dot(ray.dir);
        double b = diff.dot(ray.dir) * 2.0;
        double c = diff.dot(diff) - radius * radius;
        
        double discr = b * b - 4.0 * a * c;
        if ( discr < 0 )			// line of ray does not intersect sphere
            return Hit();			// empty hit record with negative t
        double sqrt_discr = sqrt( discr );
        double t1 = (-b + sqrt_discr)/2.0/a;
        double t2 = (-b - sqrt_discr)/2.0/a;
        
		float t = (t1<t2)?t1:t2;	// take minimum
		if(t < 0)					// minimum is behind us
			t = (t1<t2)?t2:t1;		// take maximum
		if (t < 0)					// sphere not intersected
            return Hit();			// empty hit record with negative t
        
		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		h.normal = h.position - center;
		h.normal.normalize();
		h.uv = float2( h.normal.dot(float3(0, 1, 0)) * 0.5 + 0.5,  atan2(h.normal.z, h.normal.x) / (2 * M_PI) + 0.5 );
        
		return h;
        
    }
};

class Scene
{
	Camera camera;
    std::vector<LightSource*> lightSources;
    std::vector<Intersectable*> objects;
    std::vector<Material*> materials;
public:
	Scene()
	{
        // light sources
        lightSources.push_back(new PointLight(float3(1, 3, -2), float3(5, 5, 5)));
		lightSources.push_back(new DirectionalLight(float3(1, 0, 0), float3(1, 1, 1)));
		lightSources.push_back(new DirectionalLight(float3(0, 5, 0), float3(0.5, 0.5, 0.5)));
        lightSources.push_back(new DirectionalLight(float3(0,0,1), float3(1,1,1)));
        
        
        // Materials
		Material* lambertian			= new Material();
        lambertian->kd = float3(0.6, 0.7, 0.6);
        lambertian->ks = float3(0, 0, 0);
        
		Material* phongBlinn			= new Material();
        phongBlinn->kd = float3(0.2, 0, 0);
        phongBlinn->ks = float3(1, 1, 1);
        phongBlinn->shininess = 25;
        
		Material* lambertianPhongBlinn	= new Material();
        lambertianPhongBlinn->kd = float3(0.3, 0.3, 0);
        lambertianPhongBlinn->ks = float3(1, 1, 1);
        lambertianPhongBlinn->shininess = 25;
        
		Material* idealReflector		= new Material();
        idealReflector->kd = float3(0, 0, 0);
        idealReflector->ks = float3(0, 0, 0);
        idealReflector->reflective = true;
        idealReflector->f0 = float3(0, 1, 1);
        
        Material* checkered = new Checkerboard;
        checkered->kd = float3(0.6, 0.7, 0.6);
        checkered->ks = float3(0, 0, 0);
        checkered->shininess = 25;
        
        Material* marble = new Marble;
        marble->kd = float3(1.0, 1.0, 1.0);
        marble->ks = float3(0, 0, 0);
        marble->shininess = 10;
        marble->textured = true;
        
        Material* gold = new Material();
        gold->kd = float3(0, 0, 0);
        gold->ks = float3(0, 0, 0);
        gold->reflective = true;
        gold->f0 = float3((powf(0.17-1,2)+powf(3.1,2))/(powf(0.17+1,2)+powf(3.1,2)), // Fresnel coefficient calc.
                          (powf(0.35-1,2)+powf(2.7,2))/(powf(0.35+1,2)+powf(2.7,2)),
                          (powf(1.5-1,2)+powf(1.9,2))/(powf(1.5+1,2)+powf(1.9,2)));
        
        Material* copper = new Material();
        copper->kd = float3(0, 0, 0);
        copper->ks = float3 (0 ,0, 0);
        copper->reflective = true;
        copper->shininess = 20;
        copper->f0 = float3((powf(0.2-1,2)+powf(3.6,2))/(powf(0.2+1,2)+powf(3.6,2)),
                            (powf(1.1-1,2)+powf(2.6,2))/(powf(1.1+1,2)+powf(2.6,2)),
                            (powf(1.2-1,2)+powf(2.3,2))/((powf(1.2+1,2)+powf(2.3,2))));
        
        Material* platinum = new Material();
        platinum->kd = float3(0, 0, 0);
        platinum->ks = float3 (0 ,0, 0);
        platinum->reflective = true;
        platinum->shininess = 20;
        platinum->f0 = float3((powf(2.37-1,2)+powf(4.2,2))/(powf(2.37+1,2)+powf(4.2,2)),
                              (powf(2.06-1,2)+powf(3.58,2))/(powf(2.06+1,2)+powf(3.58,2)),
                              (powf(1.83-1,2)+powf(3.1,2))/((powf(1.83+1,2)+powf(3.1,2))));
        
        Material* glass = new Material();
        glass->kd = float3(0, 0, 0);
        glass->ks = float3(0, 0, 0);
        glass->refractive = true;
        glass->reflective = true;
        glass->n = 1.1;
        glass->f0 = float3(0,0,0);
        
        Material *glow = new Material();
        glow->kd = float3(0, 0, 0);
        glow->ks = float3(1.0, 1.0, 1.0);
        //glow->refractive = true;
        glow->glowing = true;
        glow->r = float3(2.8,0.7,0.7);
        glow->n = 1.2;
        glow->shininess = 40;
        glow->f0 = float3(0,0,0);
        
		materials.push_back(lambertian				);
		materials.push_back(phongBlinn				);
		materials.push_back(lambertianPhongBlinn	);
		materials.push_back(idealReflector			);
        materials.push_back(checkered               );
        materials.push_back(marble                  );
        materials.push_back(gold                    );
        materials.push_back(glass                   );
        materials.push_back(copper                  );
        materials.push_back(platinum                );
        materials.push_back(glow                    );
        
        // add infinite ground plane
		objects.push_back(new Plane(float3(0, 1, 0), float3(0, -3, 0),       lambertian));
        
        // add cube
        objects.push_back(new Cube(float3(2,-2,-3),float3(3,-1,-2), lambertianPhongBlinn));
        
        // add egg with 2d procedural texturing
		objects.push_back(new Egg(float3(0, -1, -2), float3(0.5,.75,.5),     checkered));
        
        // add egg with 3d procedural texturing
        objects.push_back(new Egg(float3(-4, -1, -2), float3(0.5,.75,.5),     marble));
        
        // add egg with gold fresnel reflection
        objects.push_back(new Egg(float3(-3, -1, -2), float3(0.5,.75,.5), 	 gold));
        
        // add comparison copper egg
        objects.push_back(new Egg(float3(-2, -1, -2), float3(0.5,.75,.5),    copper));
        
        // add comparison platinum egg
        objects.push_back(new Egg(float3(-1, -1, -2), float3(0.5,.75,.5),  	 platinum));
        
        // add glass egg
        objects.push_back(new Egg(float3(0, -1, 0), float3(0.65, 0.9, 0.65), glass));
        
        // add glowing egg and light source
        objects.push_back(new Egg(float3(-2, 1, -2), float3(0.5, 0.75, 0.5),  glow));
        lightSources.push_back(new PointLight(float3(-2, 1, -2), float3(8.8, 2.7, 2.7)));
        
        
        
        // Bunny
        objects.push_back(new Sphere(float3(0,1,0), 0.5, lambertian));
        objects.push_back(new Egg(float3(-0.2,1.5,0.3), float3(0.1,0.45,0.1), lambertian));
        objects.push_back(new Egg(float3(0.2,1.5,0.3), float3(0.1,0.45,0.1), lambertian));
        objects.push_back(new Sphere(float3(0,0.8,0.55), 0.08, lambertian));
        
        // Cloud
        objects.push_back(new Egg(float3(0.5,0,0), float3(0.45,0.1,0.45), marble));
        objects.push_back(new Egg(float3(0.3,0.2,0.1), float3(0.45,0.11,0.45), marble));
        objects.push_back(new Egg(float3(0.6,0.1,0), float3(0.45,0.1,0.45), marble));
        objects.push_back(new Egg(float3(0.7,0.15,0.05), float3(0.45,0.14,0.45), marble));
        objects.push_back(new Egg(float3(0.9,0.2,0.1), float3(0.45,0.16,0.45), marble));
        objects.push_back(new Egg(float3(0.3,0.1,0), float3(0.45,0.17,0.45), marble));
        objects.push_back(new Egg(float3(1.0,0.21,0.04), float3(0.45,0.15,0.45), marble));
        objects.push_back(new Egg(float3(1.2,0.1,0), float3(0.45,0.14,0.45), marble));
        objects.push_back(new Egg(float3(1.23,0,0), float3(0.45,0.12,0.45), marble));

        
	}
	~Scene()
	{
        // UNCOMMENT THESE WHEN APPROPRIATE
        for (std::vector<LightSource*>::iterator iLightSource = lightSources.begin(); iLightSource != lightSources.end(); ++iLightSource)
        {
            
            delete *iLightSource;
        }
        for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
            delete *iMaterial;
        for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
            delete *iObject;
	}
    
public:
	Camera& getCamera()
	{
		return camera;
	}
    
	float3 trace(const Ray& ray, int depth)
	{
        
        float3 bestLight(0,0,0);
        if (depth > 5) // base case for reflective/refractive recursion
        {
            return float3(0,0,0);
        }
        
        Hit hit = bestHit(ray);
        if(hit.t < 0)	// no intersection found
            // background color
            return float3(0,0,0);
        else {
            float cosa = fabs(ray.normal.dot(hit.normal));
            
            // is reflective
            if (hit.material->reflective)
            {
                // calculate reflected ray
                Ray reflectedRay(hit.position + hit.normal * 0.01, ray.dir - hit.normal * ray.dir.dot(hit.normal) * 2);
                // add light calculated from reflected ray and multiply by Fresnel
                bestLight += trace(reflectedRay, depth+1)* (hit.material->f0 + ((float3(1,1,1) - hit.material->f0) * pow( 1-cosa, 5 )));
                
            }
            
            // material is refractive
            if (hit.material->refractive)
            {
                // initialize both refraction index to 1
                float n1 = 1;
                float n2 = 1;
                float3 origin;
                
                //entering object
                if (hit.normal.dot(ray.dir) < 0)
                {
                    n2 = hit.material->n; // n2 has material index of refraction
                    origin = hit.position - hit.normal*0.0001;
                }
                // exit
                else {
                    n1 = hit.material->n; // n1 has material index of refraction
                    origin = hit.position + hit.normal*0.0001;
                }
                
                
                float c = -ray.normal.dot(hit.normal);
                if(c < 0) {
                    c = -c;
                }
                float r = n1/n2;
                
                // refraction direction formula from slides
                float3 direction = ray.dir * r + hit.normal * ( r*c - sqrtf( 1 - r*r * (1 - c*c) ));
                
                Ray refractedRay(origin, direction);
                // trace refracted ray
                bestLight += trace(refractedRay, depth+1) * (float3(1,1,1)- (hit.material->f0 + ((float3(1,1,1) - hit.material->f0) * pow( 1-cosa, 5 ))));
            }
            
            // iterate through light sources
            for (int i = 0; i < lightSources.size(); i++) {
                
                // calculate possible shadow
                Ray shadow = Ray(hit.position + hit.normal*0.01, lightSources.at(i)->getLightDirAt(hit.position));
                
                // see if we find an object or obstruction to light
                Hit obstruction = bestHit(shadow);

                // if no obstruction (or object is glowing)- calculate light
                if (obstruction.t < 0 || obstruction.t > lightSources[i]->getDistanceFrom(hit.position) || obstruction.material->glowing) {
                    float cosTheta = hit.normal.dot( lightSources[i]->getLightDirAt(hit.position) );
                    if(cosTheta < 0)
                        cosTheta = 0;
                    float3 half = lightSources[i]->getLightDirAt(hit.position) - ray.dir;
                    half.normalize();
                    float cosDelta = half.dot(hit.normal);
                    if(cosDelta < 0)
                        cosDelta = 0;
                    
                    // check to see is material is textured (this is for 3d)
                    if (hit.material->textured){
                        bestLight += lightSources[i]->getRadianceAt(hit.position) * hit.material->getKDNoise(hit.position) *
                        (hit.material->getKD(hit.uv) * hit.material->getKDNoise(hit.position)* cosTheta  + hit.material->ks * (1 - hit.material->getKDNoise(hit.position)) * pow(cosDelta, hit.material->shininess) );
                    } else {

                        // if material is glowing, return our glowing color
                        if (hit.material->glowing)
                        {
                            bestLight = hit.material->r;
                        } else {
                            bestLight += lightSources[i]->getRadianceAt(hit.position)  *
                            (hit.material->getKD(hit.uv) * cosTheta  + hit.material->ks * pow(cosDelta, hit.material->shininess) );
                        }
                        
                    }
                }
            }
        }
        return bestLight;
	}
    
    // iterate through objects to find best hit
    Hit bestHit(Ray ray)
    {
        Hit best;
        
        for (int i = 0; i < objects.size(); i++)
        {
            Hit hit = objects.at(i)->intersect(ray);
            if (best.t < 0 || (hit.t > 0 && hit.t < best.t))
            {
                best = hit;
            }
        }
        
        return best;
    }
};

// YOU DO NOT NEED TO CHANGE ANYTHING BELOW THIS LINE
////////////////////////////////////////////////////////////////////////////////////////////////////////

// global application data

// screen resolution
const int screenWidth = 600;
const int screenHeight = 600;

// image to be computed by ray tracing
float3 image[screenWidth*screenHeight];

//scene object
Scene scene;

// Computes every 64th scanline of the image calling scene.trace() for every pixel. Returns true if scanlines are left uncomputed.
bool computeImage()
{
	static unsigned int iPart = 0;
    
	if(iPart >= 64)
		return false;
    for(int j = iPart; j < screenHeight; j+=64)
	{
        for(int i = 0; i < screenWidth; i++)
		{
			float3 pixelColor = float3(0, 0, 0);
			float2 ndcPixelCentre( (2.0 * i - screenWidth) / screenWidth, (2.0 * j - screenHeight) / screenHeight );
            
			Camera& camera = scene.getCamera();
			Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcPixelCentre));
			
			image[j*screenWidth + i] = scene.trace(ray, 0);
		}
	}
	iPart++;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL starts here. In the ray tracing example, OpenGL just outputs the image computed to the array.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#include <windows.h>
#endif // Win32 platform

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/
#include <GLUT/glut.h>

// Displays the image.
void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear screen
    
	// Compute some of image. If true was returned, the image is not complete yet, so redraw the window, computing some more of it.
	if(computeImage())
		glutPostRedisplay();
	// Output the image.
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);
    
    glutSwapBuffers(); // drawing finished
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);						// initialize GLUT
    glutInitWindowSize(600, 600);				// startup window size
    glutInitWindowPosition(100, 100);           // where to put window on screen
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    // 8 bit R,G,B,A + double buffer + depth buffer
    
    glutCreateWindow("Ray caster");				// application window is created and displayed
    
    glViewport(0, 0, screenWidth, screenHeight);
    
    glutDisplayFunc(onDisplay);					// register callback
    
    glutMainLoop();								// launch event handling loop
    
    return 0;
}



