//
//  LightSource.h
//  RayTracing1
//
//  Created by Walker Kennedy on 3/20/14.
//  Copyright (c) 2014 ait.hu.bud.kennedy. All rights reserved.
//

#ifndef RayTracing1_LightSource_h
#define RayTracing1_LightSource_h
#include "float2.h"
#include "float3.h"
class LightSource
{
public:
    virtual float3 getRadianceAt ( float3 x )=0;
    virtual float3 getLightDirAt ( float3 x )=0;
    virtual float getDistanceFrom( float3 x )=0;
};


/*To find the radiance, we have to divide the power
 by the surface area on which it is distributed.
 The further light gets from the light source, the larger the radius of this spherical surface is.
 Thus, we need to divide by the area of the sphere with radius |s-x|. L = /4|s-x|2
 */
class DirectionalLight : public LightSource
{
    float3 dir;
    float3 radiance;
public:
    
    
    DirectionalLight(float3 d, float3 r)
    {
        dir = d.normalize();
        radiance = r;
    }
    
    virtual float3 getRadianceAt ( float3 x ) {
        return radiance;
    }
    
    virtual float3 getLightDirAt ( float3 x ) {
        return dir;
    }
    
    virtual float getDistanceFrom( float3 x ) {
        float distanceFrom = 10000000;
        return distanceFrom;
    }
};

class PointLight : public LightSource
{
    float3 position;
    float3 power;
public:
    
    PointLight(float3 pos, float3 pow)
    {
        position = pos;
        power = pow;
    }
    
    virtual float3 getRadianceAt ( float3 x )
    {
        return power * (1 / (position - x).dot(position - x));
    }
    
    float3 getLightDirAt ( float3 x )
    {
        return (position-x).normalize();
    }
    
    float getDistanceFrom( float3 x )
    {
        return (x-position).norm();
    }
    
};

#endif
