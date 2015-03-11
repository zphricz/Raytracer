#ifndef SCENE_COMPONENTS_H
#define SCENE_COMPONENTS_H

// Geometry
struct Sphere {
    Linear::Vec3f position;
    float radius;
    float reflectivity;
    float spec_power;
    SDL_Color c;
};

struct Quad {
    Linear::Vec3f p1;
    Linear::Vec3f p2;
    Linear::Vec3f p3;
    Linear::Vec3f p4;
    Linear::Vec3f normal;
    float reflectivity;
    float spec_power;
    SDL_Color c;
};

// Lights are modelled as completely transparent spheres that emit light in
// all directions from all points within the sphere
struct Light {
    Linear::Vec3f position;
    float intensity;
    float radius;
    //SDL_Color c;
};

#endif
