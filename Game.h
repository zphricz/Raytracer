#ifndef GAME_H
#define GAME_H

#include <vector>
#include <atomic>
#include "Screen.h"
#include "Vec.h"

//constexpr float SPHERE_HIT_ERR_MARGIN = 0.005;
constexpr float SPHERE_HIT_ERR_MARGIN = 0.0005;

struct Sphere {
    Linear::Vec3f position;
    float radius;
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

struct Cell {
    float r;
    float g;
    float b;
    int count;
};

class Game {
    private:
        Screen* scr;
        std::atomic<Cell *> cells;
        std::atomic<Cell *> cells_swap;
        std::atomic<uint64_t> rays_cast;
        std::atomic<bool> running;
        const float plane_width;
        const float plane_height;
        float plane_distance;
        float pitch;
        float yaw;
        Linear::Vec3f position;
        Linear::Vec3f dx;
        Linear::Vec3f dy;
        Linear::Vec3f dz;
        std::vector<Sphere> spheres;
        std::vector<Light> lights;
        uint64_t rays_at_last_swap;
        const int num_threads;
        bool mouse_enabled;

        void render_loop();
        void handle_input();
        void calculate_vectors();
        void clear_cells();

        float detect_sphere_hit(Linear::Vec3f origin, Linear::Vec3f ray,
                                int& sphere_index, float cutoff = SPHERE_HIT_ERR_MARGIN);
        void ray_trace(Linear::Vec3f origin, Linear::Vec3f ray,
                       float& r, float& g, float& b, float weight = 1.0, int depth = 0);
    public:
        Game(Screen* scr, int num_threads);
        ~Game();
        void run();
};

#endif
