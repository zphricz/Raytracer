#ifndef GAME_H
#define GAME_H

#include <vector>
#include "Screen.h"
#include "Vec.h"

struct Sphere {
    Linear::Vec3f position;
    float radius;
    float reflectivity;
    SDL_Color c;
};

struct Light {
    Linear::Vec3f position;
    float intensity;
    //SDL_Color c;
};

class Game {
    private:
        Screen* scr;
        const float plane_width;
        const float plane_height;
        float plane_distance;
        Linear::Vec3f position;
        Linear::Vec3f dx;
        Linear::Vec3f dy;
        Linear::Vec3f orientation;
        Linear::Vec3f top_left;
        float pitch;
        float yaw;
        std::vector<Sphere> spheres;
        std::vector<Light> lights;
        const int num_threads;
        bool running;
        bool mouse_enabled;

        void render_slice(int slice);
        void handle_input();
        void draw_game();
        void calculate_vectors();

        float trace(Linear::Vec3f origin, Linear::Vec3f direction, int& sphere_index, float cutoff = 0.001);
    public:
        Game(Screen* scr);
        ~Game();
        void run();
};

#endif
