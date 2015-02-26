#ifndef GAME_H
#define GAME_H

#include <vector>
#include "Screen.h"
#include "Vec.h"

struct Sphere {
    Linear::Vec3f position;
    float radius;
    SDL_Color c;
    //float reflectivity;
};

class Game {
    private:
        Screen* scr;
        float hfov;
        const float plane_width;
        const float plane_height;
        float plane_distance;
        Linear::Vec3f orientation;
        Linear::Vec3f position;
        std::vector<Sphere> spheres;
        const int num_threads;
        bool running;

        void render_slice(int slice);
        void handle_input();
        void draw_game();
    public:
        Game(Screen* scr);
        ~Game();
        void run();
};

#endif
