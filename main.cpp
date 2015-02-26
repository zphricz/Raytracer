#include <string>
#include <iostream>
#include "Screen.h"
#include "Game.h"
#include "Matrix.h"

using namespace std;

#include <random>
#include <chrono>
using namespace Linear;
using namespace chrono;
static const float PI = M_PI;
static const float TWO_PI = M_PI * 2;

static void error(char * name) {
    printf("Usage: %s [Screen_x Screen_y]\n", name);
    exit(1);
}

int main(int argc, char* argv[])
{
    SDL_DisplayMode display;
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        exit(1);
    }
    SDL_GetCurrentDisplayMode(0, &display);
    atexit(SDL_Quit);

    int screen_width = display.w;
    int screen_height = display.h;

    switch (argc) {
    case 3:
        screen_width = atoi(argv[1]);
        screen_height = atoi(argv[2]);
        if (screen_width <= 0 || screen_height <= 0) {
            error(argv[0]);
        }
    case 1:
        break;
    default:
        error(argv[0]);
    }

    bool full_screen;
    if (screen_width == display.w && screen_height == display.h) {
        full_screen = true;
    } else {
        full_screen = false;
    }

    Screen scr(screen_width, screen_height, "Raytracer", full_screen, true,
               false);
    Game g(&scr);
    g.run();
    return 0;
}

