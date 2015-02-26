#include <iostream>
#include <vector>
#include <future>
#include <algorithm>
#include <limits>
#include <random>
#include "Game.h"

using namespace std;
using namespace Linear;

constexpr float PI = M_PI;
constexpr float fov = 75.0;
constexpr float move_speed = 3.0;
constexpr float rotate_speed = PI / 2;
constexpr int samples_per_pixel = 1;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<float> dis(-0.5, 0.5);
static float jitter[2][samples_per_pixel];

void Game::handle_input() {
    static const Uint8* key_state = SDL_GetKeyboardState(NULL);
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
        case SDL_QUIT:
            running = false;
            break;
        case SDL_KEYDOWN:
            switch (event.key.keysym.sym) {
            case SDLK_ESCAPE:
                running = false;
                break;
            case SDLK_RETURN:
                scr->write_tga("image.tga");
                break;
            default:
                break;
            }
            break;
        case SDL_MOUSEMOTION: {
            break;
        }
        default:
             break;
        }
    }
    if (key_state[SDL_SCANCODE_UP]) {
        orientation.rotate_pitch(rotate_speed * scr->frame_time());
    }
    if (key_state[SDL_SCANCODE_DOWN]) {
        orientation.rotate_pitch(-rotate_speed * scr->frame_time());
    }
    if (key_state[SDL_SCANCODE_LEFT]) {
        orientation.rotate_yaw(-rotate_speed * scr->frame_time());
    }
    if (key_state[SDL_SCANCODE_RIGHT]) {
        orientation.rotate_yaw(rotate_speed * scr->frame_time());
    }
    if (key_state[SDL_SCANCODE_W]) {
        position += orientation * move_speed * scr->frame_time();
    }
    if (key_state[SDL_SCANCODE_A]) {
        Vec3f v = orientation;
        v.rotate_yaw(-PI / 2);
        position += v * move_speed * scr->frame_time();
    }
    if (key_state[SDL_SCANCODE_D]) {
        Vec3f v = orientation;
        v.rotate_yaw(PI / 2);
        position += v * move_speed * scr->frame_time();
    }
    if (key_state[SDL_SCANCODE_S]) {
        position -= orientation * move_speed * scr->frame_time();
    }
}

void Game::render_slice(int slice) {
    float pitch = orientation.pitch();
    float yaw = orientation.yaw();
    for (int i = slice * scr->height / num_threads; i < (slice + 1) * scr->height / num_threads; ++i) {
        for (int j = 0; j < scr->width; ++j) {
            float r = 0;
            float g = 0;
            float b = 0;
            for (int sample = 0; sample < samples_per_pixel; ++sample) {
                float plane_dist_y = plane_height / 2.0 - (i + 0.5 + jitter[0][sample]) * plane_height / scr->height;
                float pitch_diff = atan(plane_dist_y / plane_distance);
                float plane_dist_x = (j + 0.5 + jitter[1][sample]) * plane_width / scr->width - plane_width / 2.0;
                float yaw_diff = atan(plane_dist_x / plane_distance);
                float pitch_prime = pitch - pitch_diff;
                float yaw_prime = yaw - yaw_diff;
                Vec3f ray(yaw_prime, pitch_prime);
                
                float min_d = numeric_limits<float>::max();
                int sphere_index = -1;
                for (int s = 0; s < spheres.size(); ++s) {
                    const Sphere& sphere = spheres[s];
                    auto distance = sphere.position - position;
                    float B = ray.dot(distance);
                    float D = B * B - distance.dot(distance) + sphere.radius * sphere.radius;
                    if (D < 0.0) {
                        continue;
                    }
                    float t0 = B - sqrt(D);
                    float t1 = B + sqrt(D);
                    float d = numeric_limits<float>::max();
                    if ((t0 > 0.0)) {
                        d = t0;
                    }
                    if ((t1 > 0.0) && t1 < d) {
                        d = t1;
                    }
                    if (d < min_d) {
                        sphere_index = s;
                        min_d = d;
                    }
                }
                if (sphere_index != -1) {
                    const Sphere& sphere = spheres[sphere_index];
                    Vec3f hit = position + ray * min_d;
                    Vec3f sphere_normal = (hit - sphere.position) / sphere.radius;
                    Vec3f light_vec = position - hit;
                    light_vec.normalize();
                    float shade = sphere_normal.dot(light_vec);
                    r += sphere.c.r * shade;
                    g += sphere.c.g * shade;
                    b += sphere.c.b * shade;
                }
            }
            scr->draw_pixel(j, i, {Uint8(rint(r / samples_per_pixel)), Uint8(rint(g / samples_per_pixel)), Uint8(rint(b / samples_per_pixel))});
        }
    }
}

void Game::draw_game() {
    scr->cls();
    vector<future<void>> futures;
    futures.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        futures.push_back(async(launch::async, &Game::render_slice, this, i));
    }
    for (auto& f: futures) {
        f.get();
    }
}

Game::Game(Screen* scr) :
    scr(scr),
    hfov(fov),
    plane_width(scr->width),
    plane_height(scr->height),
    plane_distance(scr->width / (tan(fov / 2.0) * 2.0)),
    orientation(0.0, 0.0, 1.0), // looks 1.0 in z direction
    position(0.0, 0.0, 0.0),
    num_threads(4),
    running(true) {
    scr->set_recording_style("images", 5);
    //SDL_SetRelativeMouseMode(SDL_TRUE);
    //spheres.push_back({{-1.0, 0.0, 10.0}, 1.0, {255, 0, 0}});
    //spheres.push_back({{0.0, 0.0, 10.0}, 1.0, {0, 255, 0}});
    //spheres.push_back({{1.0, 0.0, 10.0}, 1.0, {0, 0, 255}});
    spheres.push_back({{-4.0/3.0, 0.0, 10.0}, 2.0/3.0, {255, 0, 0}});
    spheres.push_back({{0.0, 0.0, 10.0}, 2.0/3.0, {0, 255, 0}});
    spheres.push_back({{4.0/3.0, 0.0, 10.0}, 2.0/3.0, {0, 0, 255}});
    scr->set_color(255, 0, 0);
    jitter[0][0] = 0.0;
    jitter[1][0] = 0.0;
    for (int i = 1; i < samples_per_pixel; ++i) {
        jitter[0][i] = dis(gen);
        jitter[1][i] = dis(gen);
    }
}

Game::~Game() {
}

void Game::run() {
    while (running) {
        handle_input();
        draw_game();
        scr->commit();
        static int i = 0;
        if(i++ == 10) {
            cout << "FPS: " << scr->fps() << endl;
            i = 0;
        }
    }
}

