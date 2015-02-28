#include <iostream>
#include <vector>
#include <future>
#include <algorithm>
#include <limits>
#include <random>
#include "Game.h"

#define THREADED

using namespace std;
using namespace Linear;

constexpr float PI = M_PI;
constexpr float TWO_PI = M_PI * 2.0;

static constexpr float rad(float theta) {
    return theta * PI / 180.0;
}

static constexpr float deg(float theta) {
    return theta * 180.0 / PI;
}

static float clip_yaw(float yaw) {
    if (yaw < 0.0) {
        return yaw + TWO_PI;
    }
    if (yaw >= TWO_PI) {
        return yaw - TWO_PI;
    }
    return yaw;
}

static float clip_pitch(float pitch) {
    if (pitch > PI / 2) {
        return PI / 2;
    }
    if (pitch < -PI / 2) {
        return -PI / 2;
    }
    return pitch;
}

constexpr float fov = rad(75.0);
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
        pitch += rotate_speed * scr->frame_time();
        pitch = clip_pitch(pitch);
    }
    if (key_state[SDL_SCANCODE_DOWN]) {
        pitch -= rotate_speed * scr->frame_time();
        pitch = clip_pitch(pitch);
    }
    if (key_state[SDL_SCANCODE_LEFT]) {
        yaw -= rotate_speed * scr->frame_time();
        yaw = clip_yaw(yaw);
    }
    if (key_state[SDL_SCANCODE_RIGHT]) {
        yaw += rotate_speed * scr->frame_time();
        yaw = clip_yaw(yaw);
    }
    if (key_state[SDL_SCANCODE_W]) {
        position += Vec3f(pitch, yaw) * move_speed * scr->frame_time();
    }
    if (key_state[SDL_SCANCODE_A]) {
        position += Vec3f(0.0, yaw - PI / 2) * move_speed * scr->frame_time();
    }
    if (key_state[SDL_SCANCODE_D]) {
        position += Vec3f(0.0, yaw + PI / 2) * move_speed * scr->frame_time();
    }
    if (key_state[SDL_SCANCODE_S]) {
        position -= Vec3f(pitch, yaw) * move_speed * scr->frame_time();
    }
}

void Game::render_slice(int slice) {
    Vec3f orientation(pitch, yaw);
    Vec3f dx(0.0, yaw + PI / 2);
    Vec3f dy(pitch - PI / 2, yaw);
    Vec3f top_left = position + orientation * plane_distance -
                     dy * (plane_height / 2.0 + 0.5) -
                     dx * (plane_width / 2.0 + 0.5);
    for (int i = slice * scr->height / num_threads; i < (slice + 1) * scr->height / num_threads; ++i) {
        for (int j = 0; j < scr->width; ++j) {
            float r = 0;
            float g = 0;
            float b = 0;
            for (int sample = 0; sample < samples_per_pixel; ++sample) {
                Vec3f pixel_center = top_left + dx * (j + jitter[0][sample]) +
                                                dy * (i + jitter[1][sample]);
                Vec3f ray = pixel_center - position;
                ray.normalize();
                
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
#if 1
                    r += sphere.c.r * shade;
                    g += sphere.c.g * shade;
                    b += sphere.c.b * shade;
#else
                    r += sphere.c.r;
                    g += sphere.c.g;
                    b += sphere.c.b;
#endif
                }
            }
            scr->draw_pixel(j, i, {Uint8(rint(r / samples_per_pixel)), Uint8(rint(g / samples_per_pixel)), Uint8(rint(b / samples_per_pixel))});
        }
    }
}

void Game::draw_game() {
    scr->cls();
#ifdef THREADED
    vector<future<void>> futures;
    futures.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        futures.push_back(async(launch::async, &Game::render_slice, this, i));
    }
    for (auto& f: futures) {
        f.get();
    }
#else
    render_slice(0);
#endif
}

Game::Game(Screen* scr) :
    scr(scr),
    plane_width(scr->width),
    plane_height(scr->height),
    plane_distance(scr->width / (tan(fov / 2.0) * 2.0)),
    position(0.0, 0.0, 0.0),
    pitch(0.0),
    yaw(0.0),
#ifdef THREADED
    num_threads(4),
#else
    num_threads(1),
#endif
    running(true) {
    scr->set_recording_style("images", 5);
    //SDL_SetRelativeMouseMode(SDL_TRUE);
    //spheres.push_back({{-1.0, 0.0, 10.0}, 1.0, {255, 0, 0}});
    //spheres.push_back({{0.0, 0.0, 10.0}, 1.0, {0, 255, 0}});
    //spheres.push_back({{1.0, 0.0, 10.0}, 1.0, {0, 0, 255}});
    spheres.push_back({{-4.0/3.0, 0.0, 8.0 / 3.0}, 2.0/3.0, {255, 0, 0}});
    spheres.push_back({{0.0, 0.0, 8.0 / 3.0}, 2.0/3.0, {0, 255, 0}});
    spheres.push_back({{4.0/3.0, 0.0, 8.0 / 3.0}, 2.0/3.0, {0, 0, 255}});
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

