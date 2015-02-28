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
static float ambient = 0.08;
static float spec_power = 100.0;

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
                if (scr->full_screen) {
                    running = false;
                } else {
                    mouse_enabled = !mouse_enabled;
                    if (mouse_enabled) {
                        SDL_SetRelativeMouseMode(SDL_FALSE);
                    } else {
                        SDL_SetRelativeMouseMode(SDL_TRUE);
                    }
                }
                break;
            case SDLK_RETURN:
                scr->write_tga("image.tga");
                break;
            default:
                break;
            }
            break;
        case SDL_MOUSEMOTION: {
            if (!mouse_enabled) {
                yaw += event.motion.xrel * rad(0.75);
                yaw = clip_yaw(yaw);
                pitch -= event.motion.yrel * rad(0.23);
                pitch = clip_pitch(pitch);
                calculate_vectors();
            }
            break;
        }
        default:
             break;
        }
    }
    if (key_state[SDL_SCANCODE_UP]) {
        pitch += rotate_speed * scr->frame_time();
        pitch = clip_pitch(pitch);
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_DOWN]) {
        pitch -= rotate_speed * scr->frame_time();
        pitch = clip_pitch(pitch);
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_LEFT]) {
        yaw -= rotate_speed * scr->frame_time();
        yaw = clip_yaw(yaw);
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_RIGHT]) {
        yaw += rotate_speed * scr->frame_time();
        yaw = clip_yaw(yaw);
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_W]) {
        position += Vec3f(pitch, yaw) * move_speed * scr->frame_time();
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_A]) {
        position += Vec3f(0.0, yaw - PI / 2) * move_speed * scr->frame_time();
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_D]) {
        position += Vec3f(0.0, yaw + PI / 2) * move_speed * scr->frame_time();
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_S]) {
        position -= Vec3f(pitch, yaw) * move_speed * scr->frame_time();
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_1]) {
        ambient -= 0.01;
        if (ambient < 0.0) {
            ambient = 0.0;
        }
    }
    if (key_state[SDL_SCANCODE_2]) {
        ambient += 0.01;
        if (ambient > 1.0) {
            ambient = 1.0;
        }
    }
    if (key_state[SDL_SCANCODE_3]) {
        spec_power -= 1.0;
        if (spec_power < 1.0) {
            spec_power = 1.0;
        }
    }
    if (key_state[SDL_SCANCODE_4]) {
        spec_power += 1.0;
    }
}

/* direction should be normalized */
float Game::trace(Vec3f origin, Vec3f direction, int& sphere_index, float cutoff) {
    float d = numeric_limits<float>::max();
    sphere_index = -1;
    int size = spheres.size();
    for (int i = 0; i < size; ++i) {
        Sphere& sphere = spheres[i];
        Vec3f distance = sphere.position - origin;
        float B = direction.dot(distance);
        float D = B * B - distance.dot(distance) + sphere.radius * sphere.radius;
        if (D < 0.0) {
            continue;
        }
        float t1 = B + sqrt(D);
        if ((t1 > cutoff && t1 < d)) {
            float t0 = B - sqrt(D);
            if (t0 > cutoff) {
                d = t0;
            } else {
                d = t1;
            }
            sphere_index = i;
        }
    }
    return d;
}

void Game::calculate_vectors() {
    orientation = Vec3f(pitch, yaw);
    dx = Vec3f(0.0, yaw + PI / 2);
    dy = Vec3f(pitch - PI / 2, yaw);
    top_left = position + orientation * plane_distance -
                     dy * (plane_height / 2.0 + 0.5) -
                     dx * (plane_width / 2.0 + 0.5);
}

void Game::render_slice(int slice) {
    int num_samples = samples_per_pixel * lights.size();
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
                int sphere_index;
                float hit_distance = trace(position, ray, sphere_index);
                if (sphere_index != -1) { // We hit a sphere
                    Sphere& sphere = spheres[sphere_index];
                    Vec3f hit = position + ray * hit_distance;
                    for (Light& light: lights) {
                        Vec3f light_vec = light.position - hit;
                        float light_distance = light_vec.magnitude();
                        light_vec.normalize();
                        int other_sphere_index;
                        float shadow_hit_distance = trace(hit, light_vec, other_sphere_index);
                        if (shadow_hit_distance > light_distance) {
                            Vec3f sphere_normal = (hit - sphere.position) / sphere.radius;
                            // This value should be positive, but for sitations when
                            // the sphere normal is orthogonal to the light vec,
                            // it may not be. Thus the abs()
                            float shade = abs(sphere_normal.dot(light_vec));
                            float square = light_distance * light_distance;

                            Vec3f ref = ray.reflect(sphere_normal);
                            float specularity = 255.0 * pow(max<float>(ref.dot(light_vec), 0.0), spec_power);
                            
                            r += sphere.c.r * (shade * light.intensity / square + ambient) + specularity;
                            g += sphere.c.g * (shade * light.intensity / square + ambient) + specularity;
                            b += sphere.c.b * (shade * light.intensity / square + ambient) + specularity;
                        } else {
                            r += sphere.c.r * (ambient);
                            g += sphere.c.g * (ambient);
                            b += sphere.c.b * (ambient);
                        }
                    }
                }
            }
            r = min<int>(rint(r / samples_per_pixel), 255);
            g = min<int>(rint(g / samples_per_pixel), 255);
            b = min<int>(rint(b / samples_per_pixel), 255);
            scr->draw_pixel(j, i, {Uint8(r), Uint8(g), Uint8(b)});
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
    running(true),
    mouse_enabled(false) {
    calculate_vectors();
    scr->set_recording_style("images", 5);

    spheres.push_back({{-4.0/3.0, 0.0, 8.0 / 3.0}, 2.0/3.0, 0.0, {150, 0, 0}});
    spheres.push_back({{0.0, 0.5, 8.0 / 3.0}, 2.0/3.0, 0.0, {0, 150, 0}});
    spheres.push_back({{4.0/3.0, 0.0, 8.0 / 3.0}, 2.0/3.0, 0.0, {0, 0, 150}});
    lights.push_back({{1.0, 2.5, 8.0 / 3.0}, 5.0});
    lights.push_back({{-10.0/3.0, -0.2, 8.0 / 3.0}, 6.0});

    scr->set_color(255, 0, 0);
    jitter[0][0] = 0.0;
    jitter[1][0] = 0.0;
    for (int i = 1; i < samples_per_pixel; ++i) {
        jitter[0][i] = dis(gen);
        jitter[1][i] = dis(gen);
    }
    SDL_SetRelativeMouseMode(SDL_TRUE);
}

Game::~Game() {
}

void Game::run() {
    int i = 0;
    while (running) {
        handle_input();
        draw_game();
        scr->commit();
        if(i++ == 10) {
            cout << "FPS: " << scr->fps() << endl;
            i = 0;
        }
    }
}

