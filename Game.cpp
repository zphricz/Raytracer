#include <iostream>
#include <vector>
#include <future>
#include <algorithm>
#include <limits>
#include <random>
#include "Game.h"

#define THREADED

constexpr float PI = M_PI;
constexpr float TWO_PI = M_PI * 2.0;

static constexpr float rad(float theta) {
    return theta * PI / 180.0;
}

//static SDL_Color background_color{255, 255, 255};
static SDL_Color background_color{0, 0, 0};
static float fov = rad(75.0);
constexpr float move_speed = 3.0;
constexpr float rotate_speed = PI / 2;
constexpr int samples_per_pixel = 1;
constexpr int cutoff_depth = 20;
static float ambient = 0.08;

using namespace std;
using namespace Linear;

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
    float speed_factor;
    if (key_state[SDL_SCANCODE_LSHIFT] || key_state[SDL_SCANCODE_LSHIFT]) {
        speed_factor = 4.0;
    } else {
        speed_factor = 1.0;
    }
    if (key_state[SDL_SCANCODE_W]) {
        position += orientation * move_speed * scr->frame_time() * speed_factor;
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_A]) {
        position -= dx * move_speed * scr->frame_time() * speed_factor;
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_D]) {
        position += dx * move_speed * scr->frame_time() * speed_factor;
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_S]) {
        position -= orientation * move_speed * scr->frame_time() * speed_factor;
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
        fov += rad(1.0);
        if (fov > rad(179.0)) {
            fov -= rad(1.0);
        }
        plane_distance = scr->width / (tan(fov / 2.0) * 2.0);
        calculate_vectors();
    }
    if (key_state[SDL_SCANCODE_4]) {
        fov -= rad(1.0);
        if (fov < rad(10.0)) {
            fov += rad(1.0);
        }
        plane_distance = scr->width / (tan(fov / 2.0) * 2.0);
        calculate_vectors();
    }

}

/* direction should be normalized */
float Game::detect_sphere_hit(Vec3f origin, Vec3f direction, int& sphere_index, float cutoff) {
    float scalar_distance = numeric_limits<float>::max();
    sphere_index = -1;
    int size = spheres.size();
    for (int i = 0; i < size; ++i) {
        Sphere& sphere = spheres[i];
        Vec3f distance = sphere.position - origin;
        float b = direction.dot(distance);
        float discriminant = b * b - distance.dot(distance) + sphere.radius * sphere.radius;
        if (discriminant < 0.0) {
            continue;
        }
        float t1 = b + sqrt(discriminant);
        if ((t1 > cutoff && t1 < scalar_distance)) {
            float t0 = b - sqrt(discriminant);
            if (t0 > cutoff) {
                scalar_distance = t0;
            } else {
                scalar_distance = t1;
            }
            sphere_index = i;
        }
    }
    return scalar_distance;
}

void Game::calculate_vectors() {
    orientation = Vec3f(pitch, yaw);
    dx = Vec3f(0.0, yaw + PI / 2);
    dy = Vec3f(pitch - PI / 2, yaw);
    top_left = position + orientation * plane_distance -
                     dy * (plane_height / 2.0 + 0.5) -
                     dx * (plane_width / 2.0 + 0.5);
}

void Game::ray_trace(Vec3f origin, Vec3f ray, float amount, float& r, float& g, float& b, int depth) {
    int sphere_index;
    float hit_distance = detect_sphere_hit(origin, ray, sphere_index);
    if (sphere_index != -1) { // We hit a sphere
        Sphere& sphere = spheres[sphere_index];
        Vec3f hit_spot = origin + ray * hit_distance;
        Vec3f sphere_normal = (hit_spot - sphere.position) / sphere.radius;
        Vec3f ref = ray.reflect(sphere_normal);
        for (Light& light: lights) {
            Vec3f light_vec = light.position - hit_spot;
            float light_distance = light_vec.magnitude();
            light_vec.normalize();
            int other_sphere_index; // I don't actually care about this
            float shadow_hit_distance = detect_sphere_hit(hit_spot, light_vec, other_sphere_index);
            // If there is nothing between this sphere and the light
            if (shadow_hit_distance > light_distance) {
                // This value should be positive, but for sitations when
                // the sphere normal is orthogonal to the light vec,
                // it may not be. Thus the abs()
                float light_amount = abs(sphere_normal.dot(light_vec));
                float square = light_distance * light_distance;
                float specularity = 0.0;
                if (sphere.spec_power > 0.0) {
                    specularity = 255.0 * light.intensity * pow(max<float>(ref.dot(light_vec), 0.0), sphere.spec_power) / square;
                }
                r += (sphere.c.r * (light_amount * light.intensity / square + ambient) + specularity) * amount * (1.0 - sphere.reflectivity);
                g += (sphere.c.g * (light_amount * light.intensity / square + ambient) + specularity) * amount * (1.0 - sphere.reflectivity);
                b += (sphere.c.b * (light_amount * light.intensity / square + ambient) + specularity) * amount * (1.0 - sphere.reflectivity);
            } else {
                r += amount * sphere.c.r * (ambient) * (1.0 - sphere.reflectivity);
                g += amount * sphere.c.g * (ambient) * (1.0 - sphere.reflectivity);
                b += amount * sphere.c.b * (ambient) * (1.0 - sphere.reflectivity);
            }
        }
        // Perform reflection ray trace
        if (depth + 1 < cutoff_depth && amount * sphere.reflectivity > 1/256.0) {
            ray_trace(hit_spot, ref, amount * sphere.reflectivity, r, g, b, depth + 1);
        }
    } else {
        r += amount * background_color.r;
        g += amount * background_color.g;
        b += amount * background_color.b;
    }
}

void Game::render_slice(int slice) {
    for (int i = slice * scr->height / num_threads; i < (slice + 1) * scr->height / num_threads; ++i) {
        for (int j = 0; j < scr->width; ++j) {
            float r = 0.0;
            float g = 0.0;
            float b = 0.0;
            for (int sample = 0; sample < samples_per_pixel; ++sample) {
                Vec3f sample_point = top_left + dx * (j + jitter[0][sample]) +
                                                dy * (i + jitter[1][sample]);
                Vec3f ray = sample_point - position;
                ray.normalize();
                ray_trace(position, ray, 1.0, r, g, b);
            }
            r = min<int>(rint(r / samples_per_pixel), 255);
            g = min<int>(rint(g / samples_per_pixel), 255);
            b = min<int>(rint(b / samples_per_pixel), 255);
            scr->draw_pixel(j, i, {Uint8(r), Uint8(g), Uint8(b)});
        }
    }
}

void Game::draw_game() {
    scr->fill_screen(background_color);
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

    spheres.push_back({{-4.0/3.0, 0.0, 8.0 / 3.0},
                       2.0/3.0,
                       0.3,
                       100.0,
                       {150, 0, 0}});
    spheres.push_back({{0.0, 0.5, 8.0 / 3.0},
                       2.0/3.0,
                       0.0,
                       100.0,
                       {0, 150, 0}});
    spheres.push_back({{4.0/3.0, 0.0, 8.0 / 3.0},
                       2.0/3.0,
                       0.3,
                       100.0,
                       {0, 0, 150}});
    lights.push_back({{1.0, 100.0, 8.0 / 3.0}, 10000.0});
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

