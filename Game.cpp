#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include <limits>
#include <random>
#include "Game.h"

constexpr float PI = M_PI;
constexpr float TWO_PI = M_PI * 2.0;
static float focus_length = 1.5;

static constexpr float rad(float theta) {
    return theta * PI / 180.0;
}

static SDL_Color background_color{255, 255, 255};
//static SDL_Color background_color{0, 0, 0};
//static SDL_Color background_color{0, 0, 255};
static float fov = rad(75.0);
static float move_speed = 3.0;
static int cutoff_depth = 20;
static float ambient_light_quantity = 0.08;
static float dof_blur_amount = 0.0;

using namespace std;
using namespace Linear;
using namespace chrono;


static float rand_float(float range_start, float range_end) {
    static thread_local mt19937 gen;
    uniform_real_distribution<float> dis(range_start, range_end);
    return dis(gen);
}

void Game::clear_cells() {
    for (int i = 0; i < scr->width * scr->height; ++i) {
        cells_swap[i].r = 0.0;
        cells_swap[i].g = 0.0;
        cells_swap[i].b = 0.0;
        cells_swap[i].count = 0;
    }
    Cell* swap = cells.load();
    cells.store(cells_swap.load());
    cells_swap.store(swap);
    rays_at_last_swap = rays_cast;
}

static float change_yaw(float yaw, float theta) {
    yaw = fmod(yaw + theta, TWO_PI);
    if (yaw < TWO_PI) {
        yaw += TWO_PI;
    }
    return yaw;
}

static float change_pitch(float pitch, float theta) {
    return max<float>(min<float>(pitch + theta, PI / 2), -PI / 2);
}

void Game::handle_input() {
    static const Uint8* key_state = SDL_GetKeyboardState(NULL);
    SDL_Event event;
    bool recalculate_stuff = false;
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
                mouse_enabled = !mouse_enabled;
                if (mouse_enabled) {
                    SDL_SetRelativeMouseMode(SDL_FALSE);
                } else {
                    SDL_SetRelativeMouseMode(SDL_TRUE);
                }
                break;
            case SDLK_SPACE:
                scr->write_tga("image.tga");
                break;
            default:
                break;
            }
            break;
        case SDL_MOUSEMOTION: {
            if (!mouse_enabled) {
                yaw = change_yaw(yaw, event.motion.xrel * rad(0.75));
                pitch = change_pitch(pitch, -event.motion.yrel * rad(0.23));
                recalculate_stuff = true;
            }
            break;
        }
        default:
             break;
        }
    }
    float speed_factor;
    if (key_state[SDL_SCANCODE_LSHIFT] || key_state[SDL_SCANCODE_LSHIFT]) {
        speed_factor = 4.0;
    } else {
        speed_factor = 1.0;
    }
    if (key_state[SDL_SCANCODE_W]) {
        position += dz * move_speed * scr->frame_time() * speed_factor;
        recalculate_stuff = true;
    }
    if (key_state[SDL_SCANCODE_A]) {
        position -= dx * move_speed * scr->frame_time() * speed_factor;
        recalculate_stuff = true;
    }
    if (key_state[SDL_SCANCODE_D]) {
        position += dx * move_speed * scr->frame_time() * speed_factor;
        recalculate_stuff = true;
    }
    if (key_state[SDL_SCANCODE_S]) {
        position -= dz * move_speed * scr->frame_time() * speed_factor;
        recalculate_stuff = true;
    }
    if (key_state[SDL_SCANCODE_1]) {
        ambient_light_quantity -= 0.01;
        if (ambient_light_quantity < 0.0) {
            ambient_light_quantity = 0.0;
        }
        recalculate_stuff = true;
    }
    if (key_state[SDL_SCANCODE_2]) {
        ambient_light_quantity += 0.01;
        if (ambient_light_quantity > 1.0) {
            ambient_light_quantity = 1.0;
        }
        recalculate_stuff = true;
    }
    if (key_state[SDL_SCANCODE_3]) {
        fov += rad(1.0);
        if (fov > rad(179.0)) {
            fov -= rad(1.0);
        }
        plane_distance = scr->width / (tan(fov / 2.0) * 2.0);
        recalculate_stuff = true;
    }
    if (key_state[SDL_SCANCODE_4]) {
        fov -= rad(1.0);
        if (fov < rad(10.0)) {
            fov += rad(1.0);
        }
        plane_distance = scr->width / (tan(fov / 2.0) * 2.0);
        recalculate_stuff = true;
    }
    if (key_state[SDL_SCANCODE_5]) {
        focus_length -= 0.1;
        if (focus_length < 0.1) {
            focus_length = 0.1;
        }
        clear_cells();
    }
    if (key_state[SDL_SCANCODE_6]) {
        focus_length += 0.1;
        clear_cells();
    }
    if (key_state[SDL_SCANCODE_7]) {
        dof_blur_amount -= 1.0;
        if (dof_blur_amount < 0.0) {
            dof_blur_amount = 0.0;
        }
        clear_cells();
    }
    if (key_state[SDL_SCANCODE_8]) {
        dof_blur_amount += 1.0;
        clear_cells();
    }
    if (recalculate_stuff) {
        calculate_vectors();
        clear_cells();
    }
}

/* ray should be normalized */
float Game::detect_sphere_hit(Vec3f origin, Vec3f ray, int& sphere_index, float cutoff) {
    float scalar_distance = numeric_limits<float>::max();
    sphere_index = -1;
    int size = spheres.size();
    for (int i = 0; i < size; ++i) {
        Sphere& sphere = spheres[i];
        Vec3f distance = sphere.position - origin;
        float b = ray.dot(distance);
        float discriminant = b * b - distance.dot(distance) + sphere.radius * sphere.radius;
        if (discriminant < 0.0) {
            continue;
        }
        float t1 = b + sqrt(discriminant);
        float t0 = b - sqrt(discriminant);
        if (t1 > cutoff && t1 < scalar_distance) {
            scalar_distance = t1;
            sphere_index = i;
        }
        if (t0 > cutoff && t0 < scalar_distance) {
            scalar_distance = t0;
            sphere_index = i;
        }
    }
    return scalar_distance;
}

void Game::calculate_vectors() {
    dz = Vec3f(pitch, yaw);
    dx = Vec3f(0.0, yaw + PI / 2);
    dy = Vec3f(pitch - PI / 2, yaw);
}

void Game::ray_trace(Vec3f origin, Vec3f ray, float& r, float& g, float& b, float weight, int depth) {
    int sphere_index;
    float hit_distance = detect_sphere_hit(origin, ray, sphere_index);
    if (sphere_index != -1) { // We hit a sphere
        Sphere& sphere = spheres[sphere_index];
        Vec3f hit_spot = origin + ray * hit_distance;
        Vec3f sphere_normal = (hit_spot - sphere.position) / sphere.radius;
        Vec3f ref = ray.reflect(sphere_normal);
        for (Light& light: lights) {
            // Sample the light somewhere within its volume to allow for soft
            // shadows
            Vec3f light_sample_point = light.position + rand_float(0.0, 1.0) * light.radius * Vec3f(rand_float(0.0, 2.0 * PI), rand_float(0.0, 2.0 * PI));
            Vec3f light_vec = light_sample_point - hit_spot;
            float light_distance = light_vec.magnitude();
            light_vec.normalize();
            int other_sphere_index; // I don't actually care about this
            float shadow_hit_distance = detect_sphere_hit(hit_spot, light_vec, other_sphere_index);
            // If there is nothing between this sphere and the light
            if (shadow_hit_distance > light_distance) {
                // This value should be positive, but for sitations when
                // the sphere normal is orthogonal to the light vec,
                // it may not be. Thus the abs()
                float diffuse_factor = abs(sphere_normal.dot(light_vec));
                float distance_square = light_distance * light_distance;
                float diffuse_light_quantity = diffuse_factor * light.intensity / distance_square;
                float specular_light_quantity = 0.0;
                if (sphere.spec_power > 0.0) {
                    specular_light_quantity = 255.0 * light.intensity * pow(max<float>(ref.dot(light_vec), 0.0), sphere.spec_power) / distance_square;
                }
                r += (sphere.c.r * diffuse_light_quantity + specular_light_quantity) * weight;
                g += (sphere.c.g * diffuse_light_quantity + specular_light_quantity) * weight;
                b += (sphere.c.b * diffuse_light_quantity + specular_light_quantity) * weight;
            }
        }
        r += sphere.c.r * ambient_light_quantity * weight;
        g += sphere.c.g * ambient_light_quantity * weight;
        b += sphere.c.b * ambient_light_quantity * weight;
        // Perform reflection ray trace
        if (depth + 1 < cutoff_depth && weight * sphere.reflectivity > 1/256.0) {
            ray_trace(hit_spot, ref, r, g, b, weight * sphere.reflectivity, depth + 1);
        }
    } else {
        // Nothing hit
        r += background_color.r * ambient_light_quantity * weight;
        g += background_color.g * ambient_light_quantity * weight;
        b += background_color.b * ambient_light_quantity * weight;
    }
}


void Game::render_loop() {
    uniform_real_distribution<float> dis_y(0.0, plane_height);
    while (running) {
        Cell * target_buf = cells.load();
        uint64_t sample = rays_cast.fetch_add(1);
        sample %= scr->width * scr->height;
        int i = sample / scr->width;
        int j = sample % scr->width;
        float ratio = focus_length / plane_distance;
        float theta = rand_float(0.0, 2.0 * PI);
        float mag = rand_float(0.0, dof_blur_amount);
        Vec3f sample_start = position + mag * ratio * (dx * cos(theta) + dy * sin(theta));
        constexpr float jitter_amount = 1.0;
        float jitter_x = rand_float(0.5 - jitter_amount / 2.0, 0.5 + jitter_amount / 2.0);
        float jitter_y = rand_float(0.5 - jitter_amount / 2.0, 0.5 + jitter_amount / 2.0);
        Vec3f sample_point = position + dz * focus_length +
            ratio * (dx * (jitter_x + j - plane_width / 2.0) +
                     dy * (jitter_y + i - plane_height / 2.0));
        Vec3f ray = sample_point - sample_start;
        ray.normalize();
        float r = 0.0;
        float g = 0.0;
        float b = 0.0;
        ray_trace(sample_start, ray, r, g, b);
        target_buf[sample].r += r;
        target_buf[sample].g += g;
        target_buf[sample].b += b;
        target_buf[sample].count++;
    }
}

Game::Game(Screen* scr, int num_threads) :
    scr(scr),
    cells(new Cell[scr->width * scr->height]),
    cells_swap(new Cell[scr->width * scr->height]),
    rays_cast(0),
    running(true),
    plane_width(scr->width),
    plane_height(scr->height),
    plane_distance(scr->width / (tan(fov / 2.0) * 2.0)),
    pitch(0.0),
    yaw(0.0),
    position(0.0, 0.0, 0.0),
    rays_at_last_swap(0),
    num_threads(num_threads),
    mouse_enabled(false) {
    scr->set_recording_style("images", 5);
    spheres.push_back({{-4.0/3.0, 0.0, 8.0/3.0},
                       2.0/3.0,
                       0.5,
                       100.0,
                       {150, 0, 0}});
    spheres.push_back({{0.0, 0.5, 8.0 / 3.0},
                       2.0/3.0,
                       0.5,
                       100.0,
                       {0, 150, 0}});
    spheres.push_back({{4.0/3.0, 0.0, 8.0 / 3.0},
                       2.0/3.0,
                       0.5,
                       100.0,
                       {0, 0, 150}});
    spheres.push_back({{0.0, 0.0, 7.0},
                       3.0,
                       0.85,
                       15.0,
                       {175, 175, 175}});
    lights.push_back({{0.0, 100.0, 0.0}, 10000.0, 0.5});
    lights.push_back({{-10.0/3.0, 0.0, 0.0}, 6.0, 0.05});

    calculate_vectors();
    clear_cells();
    SDL_SetRelativeMouseMode(SDL_TRUE);
}

Game::~Game() {
    cout << "TOTAL NUMBER OF RAYS CAST: " << rays_cast - rays_at_last_swap << endl;
    delete [] cells.load();
    delete [] cells_swap.load();
}

void Game::run() {
    vector<thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.push_back(thread(&Game::render_loop, this));
    }
    int frame_count = 0;
    uint64_t last_rays_cast = 0;
    high_resolution_clock::time_point last_time = high_resolution_clock::now();
    high_resolution_clock::time_point current_time;
    while (running) {
        for (int i = 0; i < scr->height; ++i) {
            for (int j = 0; j < scr->width; ++j) {
                if (cells[i * scr->width + j].count > 0) {
                    Uint8 r = Uint8(min<int>(rint(cells[i * scr->width + j].r / cells[i * scr->width + j].count), 255));
                    Uint8 g = Uint8(min<int>(rint(cells[i * scr->width + j].g / cells[i * scr->width + j].count), 255));
                    Uint8 b = Uint8(min<int>(rint(cells[i * scr->width + j].b / cells[i * scr->width + j].count), 255));
                    scr->draw_pixel(j, i, {r, g, b});
                }
            }
        }
        handle_input();
        frame_count++;
        if (frame_count == 10) {
            uint64_t this_rays_cast = rays_cast;
            current_time = high_resolution_clock::now();
            cout << "RAYS PER SECOND: " << (this_rays_cast - last_rays_cast) / duration_cast<duration<float>>(current_time - last_time).count() << endl;
            last_rays_cast = this_rays_cast;
            last_time = current_time;
            frame_count = 0;
        }
        scr->commit();
    }
    for (auto& t: threads) {
        t.join();
    }
}

