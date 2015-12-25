#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include <random>
#include <unistd.h>
#include "Game.h"

constexpr float PI = M_PI;
constexpr float TWO_PI = M_PI * 2.0;
static float focus_length = 1.5;

static constexpr float rad(float theta) { return theta * PI / 180.0; }

// static SDL_Color background_color{255, 255, 255};
static SDL_Color background_color{100, 100, 100};
// static SDL_Color background_color{0, 0, 255};
constexpr float move_speed = 3.0;
constexpr int cutoff_depth = 20;

using namespace std;
using namespace Linear;
using namespace chrono;

static bool use_random = true;
static float rand_float(float range_start, float range_end) {
  static thread_local mt19937 gen;
  uniform_real_distribution<float> dis(range_start, range_end);
  if (use_random) {
    return dis(gen);
  } else {
    return range_start;
    // return (range_start + range_end) / 2.0;
  }
}

static float change_yaw(float yaw, float theta) {
  yaw = fmod(yaw + theta, TWO_PI);
  if (yaw < TWO_PI) {
    yaw += TWO_PI;
  }
  return yaw;
}

static float change_pitch(float pitch, float theta) {
  return max<float>(min<float>(pitch + theta, PI / 2.0), -PI / 2.0);
}

void Game::handle_input() {
  static const Uint8 *key_state = SDL_GetKeyboardState(NULL);
  SDL_Event event;
  bool change_model = false;
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
        input_disabled = !input_disabled;
        if (input_disabled) {
          SDL_SetRelativeMouseMode(SDL_FALSE);
        } else {
          SDL_SetRelativeMouseMode(SDL_TRUE);
        }
        break;
      case SDLK_SPACE:
        scr->write_bmp("image.tga");
        break;
      case SDLK_9:
        use_random = !use_random;
        change_model = true;
        break;
      default:
        break;
      }
      break;
    case SDL_MOUSEMOTION: {
      if (!input_disabled) {
        yaw = change_yaw(yaw, event.motion.xrel * rad(0.75));
        pitch = change_pitch(pitch, -event.motion.yrel * rad(0.23));
        change_model = true;
      }
      break;
    }
    default:
      break;
    }
  }
  Vec3f unit_dz(pitch, yaw);
  Vec3f unit_dx(0.0, yaw + PI / 2);
  Vec3f unit_dy(pitch - PI / 2, yaw);
  float speed_factor;
  if (!input_disabled) {
    if (key_state[SDL_SCANCODE_LSHIFT] || key_state[SDL_SCANCODE_LSHIFT]) {
      speed_factor = 4.0;
    } else {
      speed_factor = 1.0;
    }
    if (key_state[SDL_SCANCODE_W]) {
      position += unit_dz * move_speed * scr->frame_time() * speed_factor;
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_A]) {
      position -= unit_dx * move_speed * scr->frame_time() * speed_factor;
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_D]) {
      position += unit_dx * move_speed * scr->frame_time() * speed_factor;
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_S]) {
      position -= unit_dz * move_speed * scr->frame_time() * speed_factor;
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_1]) {
      ambient_light_quantity -= 0.01;
      if (ambient_light_quantity < 0.0) {
        ambient_light_quantity = 0.0;
      }
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_2]) {
      ambient_light_quantity += 0.01;
      if (ambient_light_quantity > 1.0) {
        ambient_light_quantity = 1.0;
      }
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_3]) {
      fov += rad(1.0);
      if (fov > rad(179.0)) {
        fov -= rad(1.0);
      }
      plane_distance = scr->width / (tan(fov / 2.0) * 2.0);
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_4]) {
      fov -= rad(1.0);
      if (fov < rad(10.0)) {
        fov += rad(1.0);
      }
      plane_distance = scr->width / (tan(fov / 2.0) * 2.0);
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_5]) {
      focus_length -= 0.1;
      if (focus_length < 0.1) {
        focus_length = 0.1;
      }
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_6]) {
      focus_length += 0.1;
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_7]) {
      dof_blur_amount -= 1.0;
      if (dof_blur_amount < 0.0) {
        dof_blur_amount = 0.0;
      }
      change_model = true;
    }
    if (key_state[SDL_SCANCODE_8]) {
      dof_blur_amount += 1.0;
      change_model = true;
    }
  }
  if (change_model) {
    ratio = focus_length / plane_distance;
    dz = Vec3f(pitch, yaw) * ratio;
    dx = Vec3f(0.0, yaw + PI / 2) * ratio;
    dy = Vec3f(pitch - PI / 2, yaw) * ratio;
    top_left = position + dz * plane_distance - dx * scr->width / 2.0 -
               dy * scr->height / 2.0;
    for (int sample = 0; sample < scr->width * scr->height; ++sample) {
      auto &acc = accumulators[sample];
      acc.lock();
      if (acc.count != 0) {
        acc.r = acc.r / acc.count;
        acc.g = acc.g / acc.count;
        acc.b = acc.b / acc.count;
        acc.count = 0;
      }
      acc.unlock();
    }
  }
}

void Game::calculate_quad_normals() {
  for (Quad &quad : quads) {
    quad.normal = quad.p1.cross(quad.p2);
    quad.normal += quad.p2.cross(quad.p3);
    quad.normal += quad.p3.cross(quad.p4);
    quad.normal += quad.p4.cross(quad.p1);
    quad.normal.normalize();
    planes.push_back(
        {quad.p1, quad.normal, quad.reflectivity, quad.spec_power, quad.c});
  }
  quads.clear();
}

float Game::detect_sphere_hit(Vec3f origin, Vec3f ray, int &sphere_index,
                              float max_distance, float cutoff) {
  sphere_index = -1;
  int size = spheres.size();
  for (int i = 0; i < size; ++i) {
    Sphere &sphere = spheres[i];
    Vec3f distance = sphere.position - origin;
    float b = ray.dot(distance);
    float discriminant =
        b * b - distance.dot(distance) + sphere.radius * sphere.radius;
    if (discriminant < 0.0) {
      continue;
    }
    float t1 = b + sqrt(discriminant);
    float t0 = b - sqrt(discriminant);
    if (t1 > cutoff && t1 < max_distance) {
      max_distance = t1;
      sphere_index = i;
    }
    if (t0 > cutoff && t0 < max_distance) {
      max_distance = t0;
      sphere_index = i;
    }
  }
  return max_distance;
}

/* Taken from wikipedia */
static bool triangle_hit(Vec3f V1, Vec3f V2, Vec3f V3, Vec3f O, Vec3f D,
                         float &dist) {
  Vec3f e1, e2; // Edge1, Edge2
  Vec3f P, Q, T;
  float det, inv_det, u, v;
  float t;

#define EPSILON 0.0005
  e1 = V2 - V1;
  e2 = V3 - V1;
  P = D.cross(e2);
  det = e1.dot(P);
  // NOT CULLING
  if (det > -EPSILON && det < EPSILON)
    return false;
  inv_det = 1.f / det;

  T = O - V1;

  u = T.dot(P) * inv_det;
  if (u < 0.f || u > 1.f)
    return false;

  Q = T.cross(e1);

  v = D.dot(Q) * inv_det;
  if (v < 0.f || u + v > 1.f)
    return false;

  t = e2.dot(Q) * inv_det;

  if (t > EPSILON) { // ray intersection
    dist = t;
    return true;
  }

  // No hit, no win
  return false;
}

float Game::detect_quad_hit(Vec3f origin, Vec3f ray, int &quad_index,
                            float max_distance, float cutoff) {
  quad_index = -1;
  int size = quads.size();
  for (int i = 0; i < size; ++i) {
    Quad &quad = quads[i];
    float distance = numeric_limits<float>::max();
    if (triangle_hit(quad.p1, quad.p2, quad.p3, origin, ray, distance)) {
      if (distance < max_distance) {
        max_distance = distance;
        quad_index = i;
      }
      continue;
    }
    if (triangle_hit(quad.p1, quad.p3, quad.p4, origin, ray, distance)) {
      if (distance < max_distance) {
        max_distance = distance;
        quad_index = i;
      }
      continue;
    }
  }
  return max_distance;
}

float Game::detect_plane_hit(Vec3f origin, Vec3f ray, int &plane_index,
                             float max_distance, float cutoff) {
  plane_index = -1;
  int size = planes.size();
  for (int i = 0; i < size; ++i) {
    Plane &plane = planes[i];
    float t = plane.normal.dot(plane.p - origin) / plane.normal.dot(ray);
    if (t > cutoff && t < max_distance) {
      max_distance = t;
      plane_index = i;
    }
  }
  return max_distance;
}

/* ray should be normalized */
float Game::detect_hit(Vec3f origin, Vec3f ray, int &index,
                       GeometryType &hit_type) {
  index = -1;
  float distance = numeric_limits<float>::max();
  hit_type = GeometryType::NONE;
#if 1
  int sphere_index;
  distance = detect_sphere_hit(origin, ray, sphere_index, distance);
  if (sphere_index != -1) {
    index = sphere_index;
    hit_type = GeometryType::SPHERE;
  }
#endif
#if 0
  int quad_index;
  distance = detect_quad_hit(origin, ray, quad_index, distance);
  if (quad_index != -1) {
    index = quad_index;
    hit_type = GeometryType::QUAD;
  }
#endif
#if 1
  int plane_index;
  distance = detect_plane_hit(origin, ray, plane_index, distance);
  if (plane_index != -1) {
    index = plane_index;
    hit_type = GeometryType::PLANE;
  }
#endif
  return distance;
}

float Game::detect_hit(Vec3f origin, Vec3f ray) {
  int dont_care_index;
  GeometryType dont_care_hit_type;
  return detect_hit(origin, ray, dont_care_index, dont_care_hit_type);
}

void Game::ray_trace(Vec3f origin, Vec3f ray, float &r, float &g, float &b,
                     float weight, int depth) {
  int index;
  GeometryType hit_type;
  // Search through all geometry to dermine if this ray intersects something
  float hit_distance = detect_hit(origin, ray, index, hit_type);
  if (hit_type == GeometryType::NONE) {
    // Nothing was hit, add background color and terminate recursion
    r += background_color.r * ambient_light_quantity * weight;
    g += background_color.g * ambient_light_quantity * weight;
    b += background_color.b * ambient_light_quantity * weight;
    return;
  }
  float spec_power;
  float reflectivity;
  SDL_Color c;
  Vec3f normal;
  Vec3f hit_spot = origin + ray * hit_distance;
  switch (hit_type) {
  case GeometryType::SPHERE:
    normal = (hit_spot - spheres[index].position) / spheres[index].radius;
    spec_power = spheres[index].spec_power;
    reflectivity = spheres[index].reflectivity;
    c = spheres[index].c;
    break;
  case GeometryType::QUAD:
    normal = quads[index].normal;
    spec_power = quads[index].spec_power;
    reflectivity = quads[index].reflectivity;
    c = quads[index].c;
    break;
  case GeometryType::PLANE:
    normal = planes[index].normal;
    spec_power = planes[index].spec_power;
    reflectivity = planes[index].reflectivity;
    c = planes[index].c;
    break;
  default:
    return; // We should never get here
  }
  // sample_vec speeds up the process by putting a hard limit of only three
  // calls to rand_float() per sample taken. This has the effect of making all
  // lights get sampled at the same point relative to their center and size
  static thread_local Vec3f sample_vec;
  if (depth == 0) {
    sample_vec =
        cbrt(rand_float(0.0, 1.0)) *
        Vec3f(rand_float(-PI / 2.0, PI / 2.0), rand_float(0.0, 2.0 * PI));
  }
  Vec3f reflection = ray.reflection(normal);
  for (Light &light : lights) {
    // Sample the light somewhere within its volume to allow for soft
    // shadows
    Vec3f light_sample_point = light.position + light.radius * sample_vec;
    Vec3f light_vec = light_sample_point - hit_spot;
    float light_distance = light_vec.magnitude();
    light_vec.normalize();
    // Determine if there is geometry between the point that we hit and the
    // light that we are sampling
    float shadow_hit_distance = detect_hit(hit_spot, light_vec);
    if (shadow_hit_distance > light_distance) {
      // If there is nothing between this object and the light, illuminate
      // this pixel with the light
      //
      // This value should be positive, but for sitations when
      // the normal is orthogonal to the light vec,
      // it may not be. Thus the abs()
      float diffuse_factor = abs(normal.dot(light_vec));
      float distance_square = light_distance * light_distance;
      float diffuse_light_quantity =
          diffuse_factor * light.intensity / distance_square;
      float specular_light_quantity = 0.0;
      if (spec_power > 0.0) {
        specular_light_quantity =
            255.0 * light.intensity *
            pow(max<float>(reflection.dot(light_vec), 0.0), spec_power) /
            distance_square;
      }
#if 1
      r += (c.r * diffuse_light_quantity + specular_light_quantity) * weight;
      g += (c.g * diffuse_light_quantity + specular_light_quantity) * weight;
      b += (c.b * diffuse_light_quantity + specular_light_quantity) * weight;
#else
      r += (c.r * diffuse_light_quantity * (1.0 - reflectivity) +
            specular_light_quantity) *
           weight;
      g += (c.g * diffuse_light_quantity * (1.0 - reflectivity) +
            specular_light_quantity) *
           weight;
      b += (c.b * diffuse_light_quantity * (1.0 - reflectivity) +
            specular_light_quantity) *
           weight;
#endif
    }
  }
// Add the ambient term
#if 1
  r += c.r * ambient_light_quantity * weight;
  g += c.g * ambient_light_quantity * weight;
  b += c.b * ambient_light_quantity * weight;
#else
  r += c.r * ambient_light_quantity * weight * (1.0 - reflectivity);
  g += c.g * ambient_light_quantity * weight * (1.0 - reflectivity);
  b += c.b * ambient_light_quantity * weight * (1.0 - reflectivity);
#endif
  // If the material is reflective enough, and we haven't already recursed
  // too much, trace again from the hit point with the reflection vector
  if (depth + 1 < cutoff_depth && weight * reflectivity > 1 / 256.0) {
#if 1
#if 1
    ray_trace(hit_spot, reflection, r, g, b, weight * reflectivity, depth + 1);
#else
    float rr = 0.0;
    float rg = 0.0;
    float rb = 0.0;
    ray_trace(hit_spot, reflection, rr, rg, rb, weight * reflectivity,
              depth + 1);
    float luminosity = (0.2126 * rr + 0.7152 * rg + 0.0722 * rb) / 255.0;
    r += rr;
    g += rg;
    b += rb;
    r += c.r * luminosity * (1.0 - reflectivity) * weight;
    g += c.g * luminosity * (1.0 - reflectivity) * weight;
    b += c.b * luminosity * (1.0 - reflectivity) * weight;
#endif
#endif
  }
}

void Game::render_slice(int slice) {
  while (running) {
    for (int i = slice * scr->height / num_threads;
         i < (slice + 1) * scr->height / num_threads; ++i) {
      for (int j = 0; j < scr->width; ++j) {
        int sample = i * scr->width + j;
        float theta = rand_float(0.0, 2.0 * PI);
        // float dof_mag = rand_float(0.0, dof_blur_amount);
        float dof_mag = sqrt(rand_float(0.0, 1.0)) * dof_blur_amount;
        Vec3f sample_start =
            position + dof_mag * (dx * cos(theta) + dy * sin(theta));
        constexpr float jitter_amount = 1.0;
        float jitter_x =
            rand_float(0.5 - jitter_amount / 2.0, 0.5 + jitter_amount / 2.0);
        float jitter_y =
            rand_float(0.5 - jitter_amount / 2.0, 0.5 + jitter_amount / 2.0);
        Vec3f sample_end = top_left + dx * (jitter_x + j) + dy * (jitter_y + i);
        Vec3f ray = (sample_end - sample_start).normal();
        float r = 0.0;
        float g = 0.0;
        float b = 0.0;
        ray_trace(sample_start, ray, r, g, b);
#if 1
        /*r = min<float>(r, 255.0);
        g = min<float>(g, 255.0);
        b = min<float>(b, 255.0);*/
        float m = max<float>(max<float>(r, g), b);
        if (m > 255.0) {
          r *= 255.0 / m;
          g *= 255.0 / m;
          b *= 255.0 / m;
        }
#endif
        auto &acc = accumulators[sample];
        acc.lock();
        if (acc.count == 0) {
          acc.r = r;
          acc.g = g;
          acc.b = b;
        } else {
          acc.r += r;
          acc.g += g;
          acc.b += b;
        }
        acc.count++;
        acc.unlock();
      }
    }
    num_rays_cast += scr->width * ((slice + 1) * scr->height / num_threads -
                                   slice * scr->height / num_threads);
  }
}

Game::Game(PerfSoftScreen *scr, int num_threads)
    : scr(scr), running(true),
      plane_distance(scr->width / (tan(rad(80.0) / 2.0) * 2.0)), pitch(0.0),
      yaw(0.0), fov(rad(80.0)), ambient_light_quantity(0.1),
      dof_blur_amount(0.0), num_rays_cast(0), position(0.0, 0.0, 0.0),
      accumulators(scr->width * scr->height), input_disabled(false),
      num_threads(num_threads) {
  scr->set_recording_style("images", 5);
  spheres.push_back(
      {{-4.0 / 3.0, 0.0, 8.0 / 3.0}, 2.0 / 3.0, 0.25, 100.0, {150, 0, 0}});
  spheres.push_back(
      {{0.0, -0.5, 8.0 / 3.0}, 2.0 / 3.0, 0.25, 100.0, {0, 150, 0}});
  spheres.push_back(
      {{4.0 / 3.0, 0.0, 8.0 / 3.0}, 2.0 / 3.0, 0.25, 100.0, {0, 0, 150}});
  spheres.push_back({{0.0, 1.0, 7.0}, 2.5, 0.9, 15.0, {80, 80, 80}});
  quads.push_back({{-10.0, -3.0, -2.0},
                   {-10.0, -3.0, 18.0},
                   {10.0, -3.0, 18.0},
                   {10.0, -3.0, -2.0},
                   {0.0, 0.0, 0.0},
                   0.0,
                   -1.0,
                   {100, 100, 100}});
#if 0
  quads.push_back({{-10.0, -3.0, -2.0},
                   {-10.0, -3.0, 18.0},
                   {-10.0, 17.0, 18.0},
                   {-10.0, 17.0, -2.0},
                   {0.0, 0.0, 0.0},
                   0.0,
                   -1.0,
                   {0, 100, 100}});
  quads.push_back({{-10.0, 17.0, -2.0},
                   {-10.0, 17.0, 18.0},
                   {10.0, 17.0, 18.0},
                   {10.0, 17.0, -2.0},
                   {0.0, 0.0, 0.0},
                   0.0,
                   -1.0,
                   {0, 0, 0}});
  //{100, 0, 100}});
  quads.push_back({{10.0, 17.0, -2.0},
                   {10.0, 17.0, 18.0},
                   {10.0, -3.0, 18.0},
                   {10.0, -3.0, -2.0},
                   {0.0, 0.0, 0.0},
                   0.0,
                   -1.0,
                   {100, 100, 0}});
  quads.push_back({{-10.0, -3.0, -2.0},
                   {10.0, -3.0, -2.0},
                   {10.0, 17.0, -2.0},
                   {-10.0, 17.0, -2.0},
                   {0.0, 0.0, 0.0},
                   0.0,
                   -1.0,
                   {100, 0, 0}});
  quads.push_back({{-10.0, -3.0, 18.0},
                   {10.0, -3.0, 18.0},
                   {10.0, 17.0, 18.0},
                   {-10.0, 17.0, 18.0},
                   {0.0, 0.0, 0.0},
                   0.0,
                   -1.0,
                   {0, 0, 100}});
#endif
  lights.push_back({{-20.0, 16.0, 8.0}, 300.0, 1.0});
  lights.push_back({{0.0, 16.0, 8.0}, 150.0, 1.0});

  calculate_quad_normals();
  SDL_SetRelativeMouseMode(SDL_TRUE);
  handle_input();
  pitch = 0.0;
  yaw = 0.0;
  for (int sample = 0; sample < scr->width * scr->height; ++sample) {
    accumulators[sample].r = 0.0;
    accumulators[sample].g = 0.0;
    accumulators[sample].b = 0.0;
    accumulators[sample].count = 0;
    accumulators[sample].flag.clear();
  }
  ratio = focus_length / plane_distance;
  dz = Vec3f(pitch, yaw) * ratio;
  dx = Vec3f(0.0, yaw + PI / 2) * ratio;
  dy = Vec3f(pitch - PI / 2, yaw) * ratio;
  top_left = position + dz * plane_distance - dx * scr->width / 2.0 -
             dy * scr->height / 2.0;
}

Game::~Game() { scr->write_bmp("exit.bmp"); }

void Game::run() {
  high_resolution_clock::time_point last_time = high_resolution_clock::now();
  high_resolution_clock::time_point current_time;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (int i = 0; i < num_threads; ++i) {
    threads.emplace_back(std::thread(&Game::render_slice, this, i));
  }
  uint64_t num_rays_cast_last = 0;
  while (running) {
    for (int i = 0; i < scr->height; ++i) {
      for (int j = 0; j < scr->width; ++j) {
        int sample = i * scr->width + j;
        auto &acc = accumulators[sample];
        Uint8 r;
        Uint8 g;
        Uint8 b;
        acc.lock();
        if (acc.count != 0) {
          r = Uint8(rint(min<float>(acc.r / acc.count, 255.0)));
          g = Uint8(rint(min<float>(acc.g / acc.count, 255.0)));
          b = Uint8(rint(min<float>(acc.b / acc.count, 255.0)));
        } else {
          r = Uint8(rint(min<float>(acc.r, 255.0)));
          g = Uint8(rint(min<float>(acc.g, 255.0)));
          b = Uint8(rint(min<float>(acc.b, 255.0)));
        }
        acc.unlock();
        scr->draw_pixel(j, i, {r, g, b});
      }
    }
    handle_input();
    if (!running) {
      break;
    }
    scr->commit();
    current_time = high_resolution_clock::now();
    uint64_t this_num_rays_cast = num_rays_cast;
    static int i = 0;
    if (++i == 20) {
      cout << "RAYS PER SECOND: "
           << (this_num_rays_cast - num_rays_cast_last) /
                  duration_cast<duration<float>>(current_time - last_time)
                      .count()
           << endl;
      num_rays_cast_last = this_num_rays_cast;
      last_time = current_time;
      i = 0;
    }
  }
  for (auto & thread: threads) {
    thread.join();
  }
}
