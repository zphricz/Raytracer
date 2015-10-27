#ifndef GAME_H
#define GAME_H

#include <vector>
#include <limits>
#include <future>
#include "Screen.h"
#include "Vec.h"
#include "SceneComponents.h"

constexpr float ERR_MARGIN = 0.0005;

enum class GeometryType { NONE, SPHERE, QUAD, PLANE };

struct ColorAccumulator {
  float r;
  float g;
  float b;
  uint32_t count;
  std::atomic_flag flag;
  void lock() {
    while (flag.test_and_set(std::memory_order_acquire));
  }
  void unlock() {
    flag.clear();
  }
};

class Game {
private:
  PerfSoftScreen *scr;
  Linear::Vec3f top_left;
  Linear::Vec3f dz;
  Linear::Vec3f dx;
  Linear::Vec3f dy;
  float ratio;
  bool running;
  float plane_distance;
  float pitch;
  float yaw;
  float fov;
  float ambient_light_quantity;
  float dof_blur_amount;
  std::atomic<uint64_t> num_rays_cast;
  Linear::Vec3f position;
  std::vector<ColorAccumulator> accumulators;
  std::vector<Sphere> spheres;
  std::vector<Quad> quads;
  std::vector<Plane> planes;
  std::vector<Light> lights;
  bool input_disabled;

  void render_slice(int slice);
  void handle_input();
  void calculate_quad_normals();

  float
  detect_triangle_hit(Linear::Vec3f origin, Linear::Vec3f ray,
                      int &triangle_index,
                      float max_distance = std::numeric_limits<float>::max(),
                      float cutoff = ERR_MARGIN);
  float detect_quad_hit(Linear::Vec3f origin, Linear::Vec3f ray,
                        int &quad_index,
                        float max_distance = std::numeric_limits<float>::max(),
                        float cutoff = ERR_MARGIN);

  float detect_plane_hit(Linear::Vec3f origin, Linear::Vec3f ray,
                         int &plane_index,
                         float max_distance = std::numeric_limits<float>::max(),
                         float cutoff = ERR_MARGIN);

  float
  detect_sphere_hit(Linear::Vec3f origin, Linear::Vec3f ray, int &sphere_index,
                    float max_distance = std::numeric_limits<float>::max(),
                    float cutoff = ERR_MARGIN);
  float detect_hit(Linear::Vec3f origin, Linear::Vec3f ray, int &index,
                   GeometryType &hit_type);
  float detect_hit(Linear::Vec3f origin, Linear::Vec3f ray);
  void ray_trace(Linear::Vec3f origin, Linear::Vec3f ray, float &r, float &g,
                 float &b, float weight = 1.0, int depth = 0);

public:
  Game(PerfSoftScreen *scr);
  ~Game();
  void run();
};

#endif
