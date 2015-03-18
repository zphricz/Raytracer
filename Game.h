#ifndef GAME_H
#define GAME_H

#include <vector>
#include <atomic>
#include <limits>
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
  uint32_t model_number;
};

class Game {
private:
  SoftScreen *scr;
  std::atomic<uint32_t> model_number;
  std::atomic<uint64_t> rays_cast;
  std::atomic<bool> running;
  float plane_distance;
  float pitch;
  float yaw;
  float fov;
  float ambient_light_quantity;
  float dof_blur_amount;
  Linear::Vec3f position;
  std::vector<ColorAccumulator> accumulators;
  std::vector<Sphere> spheres;
  std::vector<Quad> quads;
  std::vector<Plane> planes;
  std::vector<Light> lights;
  uint64_t rays_at_last_model_change;
  const int num_threads;
  bool input_disabled;

  void render_loop();
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
  Game(SoftScreen *scr, int num_threads);
  ~Game();
  void run();
};

#endif
