//
// # Yocto/Model: Examples of procedural modeling
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#ifndef _YOCTO_MODEL_H_
#define _YOCTO_MODEL_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_geometry.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_scene.h>

#include <array>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE OF PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yocto {
struct branch {
    int            id;
    vec3f          start;
    vec3f          end;
    vec3f          direction;
    vec3f          parentdirection;
    vec3f          parentend;
    vector<branch> children = {};
    vector<vec3f>  attractors = {};
    float          thickness;
    bool        operator==(branch b) { return b.id == id; }
    };
struct weightedPoint {
  int    id;
  double weight;
  vec3f  position;
  vec3f  normal;
  bool   operator<(weightedPoint w) { return w.weight < weight; }
};
struct terrain_params {
  float size    = 0.1f;
  vec3f center  = zero3f;
  float height  = 0.1f;
  float scale   = 10;
  int   octaves = 8;
  vec4f bottom  = srgb_to_rgb(vec4f{154, 205, 50, 255} / 255);
  vec4f middle  = srgb_to_rgb(vec4f{205, 133, 63, 255} / 255);
  vec4f top     = srgb_to_rgb(vec4f{240, 255, 255, 255} / 255);
};

void make_terrain(shape_data& shape, const terrain_params& params);

struct displacement_params {
  float height  = 0.02f;
  float scale   = 50;
  int   octaves = 8;
  vec4f bottom  = srgb_to_rgb(vec4f{64, 224, 208, 255} / 255);
  vec4f top     = srgb_to_rgb(vec4f{244, 164, 96, 255} / 255);
  int   type   = 0;
};

void make_displacement(shape_data& shape, const displacement_params& params);

struct hair_params {
  int   num      = 100000;
  int   steps    = 1;
  float lenght   = 0.02f;
  float scale    = 250;
  float strength = 0.01f;
  float gravity  = 0.0f;
  int   samples  = 0;
  vec4f bottom   = srgb_to_rgb(vec4f{25, 25, 25, 255} / 255);
  vec4f top      = srgb_to_rgb(vec4f{244, 164, 96, 255} / 255);
  bool  elimination  = false;
};

struct tree_params {
  int   num              = 200;
  float radius           = 0.14;
  vec3f ellipsoid        = {0.15, 0.12, 0.1};
  float attraction_range = 0.052;
  float step             = 0.001;
  float kill_range       = step*10;
  int   steps            = 2000;
  float randomness       = 0.002;
  float thickness        = 0.009;
  float attack           = 0.0f;
  vec3f origin           = vec3f{-0.2, 0.245, 0.2};
};

struct points {
  int   id;
  vec3f position;
};

struct nearest {
  branch br;
  vector<vec3f> points;
};

void make_hair(
    shape_data& hair, const shape_data& shape, const hair_params& params);

struct grass_params {
  int num = 10000;
};

void make_grass(scene_data& scene, const instance_data& object,
    const vector<instance_data>& grasses, const grass_params& params);

void SampleElimination(
    shape_data& sh, vector<vec3f>& positions, vector<vec3f>& normals, int num);

static bool newcompare(weightedPoint& w1, weightedPoint& w2) {
  return w1.weight < w2.weight;
}


static float calcolateArea(shape_data& shape) {
  auto triangles  = shape.triangles;
  auto qtriangles = quads_to_triangles(shape.quads);
  triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
  return sample_triangles_cdf(triangles, shape.positions).back();
}
static int find_index(vector<weightedPoint>& map, int pos) {
  for (auto i : range(map.size()))
    if (map[i].id == pos) return i;
  return -1;
}

void calcolateWeight(std::vector<weightedPoint>& mapping, hash_grid& grid,
    double rmax, double rmin, vector<vec3f>& positions, vector<vec3f>& normals);

void generate_tree(scene_data& scene, tree_params& params);
vector<vec3f> calculate_attractors(hash_grid& grid, vec3f pos,
    vector<points>& points, vector<vec3f>& active, float r);

void kill_points(scene_data& scene, vec3f pos, vector<vec3f>& active, float k);

static void draw_circle(scene_data& scene, vec3f pos,vec3f color, float r) {
  shape_data    shape = make_sphere(32, r);
  material_data material;
  instance_data instance;
  material.color = color;
  material.type  = material_type::matte;
  scene.shapes.push_back(shape);
  scene.materials.push_back(material);
  instance.shape    = (int)scene.shapes.size() - 1;
  instance.material = (int)scene.materials.size() - 1;
  instance.frame    = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, pos};
  scene.instances.push_back(instance);
}

 vector<branch> assign_branches(scene_data& scene,
    vector<branch>& tree,vec3f pos, vector<vec3f> points, double d);
static void draw_head(scene_data& scene, vec3f pos, tree_params& params) {
  auto          shape = make_sphere(32, params.attraction_range);
  material_data material;
  instance_data instance;

  // attraction range
  material.color = {0.5, 0.8, 0};
  material.type  = material_type::transparent;
  scene.shapes.push_back(shape);
  scene.materials.push_back(material);
  instance.shape    = (int)scene.shapes.size() - 1;
  instance.material = (int)scene.materials.size() - 1;
  instance.frame    = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, pos};
  scene.instances.push_back(instance);

  // kill range
  shape          = make_sphere(32, params.kill_range);
  material.color = {0, 0, 0};
  material.type  = material_type::matte;
  scene.shapes.push_back(shape);
  scene.materials.push_back(material);
  instance.shape    = (int)scene.shapes.size() - 1;
  instance.material = (int)scene.materials.size() - 1;
  instance.frame    = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, pos};
  scene.instances.push_back(instance);
}
int index_of(vector<branch>& branches, int id);
vector<vec3f> calculate_attractors(vector<vec3f>& points, vec3f pos, float r);
static void   draw_branch(scene_data& scene, branch& b,float dL) {
  material_data material;
  instance_data instance;
  auto shape = make_uvcylinder({32, 32, 32}, {b.thickness, dL});
  material.color = vec3f{0, 0, 1};  // srgb_to_rgb(vec3f{214, 134, 56} / 255);
 
  material.type  = material_type::matte;
  scene.shapes.push_back(shape);
  scene.materials.push_back(material);
  instance.shape    = (int)scene.shapes.size() - 1;
  instance.material = (int)scene.materials.size() - 1;
  instance.frame    = frame_fromz(b.start, b.end - b.start);
  scene.instances.push_back(instance);
}
}  // namespace yocto

#endif
