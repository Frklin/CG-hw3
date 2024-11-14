//
// Implementation for Yocto/Model
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_model.h"

#include <yocto/yocto_sampling.h>

#include "ext/perlin-noise/noise1234.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto


#define alpha 8
#define beta 0.65
#define gamma 1.5
// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXAMPLE OF PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yocto {

vec3f sin3(vec3f x) { return vec3f{sin(x.x), sin(x.y), sin(x.z)}; }
vec4f sin4(vec4f x) { return vec4f{sin(x.x), sin(x.y), sin(x.z),sin(x.w)}; }
vec3f floor3f(vec3f x) { return vec3f{floor(x.x), floor(x.y), floor(x.z)}; }
vec4f floor4f(vec4f x) {
  return vec4f{floor(x.x), floor(x.y), floor(x.z), floor(x.w)};
}
vec3f fract3f(vec3f x) {return x - floor3f(x);}
vec4f fract4f(vec4f x) {return x - floor4f(x);}
vec4f hash4(vec3f p) {
  vec4f q = vec4f{dot(p, vec3f{127.1, 311.7, 411.3}),
      dot(p, vec3f{269.5, 183.3, 152.4}), 
      dot(p, vec3f{419.2, 371.9, 441.0}),
      dot(p,vec3f{512.4, 321.9, 69.420})};
  return fract4f(sin4(q) * 43758.5453f);
}
vec3f hash3(vec3f p) {
  vec3f q = vec3f{dot(p, vec3f{127.1, 311.7, 411.3}), dot(p, vec3f{269.5, 183.3,152.4}),
      dot(p, vec3f{419.2, 371.9,441.0})};
  return fract3f(sin3(q) * 43758.5453f);
  }
vec3f hash3(vec2f p) {
  vec3f q = vec3f{dot(p, vec2f{127.1, 311.7}), 
                  dot(p, vec2f{269.5, 183.3}),
                  dot(p, vec2f{419.2, 371.9})};
  vec3f f = vec3f{floor(sin(q.x) * 43758.5453f), floor(sin(q.y) * 43758.5453f),
      floor(sin(q.z) * 43758.5453f)};

  return vec3f{sin(q.x), sin(q.y), sin(q.z)} * 43758.5453 - f;
}
vec2f hash2(vec2f p) {
  p = vec2f{dot(p, vec2f{127.1, 311.7}), dot(p, vec2f{269.5, 183.3})};
  return vec2f{sin(p.x), sin(p.y)} * 43758.5453 -
         vec2f{floor(sin(p.x) * 43758.5453f), floor(sin(p.y) * 43758.5453f)};
  }


float noise(const vec3f& p) { return ::noise3(p.x, p.y, p.z); }
vec2f noise2(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11})};
}
vec3f noise3(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11}),
      noise(p + vec3f{13, 17, 19})};
}
float fbm(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float turbulence(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float ridge(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 0.5f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * (1 - fabs(noise(p * scale))) * (1 - fabs(noise(p * scale)));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float smoothvoronoi(vec3f x) 
  {
  vec3f p = floor3f(x);
  vec3f  f = fract3f(x);

    float res = 0.0;
    for (int j = -1; j <= 1; j++)
      for (int i = -1; i <= 1; i++)
        for (int n = -1; n <= 1; n++) {
        vec3f b = vec3f{float(j), float(i),float(n)};
        vec3f  r = b - f + hash3(p + b);
        float d = dot(r, r);

        res += 1.0 / std::pow(d, 8.0);
      }
    return pow(1.0 / res, 1.0 / 16.0);
  }
float voronoiDistance(vec3f x) {
  vec3f p = floor3f(x);
  vec3f  f = fract3f(x);

  vec3i mb;
  vec3f  mr;

  float res = 8.0;
  for (int j = -1; j <= 1; j++)
    for (int i = -1; i <= 1; i++)
      for (int n = -1; n <= 1; n++)
    {
        vec3i b = vec3i{i, j, n};
        auto  hash = hash4(p + vec3f{float(i), float(j), float(n)});
        vec3f r = vec3f{float(i), float(j), float(n)} + rgba_to_rgb(hash) - f;
      float d = dot(r, r);

      if (d < res) {
        res = d;
        mr  = r;
        mb  = b;
      }
    }

  res = 8.0;
  for (int j = -2; j <= 2; j++)
    for (int i = -2; i <= 2; i++) 
          for (int n = -2; n <= 2; n++) {
        vec3i b = mb + vec3i{i, j, n};
        //auto  hash = hash4(p + vec3f{float(i), float(j), float(n)});
        vec3f r = vec3f{float(i), float(j), float(n)} - f;
                  /* +rgba_to_rgb(hash)*/
        float d = dot(0.5 * (mr + r), normalize(r - mr));

      res = min(res, d);
    }

  return res;
}

float voronoi(vec3f x, float u, float v) {
    auto p = floor3f(x);
    auto f = fract3f(x);

    float k  = 1.0 + 31.0 * pow(1.0 - v, 4.0); 
    float va = 0.0;
    float wt = 0.0;
    for (int j = -2; j <= 2; j++)
      for (int i = -2; i <= 2; i++) 
         for (int n = -2; n <= 2; n++) {
        vec3f g = vec3f{float(j), float(i), float(n)};
        vec4f o = hash4(p + g) * vec4f{u, u, u,1.0f};
        vec3f r = g - f + vec3f{o.x, o.y, o.z};
        float d = dot(r, r);
        float w = pow(1.0f - smoothstep(0.0f, 1.414f, sqrt(d)), k);
        va += w * o.w;
        wt += w;
      }

    return va / wt;
  }
  void add_polyline(shape_data& shape, const vector<vec3f>& positions,
    const vector<vec4f>& colors, float thickness = 0.0001f) {
  auto offset = (int)shape.positions.size();
  shape.positions.insert(
      shape.positions.end(), positions.begin(), positions.end());
  shape.colors.insert(shape.colors.end(), colors.begin(), colors.end());
  shape.radius.insert(shape.radius.end(), positions.size(), thickness);
  for (auto idx = 0; idx < positions.size() - 1; idx++) {
    shape.lines.push_back({offset + idx, offset + idx + 1});
  }
}

void sample_shape(vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, const shape_data& shape, int num) {
  auto triangles  = shape.triangles;
  auto qtriangles = quads_to_triangles(shape.quads);
  triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
  auto cdf = sample_triangles_cdf(triangles, shape.positions);
  auto rng = make_rng(19873991);
  for (auto idx = 0; idx < num; idx++) {
    auto [elem, uv] = sample_triangles(cdf, rand1f(rng), rand2f(rng));
    auto q          = triangles[elem];
    positions.push_back(interpolate_triangle(
        shape.positions[q.x], shape.positions[q.y], shape.positions[q.z], uv));
    normals.push_back(normalize(interpolate_triangle(
        shape.normals[q.x], shape.normals[q.y], shape.normals[q.z], uv)));
    if (!texcoords.empty()) {
      texcoords.push_back(interpolate_triangle(shape.texcoords[q.x],
          shape.texcoords[q.y], shape.texcoords[q.z], uv));
    } else {
      texcoords.push_back(uv);
    }
  }
}

void make_terrain(shape_data& shape, const terrain_params& params) {
  for (int i = 0; i < shape.positions.size(); i++) {
    auto& position = shape.positions[i];
    position += shape.normals[i] *
                ridge(position * params.scale, params.octaves) * params.height *
                (1.0f - length(position - params.center) / params.size);    
    auto color = position.y/ params.height  > 0.6
                     ? params.top
                     : (position.y/params.height  > 0.3 ? params.middle
                                                           : params.bottom);
    shape.colors.push_back(color);    
  }
  shape.normals = compute_normals(shape);
}

void make_displacement(shape_data& shape, const displacement_params& params) {
  for (int i = 0; i < shape.positions.size(); i++) {
    auto& position = shape.positions[i];
    auto  oldpos   = position;
    switch (params.type) {
      case 0: { //turbolence
        position += shape.normals[i] *
                    turbulence(position * params.scale, params.octaves) *
                    params.height;
        break;
      }
      case 1: { //cell
        position +=
            shape.normals[i] *
            voronoi(position * params.scale, 0.0f, 0.0f) *
            params.height;
        break;
      }      
      case 2: { //voronoi
        position +=
            shape.normals[i] *
            voronoi(position * params.scale, 1.0f, 0.0f) *
            params.height;
        break;
      }      
      case 3: { //voronoise
        position +=
            shape.normals[i] *
            voronoi(position * params.scale, 1.0f, 1.0f) *
            params.height;
        break;
      } 
      case 4: { //smoothvoronoi 
        position +=
            shape.normals[i] *
            smoothvoronoi(position * params.scale) *
            params.height;
        break;
      }      
      case 5: {//true cell noise
        position +=
            shape.normals[i] *
            voronoiDistance(position * params.scale) *
            params.height;
        break;
      }

         
    }
    shape.colors.push_back(
        interpolate_line(params.bottom, params.top, distance(position,oldpos)/params.height));
  }
  shape.normals = compute_normals(shape);
}

void make_hair(
    shape_data& hair, const shape_data& shape, const hair_params& params) {
  auto shape2 = shape;
  vector<vec3f> pos, norm;
  auto          num = params.elimination ? params.num * 5 : params.num;
  sample_shape(pos, norm, shape2.texcoords, shape2, num); //add the condition num*5
  if(params.elimination) SampleElimination(shape2, pos,norm,params.num);
  auto n = (params.lenght / params.steps);
    auto gravity = params.gravity;
  for (int i = 0; i < pos.size(); i++) {
      std::vector<vec3f> positions;
      std::vector<vec4f> colors;
      auto     normal = norm[i];

      positions.push_back(pos[i]);
      colors.push_back(params.bottom);

      for (int k = 0; k < params.steps; k++) {
        vec3f next = positions[k] + n * normal +
                     noise3(positions[k] * params.scale) * params.strength;

        next.y -= gravity;
        normal = normalize(next - positions[k]);
        positions.push_back(next);
        colors.push_back(interpolate_line(params.bottom, params.top,
            distance(next, positions[0]) / params.lenght));
      }
      colors[params.steps] = params.top;
      add_polyline(hair, positions, colors);
  }
  auto tangents = lines_tangents(hair.lines, hair.positions);
  for (int i = 0; i < tangents.size(); i++)
    hair.tangents.push_back(rgb_to_rgba(tangents[i]));
  
}

void make_grass(scene_data& scene, const instance_data& object,
    const vector<instance_data>& grasses, const grass_params& params) {
    
    auto rng = make_rng(34);
  auto& plane = scene.shapes[object.shape];
    sample_shape(
        plane.positions, plane.normals, plane.texcoords, plane, params.num);
    for (int i = 0; i < plane.positions.size(); i++) {
      instance_data grassinstance;
      auto          grass = grasses[rand1i(rng, grasses.size())];

      grassinstance.shape    = grass.shape;
      grassinstance.material = grass.material;

      // newframe
     // grassinstance.frame = grass.frame;
     grassinstance.frame.y   = plane.normals[i];
     grassinstance.frame.x = normalize(vec3f{1, 0, 0} - dot(vec3f{1, 0, 0}, grassinstance.frame.y) * grassinstance.frame.y);
     grassinstance.frame.z = cross(grassinstance.frame.x, grassinstance.frame.y);
     grassinstance.frame.o = plane.positions[i];

      
  // TRASORMAZIONE
      // Scaling
      float rand = 0.9f + rand1f(rng) * 0.1f;
      grassinstance.frame *= scaling_frame(vec3f{rand, rand, rand});
      // Rotation y
      rand = rand1f(rng) * 2 * pi;
      grassinstance.frame *= rotation_frame(grassinstance.frame.y, rand);
      // Rotation z
      rand = 0.1f +rand1f(rng) * 0.1f;
      grassinstance.frame *= rotation_frame(grassinstance.frame.z, rand);

      
     scene.instances.push_back(grassinstance);
    }

}




void SampleElimination(
    shape_data& shape, vector<vec3f>& positions, vector<vec3f>& normals,int num) {
  vector<weightedPoint> mapping;
  // RMAX E RMIN
  auto rmax = std::pow(calcolateArea(shape) / (4.0f * std::sqrt(2) * num), 1.0 / 3);
  auto rmin = rmax * (1 - pow(num / num * 5, gamma)) * beta;

  // GRID AND SIZE
  auto grid = make_hash_grid(positions, rmax);
  auto size = positions.size();

  // INITIALIZING THE VECTOR ID WEIGHT POSITION NORMAL
  calcolateWeight(mapping, grid, rmax, rmin,positions,normals);
  printf("weights\n");

  int passo = 0;
  // START ELIMINATION
  while (size-- > num) {
    printf("%d/%d\n", size, num);
    // ORDINATE THE MAP
    make_heap(mapping.begin(), mapping.end()-passo, newcompare);
    auto eliminated      = mapping.front();  

    // SWAP THE MAX WITH THE LAST ELEMENT
    pop_heap(mapping.begin(), mapping.end()-passo);
    passo++;

    // CREATE THE LIST OF THE NEIGHBORS OF THE MAX
    std::vector<int> neighbors;
    find_neighbors(grid, neighbors, eliminated.position, 2 * rmax); 

    // FOR EACH NEIGHBOR SUBTRACT THE WEIGHT OF THE DELETED ELEMENT
    for (auto& neighbor : neighbors) {
      int i = find_index(mapping, neighbor);
      mapping[i].weight -= pow(1 - distance(positions[neighbor], eliminated.position) /(2 * rmax),alpha);
    }
  }
  printf("while\n");
  for (auto i : range(passo)) mapping.pop_back();
  positions.clear();
  normals.clear();

  // NOW THE SOLUTION ARE ELEMENT FROM 0 TO SIZE - PASSO OF MAPPING
  for (auto& el : mapping) {
    positions.push_back(el.position);
    normals.push_back(el.normal);
  }

}


void calcolateWeight(std::vector<weightedPoint>& mapping, hash_grid& grid,
    double rmax, double rmin, vector<vec3f>& positions,
    vector<vec3f>& normals) {
  weightedPoint point;
  for (auto i : range(positions.size())) {
    point.id = i;
    point.weight   = 0.0;
    point.position = positions[i];
    point.normal   = normals[i];
    vector<int> neighbors;
    find_neighbors(grid, neighbors, positions[i], 2*rmax);
    for (auto neighbor : neighbors) {
      double d    = distance(positions[i],positions[neighbor]);  // distanza sulla shape non generale
      auto   d_ij = d > 2 * rmin ? min(d, 2 * rmax) : 2 * rmin;
      point.weight += pow(1 - d_ij / (2 * rmax), alpha);
    }
    mapping.push_back(point);
  }
}



void generate_tree(scene_data& scene, tree_params& params) {
  vector<vec3f>  active_pos;
  vector<branch> tree,near;
  int            passo = 0, id = 0;
  auto          cell_size = 0.01;  
  auto          rng       = make_rng(100899);

  // generate points for leaves
   for (auto i : range(params.num)) {
    float sigma = rand1f(rng) * pi;
    float theta = rand1f(rng) * 2 * pi;
    auto  D     = vec3f{cos(theta) * sin(sigma), sin(theta) * sin(sigma), cos(sigma)};
    auto d = pow(sin(rand1f(rng) * pi / 2), 0.8) *params.radius;
    auto  pos   = D * d + params.origin;
    active_pos.push_back(pos);
   draw_circle(scene, active_pos.back(), vec3f{1, 0, 0},0.003);
  } 
   //elipse
  /* for (auto i : range(params.num)) {
    float phi = asin(rand1f(rng)*2-1.0);
    float theta = rand1f(rng) * 2 * pi;
    float x = params.ellipsoid.x* cos(theta) * sin(phi);
    float y  = params.ellipsoid.y *sin(phi)*sin(theta);
    float z     = params.ellipsoid.z * cos(phi);
    active_pos.push_back(vec3f{x, y, z}+ params.origin);
    //draw_circle(scene, active_pos.back(), vec3f{1, 0, 0}, 0.003);
  }*/

   printf("points generated\n");


  // generate the first branch int index_of(vector<branch>& branches, int id)
  branch son,parent;
  vector<vec3f> attractors; 
  parent.end = vec3f{params.origin.x, 0,params.origin.z};
  parent.direction = {0, 1, 0};
  parent.thickness = 0.003;
 // params.thickness;
  while (active_pos.empty() ||params.steps-- > 0) 
  {
       attractors = calculate_attractors(active_pos,parent.end,params.attraction_range);   
       //calcolare i vicini
       if (attractors.empty() || parent.end.y < params.attack) {
         son.id        = id++;
         son.start     = parent.end;
         son.direction = normalize(parent.direction + (rand3f(rng) * 2-one3f)*params.randomness);
         //       son.direction = parent.direction + normalize(rand3f(rng))*params.randomness;
         son.thickness =
             parent.thickness;  // max(parent.thickness  - 0.00001,0.00001);
                                // //usare la distanza dal centro?
         son.end = son.start + params.step * son.direction;
         if (distance(son.start, params.origin) < params.radius &&
             distance(son.end, params.origin) > params.radius)
           continue;
         draw_branch(scene, son, params.step);

         kill_points(scene, son.end, active_pos, params.kill_range);
         tree.push_back(son);
       } 
       else {
         near = assign_branches(scene, tree, parent.end, attractors, params.attraction_range);
         for (auto& father : near) {
           son.id = id++;

           son.start = father.end;

           auto direction = zero3f;

           for (auto& neighbor : father.attractors) {
             direction += normalize(neighbor - son.start);
             draw_circle(scene, neighbor, {0, 1, 0}, 0.0032);
           }

           son.direction = normalize(direction + (rand3f(rng) * 2-one3f) * params.randomness);
           //son.direction = normalize(direction) + normalize(rand3f(rng))*params.randomness;

           float t       = father == parent ? 1.2 : 3.5; 
           son.thickness = parent.thickness;//max(father.thickness -t * 0.00001,0.00001);

           son.end = son.start + params.step * son.direction;
            if
            (distance(son.start, params.origin) < params.radius &&
               distance(son.end, params.origin) > params.radius)
             continue;

           tree.push_back(son);

           draw_branch(scene, son, params.step);

           kill_points(scene, son.end, active_pos, params.kill_range);
         }
       }
      parent = tree[passo++];


  }
  printf("fine\n");
  draw_head(scene, son.end, params);

}

vector<branch> assign_branches(scene_data& scene, vector<branch>& tree,
    vec3f pos, vector<vec3f> points, double d) {
  vector<branch> branches, accepted;
  float           dmin,d2;
  branch&           right = tree[0];

  //take the branches around
  for (int i = tree.size() - 1; i >= 0; i--) //provare con tutto l'albero
    if (distance(pos, tree[i].end) < d) branches.push_back(tree[i]);

  for (auto& point : points) {
    dmin = INFINITY;
    for (auto& br : branches) { 
      d2    = dmin;
      dmin = min(dmin, distance(point, br.end)); 
      if (d2 != dmin) right = br;
    }
    auto index = index_of(accepted, right.id);
    if (index != -1)
        accepted[index].attractors.push_back(point);
    else{ 
      right.attractors.push_back(point);
      accepted.push_back(right);
    }

  }
  return accepted;
 }

int index_of(vector<branch>& branches, int id) {
   for (auto const& b : branches)  
     if (b.id == id)
       return std::addressof(b) - std::addressof(branches[0]);
   return -1;  
 }

void kill_points(scene_data& scene, vec3f pos,vector<vec3f>& active, float k) {
    for (auto i : range(active.size())) 
        if (distance(pos, active[i]) < k) {
   //   draw_circle(scene, active[i], zero3f,0.004);
      active.erase(active.begin() + i);
    }
  }

  vector<vec3f> calculate_attractors(vector<vec3f>& points, vec3f pos, float r) {
  vector<vec3f> res;
  for (auto& point : points) 
    if (distance(point, pos) < r) 
      res.push_back(point);
 return res;
}




}  // namespace yocto
