#include "nurbforge/mesh/triangle_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <utility>

#include "nurbforge/core/error.hpp"

namespace nf {

namespace {

// All directed edges (a, b) of the triangles, sorted.
std::vector<std::pair<std::size_t, std::size_t>> directed_edges(const TriangleMesh& m,
                                                                bool undirected) {
  std::vector<std::pair<std::size_t, std::size_t>> e;
  e.reserve(3 * m.triangles.size());
  for (const auto& t : m.triangles)
    for (std::size_t k = 0; k < 3; ++k) {
      std::size_t a = t[k], b = t[(k + 1) % 3];
      if (undirected && a > b) std::swap(a, b);
      e.emplace_back(a, b);
    }
  std::sort(e.begin(), e.end());
  return e;
}

}  // namespace

double TriangleMesh::area() const {
  double s = 0.0;
  for (const auto& t : triangles)
    s += 0.5 * norm(cross(nodes[t[1]] - nodes[t[0]], nodes[t[2]] - nodes[t[0]]));
  return s;
}

double TriangleMesh::volume() const {
  double s = 0.0;
  for (const auto& t : triangles) s += dot(nodes[t[0]], cross(nodes[t[1]], nodes[t[2]]));
  return s / 6.0;
}

std::vector<double> TriangleMesh::quality() const {
  std::vector<double> q(triangles.size());
  for (std::size_t i = 0; i < triangles.size(); ++i) {
    const Vec3 &a = nodes[triangles[i][0]], &b = nodes[triangles[i][1]],
               &c = nodes[triangles[i][2]];
    const double sides = dot(b - a, b - a) + dot(c - b, c - b) + dot(a - c, a - c);
    const double area = 0.5 * norm(cross(b - a, c - a));
    q[i] = sides > 0.0 ? 4.0 * std::sqrt(3.0) * area / sides : 0.0;
  }
  return q;
}

double TriangleMesh::min_quality() const {
  double worst = 1.0;
  for (double q : quality()) worst = std::min(worst, q);
  return worst;
}

bool TriangleMesh::is_watertight() const {
  const auto e = directed_edges(*this, true);
  for (std::size_t i = 0; i < e.size();) {
    std::size_t j = i;
    while (j < e.size() && e[j] == e[i]) ++j;
    if (j - i != 2) return false;
    i = j;
  }
  return true;
}

bool TriangleMesh::is_consistently_oriented() const {
  const auto e = directed_edges(*this, false);
  return std::adjacent_find(e.begin(), e.end()) == e.end();
}

std::size_t improve_quality(TriangleMesh& mesh, double min_quality, std::size_t max_iterations) {
  const std::size_t n = mesh.nodes.size();
  std::vector<Vec3> normal(n), centre(n);
  std::vector<double> degree(n, 0.0);
  // closed mesh: every edge occurs once per direction, so these are all neighbours
  for (const auto& t : mesh.triangles)
    for (std::size_t k = 0; k < 3; ++k) degree[t[k]] += 1.0;
  std::size_t iteration = 0;
  for (; iteration < max_iterations; ++iteration) {
    if (iteration % 5 == 0 && mesh.min_quality() >= min_quality) break;
    std::fill(normal.begin(), normal.end(), Vec3{});
    std::fill(centre.begin(), centre.end(), Vec3{});
    for (const auto& t : mesh.triangles) {
      const Vec3 &a = mesh.nodes[t[0]], &b = mesh.nodes[t[1]], &c = mesh.nodes[t[2]];
      const Vec3 face = cross(b - a, c - a);  // area-weighted
      for (std::size_t k = 0; k < 3; ++k) {
        normal[t[k]] += face;
        centre[t[k]] += mesh.nodes[t[(k + 1) % 3]];
      }
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (degree[i] == 0.0) continue;
      const double len = norm(normal[i]);
      const Vec3 nrm = len > 0.0 ? normal[i] / len : Vec3{};
      Vec3 step = centre[i] / degree[i] - mesh.nodes[i];
      step -= nrm * dot(step, nrm);
      mesh.nodes[i] += 0.6 * step;
    }
  }
  return iteration;
}

void write_stl(const TriangleMesh& mesh, const std::string& path) {
  if (mesh.triangles.size() > std::numeric_limits<std::uint32_t>::max())
    throw GeometryError("write_stl: too many triangles for an STL file");
  std::ofstream f(path, std::ios::binary);
  if (!f) throw GeometryError("write_stl: cannot open '" + path + "'");
  char header[80];
  std::memset(header, ' ', sizeof header);
  const char title[] = "Nurbforge triangle mesh";
  std::memcpy(header, title, sizeof title - 1);
  f.write(header, sizeof header);
  const auto count = static_cast<std::uint32_t>(mesh.triangles.size());
  f.write(reinterpret_cast<const char*>(&count), 4);
  for (const auto& t : mesh.triangles) {
    const Vec3 &a = mesh.nodes[t[0]], &b = mesh.nodes[t[1]], &c = mesh.nodes[t[2]];
    Vec3 nrm = cross(b - a, c - a);
    const double len = norm(nrm);
    if (len > 0.0) nrm /= len;
    float record[12];
    int k = 0;
    for (const Vec3& p : {nrm, a, b, c}) {
      record[k++] = static_cast<float>(p.x);
      record[k++] = static_cast<float>(p.y);
      record[k++] = static_cast<float>(p.z);
    }
    f.write(reinterpret_cast<const char*>(record), sizeof record);
    const std::uint16_t attribute = 0;
    f.write(reinterpret_cast<const char*>(&attribute), 2);
  }
  if (!f) throw GeometryError("write_stl: error writing '" + path + "'");
}

}  // namespace nf
