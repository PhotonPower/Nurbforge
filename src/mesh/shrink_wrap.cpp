#include "nurbforge/mesh/shrink_wrap.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>

#include "nurbforge/core/error.hpp"

namespace nf {

namespace {

// Samples lie at most kSpacing * resolution apart (along the grid lines of a face).
constexpr double kSpacing = 0.25;
// Smallest envelope distance that separates inside from outside on the grid, in units of the
// resolution: a grid edge crossing the input has an end within 0.5 of the crossing, which lies
// within kSpacing of a sample.
constexpr double kWall = 0.75;
// Triangles of the envelope are sampled with this many intervals per side.
constexpr int kSubdivision = 3;
// A part of the result must reach this deep below its surface (in units of the resolution).
constexpr double kCore = 0.1;

constexpr double kInf = std::numeric_limits<double>::infinity();

// ------------------------------------------------------------------ sampling

double max_step(const std::vector<Vec3>& p, std::size_t nu, std::size_t nv) {
  double d = 0.0;
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      if (i + 1 < nu) d = std::max(d, norm(p[(i + 1) * nv + j] - p[i * nv + j]));
      if (j + 1 < nv) d = std::max(d, norm(p[i * nv + j + 1] - p[i * nv + j]));
    }
  return d;
}

// Grid of points on the surface with neighbours at most `step` apart (step = inf: coarse).
std::vector<Vec3> sample_surface(const NurbsSurface& s, double step, std::size_t& nu,
                                 std::size_t& nv) {
  nu = nv = 9;
  std::vector<Vec3> pts;
  for (int pass = 0; pass < 8; ++pass) {
    pts = s.sample(nu, nv);
    double du = 0.0, dv = 0.0;
    for (std::size_t i = 0; i < nu; ++i)
      for (std::size_t j = 0; j < nv; ++j) {
        if (i + 1 < nu) du = std::max(du, norm(pts[(i + 1) * nv + j] - pts[i * nv + j]));
        if (j + 1 < nv) dv = std::max(dv, norm(pts[i * nv + j + 1] - pts[i * nv + j]));
      }
    if (du <= step && dv <= step) break;
    // the chords of a coarse grid underestimate the spacing of a fine one: add a margin
    if (du > step)
      nu = static_cast<std::size_t>(std::ceil(static_cast<double>(nu - 1) * du / step * 1.2)) + 1;
    if (dv > step)
      nv = static_cast<std::size_t>(std::ceil(static_cast<double>(nv - 1) * dv / step * 1.2)) + 1;
    if (static_cast<double>(nu) * static_cast<double>(nv) > 5e7)
      throw GeometryError("shrink_wrap: resolution too fine for the size of the input");
  }
  return pts;
}

void sample_curve(const NurbsCurve& c, double step, std::vector<Vec3>& out) {
  std::size_t n = 9;
  std::vector<Vec3> pts;
  for (int pass = 0; pass < 8; ++pass) {
    pts = c.sample(n);
    const double d = max_step(pts, n, 1);
    if (d <= step) break;
    n = static_cast<std::size_t>(std::ceil(static_cast<double>(n - 1) * d / step * 1.2)) + 1;
  }
  out.insert(out.end(), pts.begin(), pts.end());
}

void sample_face(const TrimmedSurface& face, double step, std::vector<Vec3>& out) {
  const NurbsSurface& s = face.surface();
  std::size_t nu = 0, nv = 0;
  const std::vector<Vec3> pts = sample_surface(s, step, nu, nv);
  if (face.is_untrimmed()) {
    out.insert(out.end(), pts.begin(), pts.end());
    return;
  }
  const auto du = s.domain_u(), dv = s.domain_v();
  for (std::size_t i = 0; i < nu; ++i) {
    const double u = du[0] + (du[1] - du[0]) * static_cast<double>(i) / static_cast<double>(nu - 1);
    for (std::size_t j = 0; j < nv; ++j) {
      const double v =
          dv[0] + (dv[1] - dv[0]) * static_cast<double>(j) / static_cast<double>(nv - 1);
      if (face.contains(u, v)) out.push_back(pts[i * nv + j]);
    }
  }
  // the grid stops up to one step short of the trim curves
  for (const auto& loop : face.boundary_curves())
    for (const NurbsCurve& c : loop) sample_curve(c, step, out);
}

std::vector<Vec3> sample(const std::vector<TrimmedSurface>& faces, double step) {
  std::vector<Vec3> pts;
  for (const TrimmedSurface& f : faces) sample_face(f, step, pts);
  if (pts.empty()) throw GeometryError("shrink_wrap: nothing to wrap");
  return pts;
}

// ------------------------------------------------------------------ grid

// Regular grid with one layer of padding; nodes are addressed by their flat index.
struct Grid {
  double h{1.0};
  Vec3 origin;
  std::array<std::size_t, 3> shape{};
  std::array<std::ptrdiff_t, 3> strides{};
  std::size_t size{0};
  std::vector<char> padding;

  Grid(const std::vector<Vec3>& points, double h_, double margin, std::size_t max_nodes,
       bool allocate)
      : h(h_) {
    Vec3 lo{kInf, kInf, kInf}, hi{-kInf, -kInf, -kInf};
    for (const Vec3& p : points) {
      lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
      hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    // +1 node of padding on every side, never touched by the algorithms below
    origin = lo - Vec3{margin + h, margin + h, margin + h};
    const std::array<double, 3> extent{hi.x - lo.x, hi.y - lo.y, hi.z - lo.z};
    double total = 1.0;
    for (std::size_t a = 0; a < 3; ++a) {
      const double count = std::ceil((extent[a] + 2.0 * margin) / h) + 3.0;
      total *= count;
      shape[a] = count < 4e9 ? static_cast<std::size_t>(count) : 0;
    }
    if (total > static_cast<double>(max_nodes)) {
      auto text = [](double v) { return std::to_string(static_cast<unsigned long long>(v)); };
      throw GeometryError(
          "shrink_wrap: the grid would have " + (total < 1e18 ? text(total) : "> 1e18") +
          " nodes (" + std::to_string(shape[0]) + " x " + std::to_string(shape[1]) + " x " +
          std::to_string(shape[2]) + "), more than max_nodes = " + std::to_string(max_nodes) +
          "; use a coarser resolution or raise max_nodes");
    }
    size = shape[0] * shape[1] * shape[2];
    strides = {static_cast<std::ptrdiff_t>(shape[1] * shape[2]),
               static_cast<std::ptrdiff_t>(shape[2]), 1};
    if (!allocate) return;
    padding.assign(size, 0);
    for (std::size_t n = 0; n < size; ++n) {
      const auto c = index(n);
      for (std::size_t a = 0; a < 3; ++a)
        if (c[a] == 0 || c[a] + 1 == shape[a]) padding[n] = 1;
    }
  }

  [[nodiscard]] std::array<std::size_t, 3> index(std::size_t n) const {
    const auto sx = static_cast<std::size_t>(strides[0]);
    return {n / sx, (n % sx) / shape[2], n % shape[2]};
  }
  [[nodiscard]] Vec3 position(std::size_t n) const {
    const auto c = index(n);
    return origin + Vec3{h * static_cast<double>(c[0]), h * static_cast<double>(c[1]),
                         h * static_cast<double>(c[2])};
  }
};

struct DistanceField {
  std::vector<double> d2;          // squared distance to the nearest sample, inf outside the band
  std::vector<std::int32_t> near;  // its index, -1 outside the band
};

// Squared distance to the nearest sample and its index for all nodes within `band` of the
// samples, by propagating nearest samples from node to node.
DistanceField distance_field(const Grid& g, const std::vector<Vec3>& points, double band) {
  if (points.size() > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max()))
    throw GeometryError("shrink_wrap: resolution too fine for the size of the input");
  DistanceField f{std::vector<double>(g.size, kInf), std::vector<std::int32_t>(g.size, -1)};
  // seeds: the eight nodes around every sample
  for (std::size_t p = 0; p < points.size(); ++p) {
    const Vec3 r = (points[p] - g.origin) / g.h;
    const std::array<double, 3> c{std::floor(r.x), std::floor(r.y), std::floor(r.z)};
    const std::size_t base =
        static_cast<std::size_t>(c[0]) * static_cast<std::size_t>(g.strides[0]) +
        static_cast<std::size_t>(c[1]) * g.shape[2] + static_cast<std::size_t>(c[2]);
    for (int corner = 0; corner < 8; ++corner) {
      const int di = corner >> 2, dj = (corner >> 1) & 1, dk = corner & 1;
      const std::size_t node =
          base + static_cast<std::size_t>(di * g.strides[0] + dj * g.strides[1] + dk);
      const Vec3 d = (Vec3{c[0] + di, c[1] + dj, c[2] + dk} - r) * g.h;
      const double dist = dot(d, d);
      if (dist < f.d2[node]) {
        f.d2[node] = dist;
        f.near[node] = static_cast<std::int32_t>(p);
      }
    }
  }
  std::vector<std::size_t> changed, next;
  for (std::size_t n = 0; n < g.size; ++n)
    if (f.near[n] >= 0) changed.push_back(n);

  std::vector<std::ptrdiff_t> shifts;
  for (int di = -1; di <= 1; ++di)
    for (int dj = -1; dj <= 1; ++dj)
      for (int dk = -1; dk <= 1; ++dk)
        if (di != 0 || dj != 0 || dk != 0)
          shifts.push_back(di * g.strides[0] + dj * g.strides[1] + dk);
  const double limit = band * band;
  std::vector<char> marked(g.size, 0);
  while (!changed.empty()) {
    next.clear();
    for (const std::ptrdiff_t shift : shifts)
      for (const std::size_t c : changed) {
        const auto target = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(c) + shift);
        if (g.padding[target]) continue;
        const std::int32_t cand = f.near[c];
        const Vec3 d = g.position(target) - points[static_cast<std::size_t>(cand)];
        const double dist = dot(d, d);
        if (dist < f.d2[target] && dist <= limit) {
          f.d2[target] = dist;
          f.near[target] = cand;
          if (!marked[target]) {
            marked[target] = 1;
            next.push_back(target);
          }
        }
      }
    for (const std::size_t n : next) marked[n] = 0;
    changed.swap(next);
  }
  return f;
}

// Free nodes connected to the start nodes through free nodes (6-neighbourhood).
std::vector<char> flood(const Grid& g, std::vector<char> free,
                        const std::vector<std::size_t>& start) {
  for (std::size_t n = 0; n < g.size; ++n)
    if (g.padding[n]) free[n] = 0;
  std::vector<char> reached(g.size, 0);
  std::vector<std::size_t> stack;
  for (const std::size_t n : start)
    if (free[n]) {
      free[n] = 0;
      reached[n] = 1;
      stack.push_back(n);
    }
  while (!stack.empty()) {
    const std::size_t n = stack.back();
    stack.pop_back();
    for (std::size_t a = 0; a < 3; ++a)
      for (const std::ptrdiff_t s : {g.strides[a], -g.strides[a]}) {
        const auto nb = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(n) + s);
        if (free[nb]) {
          free[nb] = 0;
          reached[nb] = 1;
          stack.push_back(nb);
        }
      }
  }
  return reached;
}

// ------------------------------------------------------------------ marching tetrahedra

// Surface between the outside nodes and the others. Every cube is split into the same six
// tetrahedra around its diagonal (0,0,0)-(1,1,1), so the faces of neighbouring cubes match.
// `fraction(out, other)` places the vertex of a grid edge: the part of the way from its
// outside node to the other one. The normals of the triangles point outside.
template <class Fraction>
TriangleMesh marching_tetrahedra(const Grid& g, const std::vector<char>& outside,
                                 Fraction&& fraction) {
  TriangleMesh mesh;
  std::unordered_map<std::uint64_t, std::size_t> vertex_of;
  auto vertex = [&](std::size_t a, std::size_t b) {
    const std::size_t out = outside[a] ? a : b, other = outside[a] ? b : a;
    const std::uint64_t key = static_cast<std::uint64_t>(std::min(a, b)) * g.size + std::max(a, b);
    const auto [it, inserted] = vertex_of.try_emplace(key, mesh.nodes.size());
    if (inserted) {
      const double t = std::clamp(static_cast<double>(fraction(out, other)), 1e-3, 1.0 - 1e-3);
      const Vec3 pa = g.position(out), pb = g.position(other);
      mesh.nodes.push_back(pa + t * (pb - pa));
    }
    return it->second;
  };
  auto triangle = [&](std::size_t a, std::size_t b, std::size_t c, const Vec3& towards) {
    const Vec3 n = cross(mesh.nodes[b] - mesh.nodes[a], mesh.nodes[c] - mesh.nodes[a]);
    if (dot(n, towards) < 0.0) std::swap(b, c);
    mesh.triangles.push_back({a, b, c});
  };

  // corners of the six tetrahedra as node offsets: 0, e_a, e_a + e_b, (1, 1, 1)
  std::array<std::array<std::ptrdiff_t, 4>, 6> tets{};
  {
    const std::size_t perm[6][2] = {{0, 1}, {0, 2}, {1, 0}, {1, 2}, {2, 0}, {2, 1}};
    for (std::size_t t = 0; t < 6; ++t)
      tets[t] = {0, g.strides[perm[t][0]], g.strides[perm[t][0]] + g.strides[perm[t][1]],
                 g.strides[0] + g.strides[1] + g.strides[2]};
  }
  std::array<std::ptrdiff_t, 8> corners{};
  for (std::size_t c = 0; c < 8; ++c)
    corners[c] = static_cast<std::ptrdiff_t>(c >> 2) * g.strides[0] +
                 static_cast<std::ptrdiff_t>((c >> 1) & 1) * g.strides[1] +
                 static_cast<std::ptrdiff_t>(c & 1);

  for (std::size_t i = 0; i + 1 < g.shape[0]; ++i)
    for (std::size_t j = 0; j + 1 < g.shape[1]; ++j)
      for (std::size_t k = 0; k + 1 < g.shape[2]; ++k) {
        const std::size_t base = i * static_cast<std::size_t>(g.strides[0]) + j * g.shape[2] + k;
        int count = 0;
        for (const std::ptrdiff_t c : corners)
          count += outside[base + static_cast<std::size_t>(c)] ? 0 : 1;
        if (count == 0 || count == 8) continue;
        for (const auto& tet : tets) {
          std::array<std::size_t, 4> n{};
          int inside = 0;
          for (std::size_t c = 0; c < 4; ++c) {
            n[c] = base + static_cast<std::size_t>(tet[c]);
            inside += outside[n[c]] ? 0 : 1;
          }
          if (inside == 0 || inside == 4) continue;
          if (inside == 2) {
            // the two inside corners first
            std::stable_partition(n.begin(), n.end(), [&](std::size_t x) { return !outside[x]; });
            const Vec3 towards =
                g.position(n[2]) + g.position(n[3]) - g.position(n[0]) - g.position(n[1]);
            const std::size_t ac = vertex(n[0], n[2]), ad = vertex(n[0], n[3]);
            const std::size_t bd = vertex(n[1], n[3]), bc = vertex(n[1], n[2]);
            triangle(ac, ad, bd, towards);
            triangle(ac, bd, bc, towards);
          } else {
            const bool lone_inside = inside == 1;
            std::size_t lone = 0;
            while ((outside[n[lone]] == 0) != lone_inside) ++lone;
            std::array<std::size_t, 3> v{};
            Vec3 mean;
            std::size_t m = 0;
            for (std::size_t c = 0; c < 4; ++c) {
              if (c == lone) continue;
              mean += g.position(n[c]) / 3.0;
              v[m++] = vertex(n[lone], n[c]);
            }
            const Vec3 towards = mean - g.position(n[lone]);
            triangle(v[0], v[1], v[2], lone_inside ? towards : -towards);
          }
        }
      }
  if (mesh.triangles.empty()) throw GeometryError("shrink_wrap: the envelope is empty");
  return mesh;
}

}  // namespace

TriangleMesh shrink_wrap(const std::vector<TrimmedSurface>& faces, double resolution,
                         const ShrinkWrapOptions& options) {
  if (!(resolution > 0.0)) throw GeometryError("shrink_wrap: resolution must be positive");
  if (!(options.offset >= 0.0)) throw GeometryError("shrink_wrap: offset must not be negative");
  const double h = resolution, offset = options.offset;
  const double iso = std::max(offset, kWall * h);
  double band = iso + 2.0 * h;
  // fails early if the grid would be too large
  (void)Grid(sample(faces, kInf), h, band, options.max_nodes, false);
  const std::vector<Vec3> points = sample(faces, kSpacing * h);
  const Grid grid(points, h, band, options.max_nodes, true);

  std::vector<char> outside;
  TriangleMesh mesh;
  {
    const DistanceField field = distance_field(grid, points, band);
    std::vector<char> free(grid.size);
    for (std::size_t n = 0; n < grid.size; ++n) free[n] = field.d2[n] > iso * iso ? 1 : 0;
    const auto first = static_cast<std::size_t>(grid.strides[0] + grid.strides[1] + 1);
    outside = options.fill_holes ? flood(grid, std::move(free), {first}) : std::move(free);
    for (std::size_t n = 0; n < grid.size; ++n)
      if (grid.padding[n]) outside[n] = 1;
    // Where the distance falls to iso on the way from the outside node to the other one. The
    // distance is unsigned: if the two nodes see the input in opposite directions it lies
    // between them, and the distance of the second node counts negative.
    mesh = marching_tetrahedra(grid, outside, [&](std::size_t out, std::size_t other) {
      const double d_out = std::min(std::sqrt(field.d2[out]), band);
      double d_other = std::min(std::sqrt(field.d2[other]), band);
      if (field.near[out] >= 0 && field.near[other] >= 0) {
        const Vec3 to_out = grid.position(out) - points[static_cast<std::size_t>(field.near[out])];
        const Vec3 to_other =
            grid.position(other) - points[static_cast<std::size_t>(field.near[other])];
        if (dot(to_out, to_other) < 0.0) d_other = -d_other;
      }
      const double span = d_out - d_other;
      return span > 0.0 ? (d_out - iso) / span : 0.5;
    });
  }

  if (iso > offset) {
    // Closing: take the surface at the distance iso - offset inside the envelope. Away from
    // gaps and concave corners that is the surface at the distance offset from the input.
    const double back = iso - offset;
    std::vector<Vec3> skin;
    skin.reserve(mesh.triangles.size() * 10);
    for (const auto& t : mesh.triangles)
      for (int i = 0; i <= kSubdivision; ++i)
        for (int j = 0; i + j <= kSubdivision; ++j)
          skin.push_back((mesh.nodes[t[0]] * i + mesh.nodes[t[1]] * j +
                          mesh.nodes[t[2]] * (kSubdivision - i - j)) /
                         kSubdivision);
    band = back + 2.0 * h;
    std::vector<double> value(grid.size);
    {
      const DistanceField field = distance_field(grid, skin, band);
      for (std::size_t n = 0; n < grid.size; ++n) {
        const double dist = std::min(std::sqrt(field.d2[n]), band);
        value[n] = outside[n] ? dist + back : back - dist;  // signed distance, positive outside
      }
    }
    // Keep what is connected to a node clearly below the surface: a part thinner than the
    // envelope leaves only noise around zero, which would give scattered fragments.
    std::vector<char> below(grid.size);
    std::vector<std::size_t> core;
    for (std::size_t n = 0; n < grid.size; ++n) {
      below[n] = value[n] <= 0.0 ? 1 : 0;
      if (value[n] < -kCore * h && !grid.padding[n]) core.push_back(n);
    }
    const std::vector<char> solid = flood(grid, std::move(below), core);
    bool any = false;
    for (std::size_t n = 0; n < grid.size; ++n) {
      outside[n] = solid[n] ? 0 : 1;
      if (solid[n])
        any = true;
      else
        value[n] = std::max(value[n], 1e-3 * h);
    }
    if (!any)
      throw GeometryError(
          "shrink_wrap: nothing is left of the envelope; the input is thinner than the "
          "resolution (open surfaces need an offset of at least 0.75 x resolution)");
    mesh = marching_tetrahedra(grid, outside, [&](std::size_t out, std::size_t other) {
      return value[out] / (value[out] - value[other]);
    });
  }

  if (options.min_quality > 0.0) improve_quality(mesh, options.min_quality);
  return mesh;
}

TriangleMesh shrink_wrap(const std::vector<Shell>& shells, double resolution,
                         const ShrinkWrapOptions& options) {
  std::vector<TrimmedSurface> faces;
  for (const Shell& s : shells) faces.insert(faces.end(), s.faces().begin(), s.faces().end());
  return shrink_wrap(faces, resolution, options);
}

}  // namespace nf
