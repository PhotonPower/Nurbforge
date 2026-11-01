#include "nurbforge/io/step.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

#include "nurbforge/core/error.hpp"
#include "nurbforge/topology/boolean.hpp"
#include "nurbforge/version.hpp"

namespace nf {

namespace {

// ------------------------------------------------------------------ tokens

/// STEP REAL: shortest round-trip representation, always with a decimal point ("1.", "1.E-05").
std::string real(double v) {
  if (!std::isfinite(v)) throw GeometryError("write_step: non-finite coordinate");
  char buf[64];
  const auto res = std::to_chars(buf, buf + sizeof(buf), v);
  const std::string s(buf, res.ptr);
  const auto e = s.find('e');
  std::string out = s.substr(0, e);
  if (out.find('.') == std::string::npos) out += '.';
  if (e != std::string::npos) out += "E" + s.substr(e + 1);
  return out;
}

/// STEP STRING: quotes doubled, backslashes escaped, non-printable / non-ASCII characters replaced.
std::string str(const std::string& s) {
  std::string out = "'";
  for (const char c : s) {
    const auto u = static_cast<unsigned char>(c);
    if (c == '\'') out += "''";
    else if (c == '\\') out += "\\\\";
    else if (u >= 32 && u < 127) out += c;
    else out += '_';
  }
  return out + "'";
}

std::string ref(std::size_t id) { return "#" + std::to_string(id); }

std::string logical(bool b) { return b ? ".T." : ".F."; }

template <class F>
std::string list(std::size_t n, F&& item) {
  std::string out = "(";
  for (std::size_t i = 0; i < n; ++i) {
    if (i) out += ',';
    out += item(i);
  }
  return out + ")";
}

std::string utc_now() {
  const long long secs = static_cast<long long>(
      std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch())
          .count());
  const long long rem = secs % 86400;
  // civil_from_days (H. Hinnant)
  const long long day_count = secs / 86400 + 719468;
  const long long era = (day_count >= 0 ? day_count : day_count - 146096) / 146097;
  const long long doe = day_count - era * 146097;
  const long long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  long long y = yoe + era * 400;
  const long long doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const long long mp = (5 * doy + 2) / 153;
  const long long d = doy - (153 * mp + 2) / 5 + 1;
  const long long m = mp < 10 ? mp + 3 : mp - 9;
  if (m <= 2) ++y;
  char buf[160];
  std::snprintf(buf, sizeof(buf), "%04lld-%02lld-%02lldT%02lld:%02lld:%02lld", y, m, d, rem / 3600,
                (rem / 60) % 60, rem % 60);
  return buf;
}

// ------------------------------------------------------------------ entity store

class Writer {
 public:
  std::size_t add(std::string body) {
    data_.push_back(std::move(body));
    return data_.size();
  }
  [[nodiscard]] const std::vector<std::string>& data() const noexcept { return data_; }

 private:
  std::vector<std::string> data_;
};

std::size_t point(Writer& w, const Vec3& p) {
  return w.add("CARTESIAN_POINT('',(" + real(p.x) + "," + real(p.y) + "," + real(p.z) + "))");
}

std::size_t direction(Writer& w, const Vec3& d) {
  return w.add("DIRECTION('',(" + real(d.x) + "," + real(d.y) + "," + real(d.z) + "))");
}

struct Knots {
  std::string mults, values;
};

Knots knots(const KnotVector& kv) {
  const auto unique = unique_knots(kv);
  const std::vector<double>& values = unique.first;
  const std::vector<int>& mults = unique.second;
  return {list(mults.size(), [&](std::size_t i) { return std::to_string(mults[i]); }),
          list(values.size(), [&](std::size_t i) { return real(values[i]); })};
}

bool rational(const std::vector<double>& w) {
  return std::any_of(w.begin(), w.end(), [](double x) { return x != 1.0; });
}

double extent(const std::vector<Vec3>& pts) {
  double e = 0.0;
  for (const auto& p : pts) e = std::max({e, std::abs(p.x), std::abs(p.y), std::abs(p.z)});
  return std::max(e, 1.0);
}

std::size_t write_curve(Writer& w, const NurbsCurve& c) {
  const auto& P = c.control_points();
  std::vector<std::size_t> ids;
  for (const auto& p : P) ids.push_back(point(w, p));
  const std::string pts = list(ids.size(), [&](std::size_t i) { return ref(ids[i]); });
  const Knots k = knots(c.knots());
  const std::string deg = std::to_string(c.degree());
  const bool closed = distance(c.evaluate(c.domain_start()), c.evaluate(c.domain_end())) <= 1e-12 * extent(P);
  const std::string form = ",.UNSPECIFIED.," + logical(closed) + ",.F.";
  if (!rational(c.weights()))
    return w.add("B_SPLINE_CURVE_WITH_KNOTS(''," + deg + "," + pts + form + "," + k.mults + "," + k.values +
                 ",.UNSPECIFIED.)");
  const auto& W = c.weights();
  return w.add("(BOUNDED_CURVE() B_SPLINE_CURVE(" + deg + "," + pts + form + ") B_SPLINE_CURVE_WITH_KNOTS(" +
               k.mults + "," + k.values + ",.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() " +
               "RATIONAL_B_SPLINE_CURVE(" + list(W.size(), [&](std::size_t i) { return real(W[i]); }) +
               ") REPRESENTATION_ITEM(''))");
}

std::size_t write_surface(Writer& w, const NurbsSurface& s) {
  const std::size_t nu = s.num_u(), nv = s.num_v();
  std::vector<std::size_t> ids;
  for (const auto& p : s.control_points()) ids.push_back(point(w, p));
  // control_points_list: outer index along u, inner along v (the storage order)
  const std::string net = list(nu, [&](std::size_t i) {
    return list(nv, [&](std::size_t j) { return ref(ids[i * nv + j]); });
  });
  const double tol = 1e-12 * extent(s.control_points());
  bool u_closed = true, v_closed = true;
  for (std::size_t j = 0; j < nv; ++j)
    u_closed = u_closed && distance(s.control_point(0, j), s.control_point(nu - 1, j)) <= tol;
  for (std::size_t i = 0; i < nu; ++i)
    v_closed = v_closed && distance(s.control_point(i, 0), s.control_point(i, nv - 1)) <= tol;
  const Knots ku = knots(s.knots_u()), kv = knots(s.knots_v());
  const std::string head = std::to_string(s.degree_u()) + "," + std::to_string(s.degree_v()) + "," + net +
                           ",.UNSPECIFIED.," + logical(u_closed) + "," + logical(v_closed) + ",.F.";
  const std::string kn = ku.mults + "," + kv.mults + "," + ku.values + "," + kv.values + ",.UNSPECIFIED.";
  if (!rational(s.weights())) return w.add("B_SPLINE_SURFACE_WITH_KNOTS(''," + head + "," + kn + ")");
  const std::string weights = list(nu, [&](std::size_t i) {
    return list(nv, [&](std::size_t j) { return real(s.weight(i, j)); });
  });
  return w.add("(BOUNDED_SURFACE() B_SPLINE_SURFACE(" + head + ") B_SPLINE_SURFACE_WITH_KNOTS(" + kn +
               ") GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(" + weights +
               ") REPRESENTATION_ITEM('') SURFACE())");
}

// ------------------------------------------------------------------ topology

struct OrientedUse {
  std::size_t edge;
  bool reversed;
  std::size_t index;
};

/// The edge uses of loop l of face f in loop order. Uses are sorted by trim-curve index; the
/// pieces of one trim curve (split at T-junctions while sewing) are chained by their vertices.
std::vector<OrientedUse> loop_uses(const Shell& sh, std::size_t f, std::size_t l) {
  const auto& E = sh.edges();
  std::vector<OrientedUse> uses;
  for (std::size_t e = 0; e < E.size(); ++e)
    for (const auto& u : E[e].uses)
      if (u.face == f && u.loop == l) uses.push_back({e, u.reversed, u.index});
  std::stable_sort(uses.begin(), uses.end(),
                   [](const OrientedUse& a, const OrientedUse& b) { return a.index < b.index; });
  const auto start = [&](const OrientedUse& o) { return o.reversed ? E[o.edge].v1 : E[o.edge].v0; };
  const auto end = [&](const OrientedUse& o) { return o.reversed ? E[o.edge].v0 : E[o.edge].v1; };
  std::vector<OrientedUse> out;
  for (std::size_t i = 0; i < uses.size();) {
    std::size_t j = i;
    while (j < uses.size() && uses[j].index == uses[i].index) ++j;
    std::vector<OrientedUse> group(uses.begin() + static_cast<std::ptrdiff_t>(i),
                                   uses.begin() + static_cast<std::ptrdiff_t>(j));
    while (!group.empty()) {
      std::size_t pick = 0;
      bool found = false;
      if (!out.empty())
        for (std::size_t k = 0; k < group.size() && !found; ++k)
          if (start(group[k]) == end(out.back())) {
            pick = k;
            found = true;
          }
      if (!found)  // chain start: a piece that does not continue another piece of the group
        for (std::size_t k = 0; k < group.size() && !found; ++k) {
          bool continues = false;
          for (std::size_t m = 0; m < group.size(); ++m)
            if (m != k && end(group[m]) == start(group[k])) continues = true;
          if (!continues) {
            pick = k;
            found = true;
          }
        }
      out.push_back(group[pick]);
      group.erase(group.begin() + static_cast<std::ptrdiff_t>(pick));
    }
    i = j;
  }
  return out;
}

/// Writes the vertices, edges and faces of one shell (each vertex / edge once).
class ShellWriter {
 public:
  ShellWriter(Writer& w, const Shell& sh)
      : w_(w), sh_(sh), vertex_(sh.vertices().size(), 0), edge_(sh.edges().size(), 0) {}

  /// ADVANCED_FACE for face f. `flip`: the face normal is opposite to the surface normal
  /// (same_sense .F.; the loops, oriented counter-clockwise about the surface normal, are then
  /// clockwise about the face normal: bound orientation .F.).
  std::size_t face(std::size_t f, bool flip) {
    const TrimmedSurface& t = sh_.faces()[f];
    const std::size_t surf = write_surface(w_, t.surface());
    std::vector<std::size_t> bounds;
    for (std::size_t l = 0; l <= t.holes().size(); ++l) {
      const auto uses = loop_uses(sh_, f, l);
      if (uses.empty()) continue;  // only degenerate pieces (pole)
      std::vector<std::size_t> oe;
      for (const auto& u : uses)
        oe.push_back(w_.add("ORIENTED_EDGE('',*,*," + ref(edge(u.edge)) + "," + logical(!u.reversed) + ")"));
      const std::size_t loop = w_.add("EDGE_LOOP(''," + list(oe.size(), [&](std::size_t i) { return ref(oe[i]); }) + ")");
      bounds.push_back(w_.add(std::string(l == 0 ? "FACE_OUTER_BOUND" : "FACE_BOUND") + "(''," + ref(loop) + "," +
                              logical(!flip) + ")"));
    }
    return w_.add("ADVANCED_FACE(''," + list(bounds.size(), [&](std::size_t i) { return ref(bounds[i]); }) + "," +
                  ref(surf) + "," + logical(!flip) + ")");
  }

 private:
  Writer& w_;
  const Shell& sh_;
  std::vector<std::size_t> vertex_, edge_;

  std::size_t vertex(std::size_t i) {
    if (!vertex_[i]) vertex_[i] = w_.add("VERTEX_POINT(''," + ref(point(w_, sh_.vertices()[i])) + ")");
    return vertex_[i];
  }
  std::size_t edge(std::size_t e) {
    if (!edge_[e]) {
      const ShellEdge& se = sh_.edges()[e];
      const std::size_t a = vertex(se.v0), b = vertex(se.v1);
      const std::size_t c = write_curve(w_, se.curve);
      edge_[e] = w_.add("EDGE_CURVE(''," + ref(a) + "," + ref(b) + "," + ref(c) + ",.T.)");
    }
    return edge_[e];
  }
};

/// Connected components of the faces (through edges used by two faces).
std::vector<std::vector<std::size_t>> components(const Shell& sh) {
  std::vector<std::size_t> parent(sh.faces().size());
  std::iota(parent.begin(), parent.end(), std::size_t{0});
  const auto find = [&](std::size_t x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
  };
  for (const auto& e : sh.edges())
    for (std::size_t k = 1; k < e.uses.size(); ++k) parent[find(e.uses[k].face)] = find(e.uses[0].face);
  std::vector<std::vector<std::size_t>> out;
  std::vector<std::size_t> slot(parent.size(), std::numeric_limits<std::size_t>::max());
  for (std::size_t f = 0; f < parent.size(); ++f) {
    const std::size_t r = find(f);
    if (slot[r] == std::numeric_limits<std::size_t>::max()) {
      slot[r] = out.size();
      out.emplace_back();
    }
    out[slot[r]].push_back(f);
  }
  return out;
}

/// Signed volume enclosed by a closed component (divergence theorem, as Shell::volume).
double signed_volume(const Shell& sh, const std::vector<std::size_t>& faces) {
  double size = 1e-300;
  for (std::size_t f : faces) size = std::max(size, extent(sh.faces()[f].surface().control_points()));
  const double abs_tol = 1e-9 * size * size * size / static_cast<double>(faces.size());
  double V = 0.0;
  for (std::size_t f : faces) {
    const NurbsSurface& s = sh.faces()[f].surface();
    V += sh.faces()[f].integrate(
        [&](double u, double v) {
          const auto d = s.first_derivatives(u, v);
          return dot(d[0], cross(d[1], d[2])) / 3.0;
        },
        1e-6, abs_tol);
  }
  return V;
}

/// Solids of a closed shell: MANIFOLD_SOLID_BREP per outer component, BREP_WITH_VOIDS for outer
/// components containing void components (components with inward normals).
std::vector<std::size_t> write_solids(Writer& w, const Shell& input) {
  const Shell sh = input.is_consistently_oriented() ? input : input.oriented();
  ShellWriter sw(w, sh);
  const auto comps = components(sh);
  std::vector<double> vol;
  for (const auto& c : comps) vol.push_back(signed_volume(sh, c));
  const auto closed_shell = [&](const std::vector<std::size_t>& faces, bool flip) {
    std::vector<std::size_t> ids;
    for (std::size_t f : faces) ids.push_back(sw.face(f, flip));
    return w.add("CLOSED_SHELL(''," + list(ids.size(), [&](std::size_t i) { return ref(ids[i]); }) + ")");
  };
  // assign every inward component to the smallest outward component containing it
  constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> owner(comps.size(), none);
  if (comps.size() > 1) {
    std::vector<Shell> outer(comps.size());
    for (std::size_t i = 0; i < comps.size(); ++i) {
      if (vol[i] >= 0.0) continue;
      const Vec3 probe = [&] {
        for (const auto& e : sh.edges())
          if (!e.uses.empty() && std::find(comps[i].begin(), comps[i].end(), e.uses[0].face) != comps[i].end())
            return sh.vertices()[e.v0];
        return sh.faces()[comps[i][0]].surface().control_points()[0];
      }();
      for (std::size_t j = 0; j < comps.size(); ++j) {
        if (vol[j] < 0.0 || (owner[i] != none && vol[j] >= vol[owner[i]])) continue;
        if (outer[j].faces().empty()) {
          std::vector<TrimmedSurface> faces;
          for (std::size_t f : comps[j]) faces.push_back(sh.faces()[f]);
          outer[j] = Shell::sew(std::move(faces));
        }
        if (point_in_shell(outer[j], probe)) owner[i] = j;
      }
    }
  }
  std::vector<std::size_t> solids;
  for (std::size_t j = 0; j < comps.size(); ++j) {
    if (vol[j] < 0.0 && owner[j] != none) continue;  // written as a void below
    // an inward component without container is an inside-out solid: written flipped
    const std::size_t cs = closed_shell(comps[j], vol[j] < 0.0);
    std::vector<std::size_t> voids;
    for (std::size_t i = 0; i < comps.size(); ++i)
      if (owner[i] == j)  // void shell: normals reversed to point out of the void, orientation .F.
        voids.push_back(w.add("ORIENTED_CLOSED_SHELL('',*," + ref(closed_shell(comps[i], true)) + ",.F.)"));
    if (voids.empty()) solids.push_back(w.add("MANIFOLD_SOLID_BREP(''," + ref(cs) + ")"));
    else
      solids.push_back(w.add("BREP_WITH_VOIDS(''," + ref(cs) + "," +
                             list(voids.size(), [&](std::size_t i) { return ref(voids[i]); }) + ")"));
  }
  return solids;
}

std::size_t write_open_shell(Writer& w, const Shell& sh) {
  ShellWriter sw(w, sh);
  std::vector<std::size_t> ids;
  for (std::size_t f = 0; f < sh.faces().size(); ++f) ids.push_back(sw.face(f, false));
  return w.add("OPEN_SHELL(''," + list(ids.size(), [&](std::size_t i) { return ref(ids[i]); }) + ")");
}

std::size_t write_length_unit(Writer& w, StepLengthUnit unit) {
  switch (unit) {
    case StepLengthUnit::Millimetre: return w.add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))");
    case StepLengthUnit::Centimetre: return w.add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.CENTI.,.METRE.))");
    case StepLengthUnit::Metre: return w.add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.))");
    case StepLengthUnit::Inch: {
      const std::size_t mm = w.add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))");
      const std::size_t factor = w.add("LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4)," + ref(mm) + ")");
      const std::size_t dim = w.add("DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)");
      return w.add("(CONVERSION_BASED_UNIT('INCH'," + ref(factor) + ") LENGTH_UNIT() NAMED_UNIT(" + ref(dim) + "))");
    }
  }
  throw GeometryError("write_step: unknown length unit");
}

/// Appends an entity instance, breaking long lines after commas outside strings.
void append_instance(std::string& out, std::size_t id, const std::string& body) {
  std::string line = ref(id) + "=" + body + ";";
  std::size_t col = 0;
  bool in_string = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    out += c;
    ++col;
    if (c == '\'') in_string = !in_string;
    if (c == ',' && !in_string && col >= 76 && i + 1 < line.size()) {
      out += "\n  ";
      col = 2;
    }
  }
  out += '\n';
}

}  // namespace

std::string to_step(const std::vector<Shell>& shells, const StepOptions& options) {
  Writer w;
  // product structure (AP214 conformance class 1 boilerplate)
  const std::size_t app = w.add("APPLICATION_CONTEXT('automotive design')");
  w.add("APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000," + ref(app) + ")");
  const std::size_t pctx = w.add("PRODUCT_CONTEXT(''," + ref(app) + ",'mechanical')");
  const std::size_t product = w.add("PRODUCT(" + str(options.name) + "," + str(options.name) + ",''," + "(" +
                                    ref(pctx) + "))");
  w.add("PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(" + ref(product) + "))");
  const std::size_t pdf = w.add("PRODUCT_DEFINITION_FORMATION('',''," + ref(product) + ")");
  const std::size_t pdctx = w.add("PRODUCT_DEFINITION_CONTEXT('part definition'," + ref(app) + ",'design')");
  const std::size_t pd = w.add("PRODUCT_DEFINITION('design',''," + ref(pdf) + "," + ref(pdctx) + ")");
  const std::size_t pds = w.add("PRODUCT_DEFINITION_SHAPE('',''," + ref(pd) + ")");

  // representation context: units and uncertainty
  const std::size_t len = write_length_unit(w, options.unit);
  const std::size_t ang = w.add("(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))");
  const std::size_t sr = w.add("(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT())");
  const std::size_t unc = w.add("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(" + real(options.uncertainty) + ")," +
                                ref(len) + ",'distance_accuracy_value','confusion accuracy')");
  const std::size_t ctx = w.add("(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((" +
                                ref(unc) + ")) GLOBAL_UNIT_ASSIGNED_CONTEXT((" + ref(len) + "," + ref(ang) + "," +
                                ref(sr) + ")) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and " +
                                "UNCERTAINTY'))");
  const std::size_t origin = point(w, {0, 0, 0});
  const std::size_t axis = w.add("AXIS2_PLACEMENT_3D(''," + ref(origin) + "," + ref(direction(w, {0, 0, 1})) + "," +
                                 ref(direction(w, {1, 0, 0})) + ")");

  // geometry and topology
  std::vector<std::size_t> solids, open_shells;
  for (const auto& sh : shells) {
    if (sh.faces().empty()) continue;
    if (sh.is_closed()) {
      for (std::size_t s : write_solids(w, sh)) solids.push_back(s);
    } else {
      open_shells.push_back(write_open_shell(w, sh));
    }
  }
  const auto items = [&](const std::vector<std::size_t>& ids) {
    std::string s = "(" + ref(axis);
    for (std::size_t id : ids) s += "," + ref(id);
    return s + ")";
  };
  // solids only: ADVANCED_BREP_SHAPE_REPRESENTATION; open shells only:
  // MANIFOLD_SURFACE_SHAPE_REPRESENTATION; both: a plain SHAPE_REPRESENTATION holding both items
  std::vector<std::size_t> rep_items = solids;
  if (!open_shells.empty())
    rep_items.push_back(w.add("SHELL_BASED_SURFACE_MODEL(''," +
                              list(open_shells.size(), [&](std::size_t i) { return ref(open_shells[i]); }) + ")"));
  const char* rep_type = open_shells.empty() ? "ADVANCED_BREP_SHAPE_REPRESENTATION"
                         : solids.empty()    ? "MANIFOLD_SURFACE_SHAPE_REPRESENTATION"
                                             : "SHAPE_REPRESENTATION";
  const std::size_t rep =
      w.add(std::string(rep_type) + "(" + str(options.name) + "," + items(rep_items) + "," + ref(ctx) + ")");
  w.add("SHAPE_DEFINITION_REPRESENTATION(" + ref(pds) + "," + ref(rep) + ")");

  // exchange structure
  const std::string version = "Nurbforge " NURBFORGE_VERSION;
  std::string out;
  out += "ISO-10303-21;\nHEADER;\n";
  out += "FILE_DESCRIPTION(('Nurbforge B-rep model'),'2;1');\n";
  out += "FILE_NAME(" + str(options.name) + "," + str(options.timestamp.empty() ? utc_now() : options.timestamp) +
         ",(" + str(options.author) + "),(" + str(options.organization) + ")," + str(version) + "," + str(version) +
         ",'');\n";
  out += "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n";
  for (std::size_t i = 0; i < w.data().size(); ++i) append_instance(out, i + 1, w.data()[i]);
  out += "ENDSEC;\nEND-ISO-10303-21;\n";
  return out;
}

void write_step(const std::string& path, const std::vector<Shell>& shells, const StepOptions& options) {
  const std::string text = to_step(shells, options);
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("write_step: cannot open '" + path + "' for writing");
  f << text;
  if (!f) throw std::runtime_error("write_step: error writing '" + path + "'");
}

}  // namespace nf
