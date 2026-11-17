"""The NurbsSurface class: control net, evaluation, derivatives, iso-curves and
shape-preserving operations in u and v."""
import math

import numpy as np
import nurbforge as nf

# a bicubic patch over [0, 3] x [0, 2] with a bump: 4 x 4 control points, row-major
nu, nv = 4, 4
pts = np.array([[3 * i / 3, 2 * j / 3, 1.0 if i in (1, 2) and j in (1, 2) else 0.0]
                for i in range(nu) for j in range(nv)])  # index i * nv + j, i along u
k = [0, 0, 0, 0, 1, 1, 1, 1]
s = nf.NurbsSurface(3, 3, k, k, nu, nv, pts)
print(f"degrees ({s.degree_u}, {s.degree_v}), net {s.num_u} x {s.num_v}, "
      f"domain {s.domain_u} x {s.domain_v}")
print("S(0.5, 0.5) =", s(0.5, 0.5))
p, su, sv = s.first_derivatives(0.25, 0.5)
print("S_u =", su, " S_v =", sv)
d = s.derivatives(0.25, 0.5, 2)  # d[k][l] = d^(k+l) S / du^k dv^l
print("S_uv =", d[1][1], " S_vv =", d[0][2])
print("unit normal", s.normal(0.25, 0.5))
grid = s.sample(5, 3)
print("sample shape", grid.shape, "corner", grid[4, 2])

# iso-curves: u = const runs along v, v = const along u
cu, cv = s.iso_curve_u(0.5), s.iso_curve_v(0.5)
print(f"iso_curve_u(0.5)(0.5) = {cu(0.5)}, iso_curve_v(0.5)(0.5) = {cv(0.5)}")

# knot refinement and degree elevation keep the shape
r = s.refine_knots_u([0.25, 0.5, 0.75]).elevate_degree_v(1)
dev = max(nf.distance(s(u, v), r(u, v)) for u in np.linspace(0, 1, 11) for v in np.linspace(0, 1, 11))
print(f"refined/elevated: net {r.num_u} x {r.num_v}, degree_v {r.degree_v}, deviation {dev:.1e}")
back, removed = r.try_remove_knot_u(0.5, 1)
print("knot u = 0.5 removed", removed, "time(s); degree_v back to",
      r.reduce_degree_v(1e-10).degree_v)

# splitting, sub-surfaces, Bezier patches, transposition
a, b = s.split_u(0.4)
print("split_u domains", a.domain_u, b.domain_u)
print("sub_surface domain", s.sub_surface(0.2, 0.6, 0.1, 0.9).domain_v)
print("Bezier patches of the refined surface:", len(r.decompose_bezier()))
t = s.transposed()
print("transposed: T(v, u) == S(u, v):", nf.distance(t(0.5, 0.25), s(0.25, 0.5)) < 1e-14)
joined = nf.join_u(a, b)
dev = max(nf.distance(s(u, v), joined(u, v)) for u in np.linspace(0, 1, 11)
          for v in np.linspace(0, 1, 11))
print(f"join_u of the halves: deviation {dev:.1e}")

# curvature of a sphere of radius 2: k1 = k2 = -1/2 with the outward normal
meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
kc = nf.surface_curvature(sphere, 0.3, 0.4)
print(f"sphere r = 2: K = {kc.gaussian:.10f}, H = {kc.mean:.10f}, k1 = {kc.k1:.10f}, "
      f"k2 = {kc.k2:.10f}")
print(f"area of the upper half: {nf.area(sphere, 0, 1, 0.5, 1):.10f} (8 pi = {8 * math.pi:.10f})")

# closest point on a surface
cp = nf.closest_point(sphere, (3, 0, 4))
print(f"closest point: (u, v) = ({cp.u:.6f}, {cp.v:.6f}), point {cp.point}, "
      f"distance {cp.distance:.6f}")

# integrate: adaptive Gauss-Kronrod quadrature of any Python function
print("integral of sin over [0, pi] =", nf.integrate(math.sin, 0.0, math.pi))
print("integral of sqrt(|x|) over [-1, 1] with a break at 0 =",
      nf.integrate(lambda x: math.sqrt(abs(x)), -1.0, 1.0, breaks=[0.0]))
