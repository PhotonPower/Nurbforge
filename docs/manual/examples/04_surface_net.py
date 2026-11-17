"""Surface control nets are stored row-major: index = i * num_v + j, i along u."""
import nurbforge as nf

# a bilinear patch: 2 x 2 control points, degree 1 in both directions
pts = [(0, 0, 0), (0, 1, 0),      # i = 0 (u = 0): j = 0, 1
       (2, 0, 0), (2, 1, 1)]      # i = 1 (u = 1): j = 0, 1
s = nf.NurbsSurface(1, 1, [0, 0, 1, 1], [0, 0, 1, 1], 2, 2, pts)
print("S(1, 0) =", s(1.0, 0.0), " S(0, 1) =", s(0.0, 1.0), " S(1, 1) =", s(1.0, 1.0))
print("net shape", s.num_u, "x", s.num_v, " control_points:", s.control_points.shape)
print("normal at (0.5, 0.5):", s.normal(0.5, 0.5))
su, sv = s.first_derivatives(0.5, 0.5)[1:]
print("S_u =", su, " S_v =", sv)
t = s.transposed()                 # swaps u and v (and flips the normal)
print("transposed normal:", t.normal(0.5, 0.5))
