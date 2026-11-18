"""Fillets and chamfers of solid edges: a box edge, the rim of a cylinder and the base of a boss."""
import math

import nurbforge as nf

# --- one box edge, rounded and chamfered
box = nf.make_box((0, 0, 0), (4, 3, 2))
edge = nf.edges_between(box, "top", "front")
r = 0.5
rounded = nf.fillet_edges(box, edge, r)
print(f"box edge {edge}: fillet r = {r}: {len(rounded.faces)} faces, closed={rounded.is_closed}, "
      f"volume {rounded.volume():.8f} (24 - 4 r^2 (1 - pi/4) = {24 - 4 * r * r * (1 - math.pi / 4):.8f})")
bevel = nf.chamfer_edges(box, edge, r)
print(f"chamfer d = {r}: volume {bevel.volume():.8f} (24 - 4 d^2/2 = 23.5), "
      f"labels {sorted(set(bevel.face_labels))}")

# --- a closed edge: the rim of a cylinder (one fillet face around)
cyl = nf.make_cylinder((0, 0, 0), (0, 0, 1), 2.0, 3.0)
rim = nf.fillet_edges(cyl, nf.edges_between(cyl, "top", "side"), 0.4)
fillet = rim.faces[rim.faces_labeled("fillet")[0]]
print(f"cylinder rim: {len(rim.faces)} faces, fillet face area {fillet.area():.6f}")

# --- a concave edge adds material: the base of a boss on a plate
plate = nf.make_box((-3, -3, -1), (6, 6, 1))
boss = nf.boolean(plate, nf.make_cylinder((0, 0, -0.5), (0, 0, 1), 1.0, 2.5),
                  nf.BooleanOp.Union)
base = [i for i in nf.edges_between(boss, "top", "side")
        if abs(boss.edges[i].curve(boss.edges[i].curve.domain[0]).z) < 1e-9]
strong = nf.fillet_edges(boss, base, 0.3)
print(f"boss base fillet: volume {boss.volume():.6f} -> {strong.volume():.6f} (material added)")
