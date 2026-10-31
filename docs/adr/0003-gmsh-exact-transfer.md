# ADR 0003 – gmsh integration via Python with exact NURBS transfer

**Context.** Meshing should use the exact geometry. gmsh's OCC kernel accepts rational B-splines
with arbitrary knots.

**Decision.** Integrate gmsh through its Python API in `nurbforge.mesh`, passing degree, distinct
knots, multiplicities and weights unchanged. Linking gmsh into the C++ core is deferred.

**Consequences.** No gmsh build dependency for the kernel; meshes conform to the true geometry.
Tests compare gmsh evaluations with Nurbforge. Pure-C++ users export STEP (planned) or use the
gmsh C++ API themselves.
