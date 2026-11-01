# Benchmarks

`benchmarks/bench.cpp` times representative kernel operations (best of five repetitions,
wall clock) and prints a Markdown table:

```bash
cmake --preset release -DNF_BUILD_BENCHMARKS=ON     # Windows: cmake --preset msvc -DNF_BUILD_BENCHMARKS=ON
cmake --build --preset release --target nf_bench
./build/release/benchmarks/nf_bench                 # --quick: one repetition (ctest smoke test)
```

The CI build (`ci` preset) compiles the benchmarks and runs `nf_bench --quick` as a ctest case,
so they cannot rot; timings on shared CI runners are not meaningful and are not checked.

## Boolean stress test

`benchmarks/stress.cpp` (`nf_stress [cases per family] [filter]`) runs union, intersection and
difference on random but reproducible pairs of seven surface families – sphere/sphere,
sphere/cylinder, sphere/box, cylinder/box, torus/box, torus/sphere, perpendicular cylinders.
Instead of analytic volumes (which do not exist for most pairs) it checks that every result is
closed and consistently oriented and that the volumes satisfy V(A ∪ B) + V(A ∩ B) = V(A) + V(B)
and V(A − B) = V(A) − V(A ∩ B) to 1e-5. ctest runs one case per family (`stress_smoke`); a full
run with 20 cases per family (140 cases, 420 Booleans) passes with a worst identity error of
2.6e-6.

## Reference results

Windows 11, Intel Core i7-11800H (8 cores / 16 threads), MSVC (Visual Studio 2026), Release.
The surface/surface intersection searches its seeds in parallel (see below); everything else
runs on one thread.

| case | per | time [ms] |
|---|---|---:|
| curve evaluation (cubic, 1000 ctrl) | 1e5 points | 13.7 |
| surface evaluation (rational sphere) | 1e5 points | 27.6 |
| closest point on sphere | 1000 points | 188 |
| curve interpolation (1000 points) | 1 fit | 0.21 |
| curve approximation to 1e-4 | 1 fit | 2.4 |
| surface/surface intersection (two spheres) | 1 intersection | 67 (sequential: 167) |
| sewing a cube | 1 shell | 5.4 |
| point in shell (sphere) | 100 points | 300 |
| Boolean intersection (two spheres) | 1 Boolean | 240 (sequential: 351) |
| Boolean difference (drilled cube) | 1 Boolean | 226 |
| Boolean union (coincident boxes) | 1 Boolean | 218 |
| volume of the lens (rel. 1e-6) | 1 volume | 144 |

Evaluation runs at roughly 7 million curve points or 3.6 million rational surface points per
second. The robust global operations (closest point by branch and bound, surface/surface
intersection, ray casting) dominate the Boolean pipeline.

## Parallel seed search in the surface/surface intersection

Profiling the two-sphere intersection showed 90 % of the time in the seed search, almost all of
it in the boundary seeds of the leaf patch pairs (pairs with separated normal cones: every
intersection curve of the pair crosses a patch boundary, found by intersecting the four
boundary iso-curves of each patch with the other patch) – 212 independent leaves of about
0.8 ms each. The subdivision now only collects the leaves; their boundary seeds are computed
on all hardware threads (`std::thread`, work distributed by an atomic counter) and merged at the
positions where the sequential subdivision would have produced them, so the result is
identical to the sequential one. Best of 7 runs:

| intersection | 1 thread | 4 threads | 16 threads |
|---|---:|---:|---:|
| two spheres | 174 ms | 92 ms | 77 ms |
| sphere / cylinder | 64 ms | 38 ms | 36 ms |

Scaling flattens beyond four threads (memory allocation in the small curve/surface
intersections, and the remaining sequential tracing and fitting of about 20 ms).

## Parallelisation (OpenMP) – measured, not adopted

The independent face pairs of a Boolean operation (intersection or coincidence handling per
pair) and the per-face splitting and classification were parallelised with OpenMP as an
experiment. On the machine above it gave no speed-up for the cases in the table: the two-sphere
Boolean has a single face pair, whose surface/surface intersection is about half of the total,
and became slower with 16 threads (≈ 420 ms vs. 330 ms); the coincident boxes (36 pairs)
changed within the noise. The option was therefore not added. Parallelism would pay off inside
the surface/surface intersection (independent seed marching) and for bodies with many faces;
that is left for when such workloads appear.
