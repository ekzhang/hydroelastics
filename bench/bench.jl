using BenchmarkTools
using Hydroelastics
using LinearAlgebra
using Polyhedra
using StaticArrays
import Hydroelastics: intersect_halfplanes, intersect_halfplanes_slow, HalfPlane, Point

macro run_bench(name, expr)
    quote
        printstyled("Benchmarking: ", $name, "\n"; color = :yellow)
        display(@benchmark esc($expr))
        println("\n\n")
    end
end

tet1 = Object(
    Mesh(
        [
            0.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            0.0 0.0 0.0 1.0
        ],
        reshape([1, 2, 3, 4], :, 1),
        [0.0, 0.0, 0.0, 1.0],
    ),
)
tet2 = Object(
    Mesh(
        [
            0.0 1.0 0.0 0.0
            0.0 0.0 1.0 0.0
            1.8 1.8 1.8 0.8
        ],
        reshape([1, 2, 3, 4], :, 1),
        [0.0, 0.0, 0.0, 1.0],
    ),
)

@run_bench "Tet force" begin
    Hydroelastics.tet_force(tet1, tet2, 1, 1)
end

halfplanes = [
    HalfPlane(Point(-0.16469227056401645, -0.2664776914569725), Point(-0.8506508083520399, 0.5257311121191335)),
    HalfPlane(Point(0.028287259696359788, -0.040547756753354335), Point(0.8201447309019166, 0.5721561153861965)),
    HalfPlane(Point(-0.16620458460952442, 0.09547326678562677), Point(-0.49810107091297723, -0.8671189786617204)),
    HalfPlane(Point(0.2970629430715192, -0.23999299562154927), Point(0.6284282897093981, 0.77786752387082)),
    HalfPlane(Point(-0.05197402131154934, 0.06002998491806123), Point(-0.7560128937297079, -0.6545567236797993)),
    HalfPlane(Point(0.1879473775300143, 0.30410524493997143), Point(0.85065080835204, -0.5257311121191336)),
    HalfPlane(Point(0.1619918586413753, -0.0851641599977757), Point(0.46534112719498627, 0.8851314226378474)),
    HalfPlane(Point(0.2970629430715192, -0.23999299562154938), Point(0.6284282897093982, 0.7778675238708199)),
]

@run_bench "Original halfplane algorithm using Polyhedra" begin
    intersect_halfplanes_slow(halfplanes)
end

@run_bench "New halfplane algorithm" begin
    intersect_halfplanes(halfplanes)
end
