using BenchmarkTools
using Hydroelastics
using LinearAlgebra
using Polyhedra
using GLPK
import Hydroelastics: intersect_polygons

lib = DefaultLibrary{Float64}((GLPK.Optimizer))

function slow_intersect(Apoly, Bpoly)
    PA = polyhedron(vrep(Apoly'), lib)
    PB = polyhedron(vrep(Bpoly'), lib)
    res = polyhedron(vrep(intersect(PA, PB)))
    hcat(points(res.vrep)...)
end

polygonA = [
    -1.0 1.0 1.0 -1.0
    -1.0 -1.0 1.0 1.0]
polygonB = [
    0.0 2.0 2.0 0.0
    0.0 0.0 2.0 2.0
]
println("Original Algorithm, non-differentiable GLPK")
@btime slow_intersect(polygonA, polygonB)

println("New Algorithm")
@btime intersect_polygons(polygonA, polygonB)
