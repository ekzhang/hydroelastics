module Hydroelastics

using RigidBodyDynamics
using LinearAlgebra
using Statistics
using StaticArrays
using GeometryBasics, Meshes
using CoordinateTransformations, MeshCat
using Polyhedra, Polyhedra.Core
using DiffOpt, GLPK

include("object.jl")
include("shapes.jl")
include("rendering.jl")
include("forces.jl")

end
