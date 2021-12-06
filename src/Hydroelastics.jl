module Hydroelastics

using RigidBodyDynamics
using LinearAlgebra
using Statistics
using StaticArrays
using GeometryBasics, Meshes
using CoordinateTransformations, MeshCat
using Polyhedra
using GLPK


include("object.jl")
include("shapes.jl")
include("rendering.jl")
include("forces.jl")

end
