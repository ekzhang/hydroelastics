module Hydroelastics

using RigidBodyDynamics
using LinearAlgebra
using Statistics
using StaticArrays
using GeometryBasics, Meshes, CoordinateTransformations
using Polyhedra, Polyhedra.Core
using DiffOpt, GLPK

include("object.jl")
include("shapes.jl")

end
