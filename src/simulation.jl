const World = Dict{Symbol,Object}

"""
Simulate a world of objects for some number of iterations and fixed timestep.

This repeatedly applies the hydroelastic force model to each pair of objects in
the scene, returning a vector of how the objects change over time.

The kinematics are implemented using a symplectic Euler method, which is more
numerically stable than forward Euler integration.
"""
function simulate(
    world::World;
    dt::Float64 = 1e-3,
    num_steps::Int64 = 1000,
    forces::Function = (a, α) -> (),
)::Dict{Symbol,Vector{Object}}
    ids = collect(keys(world))
    objects = Dict(id => [world[id]] for id in ids)

    for step = 1:num_steps
        # Hydroelastic computations for each pair of objects
        a = Dict(id => @SVector zeros(3) for id in ids)
        α = Dict(id => @SVector zeros(3) for id in ids)
        for i = 1:length(ids)
            for j = i+1:length(ids)
                u, v = ids[i], ids[j]
                results = compute_force(last(objects[u]), last(objects[v]))
                a[u] += 100.0 * results.F_AB
                a[v] += 100.0 * results.F_BA
                α[u] += 100.0 * results.τ_AB
                α[v] += 100.0 * results.τ_BA
            end
        end

        # Apply custom user-specified forces
        forces(a, α)

        # Semi-implicit Euler integration
        for id in ids
            obj = last(objects[id])
            obj = Object(obj.mesh, obj.pose, obj.v .+ dt .* a[id], obj.ω .+ dt .* α[id])
            obj = translate(obj, dt * obj.v)
            obj = changeRotation(obj, dt * obj.ω)
            push!(objects[id], obj)
        end
    end

    objects
end

export World, simulate
