### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ a0b3bf62-55a0-11ec-236d-51c7aef190bb
begin
    import Pkg
    Pkg.activate("..")
end

# ╔═╡ d9e6854e-2896-43ed-9d68-9f3af465aef2
begin
    using Revise
    using LinearAlgebra
    using Plots
    using PlutoUI
    using MeshCat
end

# ╔═╡ 341c3b02-7a0f-45b6-8f6c-91cd1fefd889
begin
    function sort_polygon(polygonA::Matrix{Float64})
        # sort 2xN polygon in counter-clockwise order about center of mass
        n = size(polygonA)[2]
        com = [sum(polygonA[i, :]) for i = 1:2] / n
        angles = []  # list of (angle, vtx index corresponding to angle)
        displacements = polygonA - reshape(repeat(com, n), 2, n)
        for ind = 1:n
            disp = displacements[:, ind]
            angle = -atan(disp[2], disp[1]) # negative for counterclockwise
            push!(angles, angle)
        end
        order = sortperm(angles)
        hcat([polygonA[:, order[i]] for i = 1:n]...)::Matrix{Float64}
    end

    function nextidx(idx::Int64, n::Int64)
        if idx == n
            return 1
        else
            return idx + 1
        end
    end

    struct Point
        x::Float64
        y::Float64
    end

    function Base.:+(x::Point, y::Point)
        Point(x.x + y.x, x.y + y.y)
    end

    function Base.:-(x::Point, y::Point)
        Point(x.x - y.x, x.y - y.y)
    end

    function Base.:*(x::Point, y::Float64)
        Point(x.x * y, x.y * y)
    end

    function cross(x::Vector{Float64}, y::Vector{Float64})
        x[1] * y[2] - x[2] * y[1]
    end

    function cross(x::Point, y::Point)
        x.x * y.y - x.y * y.x
    end


    struct HalfPlane
        p::Point
        pq::Point
        angle::Float64

        HalfPlane(p::Point, q::Point) = begin
            p = p
            pq = Point(q.x - p.x, q.y - p.y)
            angle = atan(pq.y, pq.x)
            new(p, pq, angle)
        end
    end

    function out(h::HalfPlane, r::Point)
        return cross(h.pq, r - h.p) > 1e-6
    end

    function Base.isless(x::HalfPlane, y::HalfPlane)
        if (abs(x.angle - y.angle) < 1e-6)
            return cross(x.pq, y.p - x.p) > 0
        end
        return x.angle < y.angle
    end

    function Base.isequal(x::HalfPlane, y::HalfPlane)
        println("is equal called", x.angle, " ", y.angle)
        if (abs(x.angle - y.angle) < 1e-6)
            return true
        else
            return false
        end
    end

    function intersect_halfplanes(s::HalfPlane, t::HalfPlane)
        alpha = cross((t.p - s.p), t.pq) / cross(s.pq, t.pq)
        return s.p + s.pq * alpha
    end

    function cust_unique(s::Vector{HalfPlane})
        res = Vector{HalfPlane}()
        for i = 1:size(s)[1]
            if size(res)[1] > 0
                if abs(res[size(res)[1]].angle - s[i].angle) > 1e-6
                    push!(res, s[i])
                end
            else
                push!(res, s[i])
            end
        end
        res
    end

    function intersect_polygons(polygonA::Matrix{Float64}, polygonB::Matrix{Float64})
        #polygon A is a 2xN matrix, polygon B is a 2xM matrix.
        #returns the intersection of polygon A and polygon B
        #if there is no intersection, returns an empty matrix
        box = [1e7 -1e7 -1e7 1e7
            1e7 1e7 -1e7 -1e7]
        n, m = size(polygonA, 2), size(polygonB, 2)
        polygonA, polygonB = sort_polygon(polygonA), sort_polygon(polygonB)
        box = sort_polygon(box)
        #println("polygon A ", polygonA)
        halfplanes::Array{HalfPlane} = []
        for i = 1:n
            pt1 = Point(polygonA[1, i], polygonA[2, i])
            pt2 = Point(polygonA[1, nextidx(i, n)], polygonA[2, nextidx(i, n)])
            s = HalfPlane(pt1, pt2)
            push!(halfplanes, s)
        end
        for i = 1:m
            t = HalfPlane(Point(polygonB[1, i], polygonB[2, i]), Point(polygonB[1, nextidx(i, m)], polygonB[2, nextidx(i, m)]))
            push!(halfplanes, t)
        end
        for i = 1:4
            t = HalfPlane(Point(box[1, i], box[2, i]), Point(box[1, nextidx(i, 4)], box[2, nextidx(i, 4)]))
            push!(halfplanes, t)
        end
        #println("orig halfplanes", halfplanes)
        sort!(halfplanes)
        #println("sort halfplanes", halfplanes)
        halfplanes = cust_unique(halfplanes)
        #unique!(halfplanes)
        println("unique halfplanes", halfplanes)
        dq = Vector{HalfPlane}()
        len = 0
        for i = 1:length(halfplanes)
            while ((len > 1) && out(halfplanes[i], intersect_halfplanes(dq[len], dq[len-1])))
                pop!(dq)
                len = len - 1
            end

            while ((len > 1) && out(halfplanes[i], intersect_halfplanes(dq[1], dq[2])))
                popfirst!(dq)
                len = len - 1
            end
            push!(dq, halfplanes[i])
            len = len + 1
        end

        while (len > 2 && out(dq[1], intersect_halfplanes(dq[len], dq[len-1])))
            pop!(dq)
            len = len - 1
        end
        while (len > 2 && out(dq[len], intersect_halfplanes(dq[1], dq[2])))
            popfirst!(dq)
            len = len - 1
        end

        if (len < 3)
            return Vector{Point}()
        end
        result = Vector{Point}()
        for i = 1:len
            push!(result, intersect_halfplanes(dq[i], dq[nextidx(i, len)]))
        end
        return result
    end

    function convert_points(points::Vector{Point})
        if size(points)[1] == 0
            return Array{Float64}(undef, 0, 2)
        end
        hcat([[points[i].x; points[i].y] for i = 1:(size(points)[1])]...)::Matrix{Float64}
    end
end

# ╔═╡ 267e9b19-c9f0-45c0-94d6-429c01b81f66
begin
    polygonA = [-1.0 1.0 1.0 -1.0
        -1.0 -1.0 1.0 1.0]
    polygonB = [0.0 2.0 2.0 0.0
        0.0 0.0 2.0 2.0]
    total_res = convert_points(intersect_polygons(polygonA, polygonB))
    with_terminal() do
        println(total_res)
        println(length(total_res))
    end
    #intersect_polygons(polygonA, polygonB)
end

# ╔═╡ c5f6f3e6-e1f7-4f1d-a4bd-ea35e0716fa7

begin
    polygonA2 = [-3.0 -1.0 -1.0 -3.0
        -3.0 -3.0 -1.0 -1.0]
    polygonB2 = [0.0 2.0 2.0 0.0
        0.0 0.0 2.0 2.0]
    total_res2 = convert_points(intersect_polygons(polygonA2, polygonB2))
    with_terminal() do
        println(total_res2)
        println(length(total_res2))
    end
    #intersect_polygons(polygonA, polygonB)
end

# ╔═╡ a0895905-5a34-40db-adb8-20c85bf0ede0

begin
    polygonA3 = [-3.0 2.0 2.0 -3.0
        -3.0 -3.0 2.0 2.0]
    polygonB3 = [0.0 -1.0 -1.0 0.0
        0.0 0.0 -1.0 -1.0]
    total_res3 = convert_points(intersect_polygons(polygonA3, polygonB3))
    with_terminal() do
        println(total_res3)
        println(length(total_res3))
    end
    #intersect_polygons(polygonA, polygonB)
end

# ╔═╡ 00d3d951-d6cc-474e-b03b-71057e9e977a

begin
    polygonA4 = [-1.0 1.0 0.0
        0.0 0.0 3.0]
    polygonB4 = [-1.0 1.0 0.0
        2.0 2.0 1.0]
    total_res4 = convert_points(intersect_polygons(polygonA4, polygonB4))
    with_terminal() do
        println(total_res4)
        println(length(total_res4))
    end
    #intersect_polygons(polygonA, polygonB)
end

# ╔═╡ f945f54a-851c-442f-aea8-b5a7fd4f837e
begin
    polygonB5 = [-3.0 2.0 2.0 -3.0
        -3.0 -3.0 2.0 2.0]
    polygonA5 = [0.0 -1.0 -1.0 0.0
        0.0 0.0 -1.0 -1.0]
    total_res5 = convert_points(intersect_polygons(polygonA5, polygonB5))
    with_terminal() do
        println(total_res5)
        println(length(total_res5))
    end
    #intersect_polygons(polygonA, polygonB)
end

# ╔═╡ Cell order:
# ╠═a0b3bf62-55a0-11ec-236d-51c7aef190bb
# ╠═d9e6854e-2896-43ed-9d68-9f3af465aef2
# ╠═341c3b02-7a0f-45b6-8f6c-91cd1fefd889
# ╠═267e9b19-c9f0-45c0-94d6-429c01b81f66
# ╠═c5f6f3e6-e1f7-4f1d-a4bd-ea35e0716fa7
# ╠═a0895905-5a34-40db-adb8-20c85bf0ede0
# ╠═00d3d951-d6cc-474e-b03b-71057e9e977a
# ╠═f945f54a-851c-442f-aea8-b5a7fd4f837e
