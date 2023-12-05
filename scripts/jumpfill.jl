
using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise
include(srcdir("HotInverse.jl"))#this line will make all the code available
using LinearAlgebra


function jump_flood_voronoi(points::AbstractArray{T}, width::Int, height::Int) where T<:AbstractVector
    # Initialize the distance and index arrays
    @show T
    dist = Inf .* ones(Float32, width, height)
    index = zeros(Int, width, height)

    # Compute the maximum distance
    max_dist = norm([width, height])

    # Initialize the seed points
    seeds = [Int.(floor.(p)) for p in points]

    # Assign the seed points to the nearest pixels
    for i in 1:length(points)
        x, y = seeds[i] .+ 1
        @show i
        dist[x, y] = norm(points[i] .- seeds[i])
        index[x, y] = i
    end

    # Propagate the distances and indices
    for d in 1:floor(Int, log2(max_dist))
        step = 2^(d-1)
        for x in step:step:width
            for y in step:step:height
                # Find the minimum distance and index in the neighborhood
                min_dist = Inf
                min_index = 0
                for dx in -step:step
                    for dy in -step:step
                        nx, ny = x + dx, y + dy
                        if nx >= 1 && nx <= width && ny >= 1 && ny <= height
                            d = dist[nx, ny] + norm([dx, dy])
                            if d < min_dist
                                min_dist = d
                                min_index = index[nx, ny]
                            end
                        end
                    end
                end

                # Update the distance and index at the current pixel
                dist[x, y] = min_dist
                index[x, y] = min_index
            end
        end
    end

    # Compute the Voronoi diagram
    voronoi = Dict()
    for i in 1:length(points)
        voronoi[i] = []
    end

    res = zeros(Int, width, height)
    for x in 1:width
        for y in 1:height
            i = index[x, y]
            if i > 0
                p = points[i]
                v = [x-1, y-1]
                push!(voronoi[i], v)
                res[x, y] = i
            end
        end
    end

    return voronoi,res
end


a = [[30,10],[10,20],]

a = [randin]

points = [[rand(0:100), rand(0:100)] for i in 1:100]

voronoi_d,res = jump_flood_voronoi(a,100,100)

heatmap(res)

rev_d = Dict(v => k for (k,v) in voronoi_d)

new_d = Dict()
for (k, v) in d
    for coord in v
        new_d[coord] = k
    end
end


x = [i[1] for i in a]
y = [i[2] for i in a]

heatmap(index)
scatter!(x,y,seriestype = :scatter)



# Define a dictionary that maps pixels to categorical variables
data = rev_d
# Convert the dictionary to a matrix
matrix = [data[(x, y)] for x in 1:10, y in 1:10]

# Define the color map for the categories
colors = Dict("A" => :red, "B" => :green, "C" => :blue)

# Plot the heatmap
heatmap(matrix, color = [colors[x] for x in matrix], xticks = false, yticks = false)rev_d



