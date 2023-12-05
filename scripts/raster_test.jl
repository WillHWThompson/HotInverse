using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise
include(srcdir("HotInverse.jl"))#this line will make all the code available
using LibGEOS
using LinearAlgebra



us_data_dir = srcdir("data","rect_data")
geo_info = init_geoinfo(us_data_dir,pop_dir = "/pop_points/pop_points.csv")

function generate_pop_points(dist;min::Float64 = 0.,max::Float64 = 1.0,inc::Float64 = 0.1)
    x = min:inc:max
    y = min:inc:max
    #pop tuples
    pop_tuples = [[i,j] for i in x, j in y]
    pop_points = Point2.(pop_tuples)
    pdf_values = map( x -> pdf(dist,x),pop_tuples)
    return vec(pop_points),pdf_values
end

function generate_pop_points_mvGaussian(;μ = [0.5,0.5], Σ = [1.0 0.0;0.0 1.0],kwargs...)
    dist  = MvNormal(μ,Σ)
    return Population(generate_pop_points(dist;kwargs...)...)
end






abstract type population end
struct Population{P<:AbstractArray{Point2{Float64}}} <: population
    pop_points::P
    pdf_values::AbstractArray{Float64}
end

pop_struct = generate_pop_points_mvGaussian()



genome = generate_genome(geo_info.border,10)
polygon = convertShapefileMBRtoRectangle(geo_info.border.MBR)
mbr = generate_MBR_from_points(geo_info.pop_points,delta = 5.0) 

genome = generate_genome(geo_info.border,10)
ind = make_voronoi_individual(genome,get_fitness,geo_info.border)

pop_struct.pop_points



idxs,dist = get_nearest_neighbors(genome,pop_struct.pop_points)#perform neartest neighbors search
df = DataFrame([idxs,dist],["idxs","dist"])
df[!,:pdf] = vec(pop_struct.pdf_values)

df[!,:dist_weighted] = df[!,:dist] .*  vec(pop_struct.pdf_values)









#function create_rasterized_tess(ind::VoronoiIndividual)
tess_poly = Polygon.(ind.tess.Cells)#convert voronoi cells to polygons
rasterized_poly =  map(poly_i -> rasterize!(A,poly_i, fill=1, progress=false),tess_poly)#rasterize voronoi polygons
A = zeros(UInt32, dimz; missingval=UInt32(0))
rasterize!(A,poly_i, fill=i, progress=false) 

plot(A)

Polygon(Point{2, Int}[(3, 1), (4, 4), (2, 4), (1, 2), (3, 1)])







function get_fitness(genome::SVector{Any,Any},pop_struct::Population{Vector{Point2{Float64}}};loss_function::Function = mean_distance)
"""
get_fitness:given an individual and a loss function, compute the fitness of the individual
inputs: 
    genome::SVector - A static vector containing the locations of the facilities 
    pop_points::Vector{Point2{Float64}} - A vector containing the locations of citizens
    loss_function::Function - a the function to calculate the fitness with - takes an a set of arrays and takes out the value of the fitness
output:
    return fitness::Float64 - the fitness value for the individual
"""
    idxs,dist = get_nearest_neighbors(genome,pop_struct.pop_points)#perform neartest neighbors search
    #convert result to dataframe
    df = DataFrame([idxs,dist],["idxs","dist"])
    fitness = loss_function(df,vec(pop_struct.pdf_values))
    return fitness
end



function mean_distance2(nn_df,weights::Vector{Any},beta = 1)
    nn_df[!,:dist_weighted] = nn_df[!,:dist] .* weights 
    grouped = groupby(nn_df,:idxs)
    fac_mean_dist = combine(grouped,:dist_weighted=> mean)#average the mean distance to the facility for each facility
    mean_fac_mean_dist = combine(fac_mean_dist,:dist_mean=>mean).dist_mean_mean[1]#average the mean facility desnity for each facility
    return mean_fac_mean_dist
end


get_fitness(genome,pop_struct;mean_distance2)

typeof(pop_struct)


genome
