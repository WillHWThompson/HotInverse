using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise
include(srcdir("HotInverse.jl"))#this line will make all the code available
##
using LibGEOS

us_data_dir = srcdir("data","rect_data")
geo_info = init_geoinfo(us_data_dir,pop_dir = "/pop_points/pop_points.csv")









shp.points
parts = Int.(ones(length(shp.points)))

# Define an array of points
my_points = [Shapefile.Point(0, 0), Shapefile.Point(0, 1), Shapefile.Point(1, 1), Shapefile.Point(1, 0)]
my_points = [Point2(0, 0), Point2(0, 1), Point2(1, 1), Point2(1, 0)]

# Create a polygon from the points
parts = [0]
polygon = Polygon(my_points, parts)


Point2(0,0)

function print_type_hierarchy(t::Type)
    println(t)
    for st in subtypes(t)
        print_type_hierarchy(st)
    end
end

print_type_hierarchy(Shapefile.Point)




lon, lat = X(0:0.01:1), Y(0:0.01:1)
#ti = Ti(DateTime(2001):Month(1):DateTime(2002))
ras = Raster(rand(lon, lat)) # this generates random numbers with the dimensions given

polygon




genome = generate_genome(geo_info.border,10)
polygon = convertShapefileMBRtoRectangle(geo_info.border.MBR)
mbr = generate_MBR_from_points(geo_info.pop_points,delta = 5.0) 

polygon


function generate_pop_points(dist;min::Float64 = 0.,max::Float64 = 1.0,inc::Float64 = 0.1)
    x = min:inc:max
    y = min:inc:max
    #pop tuples
    pop_tuples = [[i,j] for i in x, j in y]
    pop_points = Point2.(pop_tuples)
    pdf_values = map( x -> pdf(dist,x),pop_tuples)
    return pdf_values,pop_points
end

function generate_pop_points_mvGaussian(;μ = [0.5,0.5], Σ = [1.0 0.0;0.0 1.0],kwargs...)
    dist  = MvNormal(μ,Σ)
    return generate_pop_points(dist;kwargs...)
end

generate_pop_points_mvGaussian()



using Distributions
μ = [0.5,0.5]
Σ = [1.0 0.0; 0.0 1.0]
dist  = MvNormal(μ,Σ)
pdf_values = map( x -> pdf(dist,x),pop_tuples)




genome = generate_genome(geo_info.border,10)
ind = make_voronoi_individual(genome,get_fitness,geo_info.border)
sum(ind.perimeters)
sum(ind.areas)

area_constraint(ind,5.0)

ind.tess


#fitness_over_time,best_ind_over_time,info = sim_anneal(get_fitness,generate_genome,x ->  number_constraint(x,11.),(x,y,z) -> mutate(x,y,num_inds_to_change = z,remove_fac_prob = 0.1,add_fac_prob = 0.1),geo_info,total_generations = 500,delta_t = 5,genome_length = 50)

fitness_over_time,best_ind_over_time,info = sim_anneal(get_fitness,generate_genome,x ->  area_constraint(x,11.),(x,y,z) -> mutate(x,y,num_inds_to_change = z,remove_fac_prob = 0.1,add_fac_prob = 0.1),geo_info,total_generations = 500,delta_t = 5,genome_length = 50)

#fitness_over_time,best_ind_over_time,info = sim_anneal(get_fitness,generate_genome,x -> perimiter_constraint(x,5.0),mutate,geo_info,total_generations = 500,delta_t = 5)

plot(fitness_over_time)

@show length.(best_ind_over_time)





function voronoipopulation(idx_list::Vector{Int64})
    """
    calc_population: given a list of the nearest facility for each citizen in the population, returen the total population of each Voronoi Cell 
    input: 
        idx_list::Vector{Int64}: a list of the nearest facility for each individual in the population
    returns: 
        counts::Vector{Int64}: the total number of individuals in each index
    """
    counts =Int.(zeros(maximum(idx_list)))
    map(x ->  counts[x] = get(counts,x,0)+1,idx_list)
    return counts
end


function calc_perimeter(vertices)
    """
    calc_perimeter: take in a list of veritices of a polygon and returns a the perimiter of the shape
    """
    perimeter = 0.0
    # Loop through each pair of adjacent vertices and add the distance between them to the perimeter
    n = length(vertices)
    @inbounds for i in eachindex(vertices)
        j = mod(i, n) + 1
        dx = vertices[j][1] - vertices[i][1]
        dy = vertices[j][2] - vertices[i][2]
        perimeter += sqrt(dx^2 + dy^2)
    end
    return perimeter
end


function voronoiperimeters(tess::Tessellation{Point2{Float64}})
    """
    get_cell_perimeters: take in a Tesselation object from VoronoiCells.jl and returns the perimiters of each cell in a list
    input: 
        tess::Tesselation, the Voronoi Tesselation object returned by vornoitess()
    """
    map(x ->calc_perimeter(x),tess.Cells)
end



#colors = [:red,:blue,:green,:purple,:yellow]
#idx_color = Dict(i => i%length(colors) for i in collect(1:maximum(idxs)))
#
#
#scatter(geo_info.pop_points,color =[idx_color[i] for i in idxs])
#
#plot!(my_tess)
#
#scatter!(Vector([i for i in genome]),color = [idx_color[i] for i in eachindex(genome) ],markersize = 10)
#annotate!([(genome[i][1] + 0.02, genome[i][2] + 0.03, Plots.text(i)) for i in eachindex(genome)])
#savefig(plotsdir("testing_plots","test.png"))
