using Revise
Revise.revise()
using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise
include(srcdir("HotInverse.jl"))#this line will make all the code available
#using LibGEOS
using LinearAlgebra
using GeometryBasics
using Rasters

##
dist = MvNormal([0.5,0.5],[0.1 0.0;0.0 0.1])
geo_info = init_geoinfo(1,1,100,dist)
##
generate_genome(geo_info.population,geo_info,n_facs = 10)
##
my_beta = 1
my_n_facs = 10
my_num_inds_to_change = 2
#fitness closures
#dist_metric_cloj = (x,y) -> euclidean_exp(x,y;beta = my_beta) 
#fitness_function = (x,y,z) -> mean_dist(x,y,z,dist_metric_cloj)
#mean_distance_exp = x -> mean_distance(x,beta_val)#closure for mean fitness function with correct exponent
#constraint closures
constraint_func = area_constraint
const_val = 10.
alpha_val = 1.
constraint_closure = x -> constraint_func(x,const_val,alpha_val)#create closure with constraint value and exponent

generate_genome_cloj = x -> generate_genome(geo_info.population,geo_info,n_facs = x)
generate_genome_function() = generate_genome(geo_info.population,geo_info,n_facs = my_n_facs)
generate_ind_function() =  make_voronoi_individual(generate_genome_function,fitness_function,geo_info)

dist_metric_cloj = (x,y) -> euclidean_exp(x,y;beta = my_beta) 
fitness_function = (x,y,z) -> -1*mean_dist(x,y,z,dist_metric_cloj)
instant_ind_function = (x) -> make_voronoi_individual(x,fitness_function,geo_info)
mutation_function =(x) -> mutate(x,geo_info,instant_ind_function,num_inds_to_change = my_num_inds_to_change,remove_fac_prob = 0.1,add_fac_prob = 0.1)

#constraint_function = x -> constraint_func(x,const_val,alpha_val)#create closure with constraint value and exponent

##
ind = generate_ind_function()
##
mutation_function(ind)
##

#fitness_over_time,best_ind_over_time,info = sim_anneal((x,y) -> get_fitness(x,y,mean_distance_exp),,constaint_closure,(x,y,z) -> mutate(x,y,num_inds_to_change = z,remove_fac_prob = 0.1,add_fac_prob = 0.1),geo_info,total_generations = 500,delta_t = 5,genome_length = 1)
fitness_over_time,best_ind_over_time,info = sim_anneal(fitness_function,
                                                       generate_ind_function,
                                                       constraint_closure,
                                                       mutation_function,
                                                       total_generations = 500,
                                                       delta_t = 5)

##

my_savename = datadir("rasters","io_test",savename("test","jld2")) 


##
best_ind_i = best_ind_over_time[end]#
my_plot = plot_population_with_tess(best_ind_i,fitness_function,geo_info)
savefig(my_plot)
genome = best_ind_over_time[end]
v_ind = make_voronoi_individual(genome,fitness_function,geo_info)
mbr_rect = convertMBRtoRectangle(geo_info.MBR)
my_tess = voronoicells(Vector([i for i in genome]),mbr_rect)
tess_poly = GeometryBasics.Polygon.(my_tess.Cells)
tess_mp = MultiPolygon(tess_poly)
ras = geo_info.population.pop_points
ras_tess = create_rasterized_tess(tess_mp,ras)
ras_tess = ras_tess .|> i -> replace_missing(i,0.0)
ras_tess_collected = collect_voronoi_raster(ras_tess)
heatmap(ras_tess_collected)
#
#function plot_population_with_tess(genome,fitness_function,geo_info)
#    v_ind = make_voronoi_individual(genome,fitness_function,geo_info)
#    mbr_rect = convertMBRtoRectangle(geo_info.MBR)
#    my_tess = voronoicells(Vector([i for i in genome]),mbr_rect)
#    tess_poly = GeometryBasics.Polygon.(my_tess.Cells)
#    tess_mp = MultiPolygon(tess_poly)
#    heatmap(geo_info.population.pop_points)
#    scatter!(genome,color = :white)
#    plot!(my_tess,color = :white)
#end
##
final_genome = best_ind_over_time[end]
plot_voronoi_raster(final_genome,fitness_function,geo_info)
##
best_ind_i = best_ind_over_time[end]#
my_plot = plot_population_with_tess(best_ind_i,fitness_function,geo_info)
##
