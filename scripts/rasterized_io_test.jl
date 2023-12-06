@everywhere using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise

src_path = srcdir("HotInverse.jl")
@everywhere src_path
@everywhere include($src_path)#this line will make all the code available

##
#debug_logger = ConsoleLogger(stderr, Logging.Debug)
debug_logger = ConsoleLogger(stderr, Logging.Info)
global_logger(debug_logger);

#create geo info
dist = MvNormal([0.5,0.5],[0.1 0.0;0.0 0.1])
geo_info = init_geoinfo(1,1,100,dist)
generate_genome(geo_info.population,geo_info,n_facs = 10)
#model parameters
my_beta = 1
my_n_facs = 10
my_num_inds_to_change = 2

#fitness closures
dist_metric_cloj = (x,y) -> euclidean_exp(x,y;beta = my_beta) 
fitness_function = (x,y,z) -> mean_dist(x,y,z,dist_metric_cloj)
mean_distance_exp = x -> mean_distance(x,beta_val)#closure for mean fitness function with correct exponent

#create all closures
generate_genome_cloj = x -> generate_genome(geo_info.population,geo_info,n_facs = x)
generate_genome_function() = generate_genome(geo_info.population,geo_info,n_facs = my_n_facs)
generate_ind_function() =  make_voronoi_individual(generate_genome_function,fitness_function,geo_info)

dist_metric_cloj = (x,y) -> euclidean_exp(x,y;beta = my_beta) 
fitness_function = (x,y,z) -> mean_dist(x,y,z,dist_metric_cloj)
instant_ind_function = (x) -> make_voronoi_individual(x,fitness_function,geo_info)
mutation_function =(x) -> mutate(x,geo_info,instant_ind_function,num_inds_to_change = my_num_inds_to_change,remove_fac_prob = 0.1,add_fac_prob = 0.1)


#constraint_funcs = [area_constraint,perimeter_constraint]
constraint_funcs = [area_constraint]
constraint_range = [10.]
beta_value_range = [1.]
alpha_value_range = [1.]


all_params = Dict(
    "constraint" => constraint_funcs,
    "const_val" => constraint_range,
    "beta_val" => beta_value_range,#cost benefit constraint
    "alpha_val" => alpha_value_range
    )

dicts = dict_list(all_params)


function convert_point_to_vec(ind)
    return Vector.(ind)
end


function make_sim(d::Dict)

    @unpack constraint,const_val,beta_val,alpha_val = d
    println("running simulation...")

    constraint_closure = x -> constraint(x,const_val,alpha_val)#create closure with constraint value and exponent

    fitness_over_time,best_ind_over_time,info = sim_anneal(fitness_function,
                                                           generate_ind_function,
                                                           constraint_closure,
                                                           mutation_function,
                                                           total_generations = 500,
                                                           delta_t = 5)
    println("finished sim")
    fulld = copy(d)
    fulld["fitness_over_time"] = fitness_over_time
    fulld["best_ind_over_time"] = best_ind_over_time
    println("saving data")
   # push!(best_ind_by_const,best_ind_over_time[length(best_ind_over_time)])
   # push!(fitness_over_time_by_const,fitness_over_time)
   # push!(genome_length_over_time_by_const,length.(best_ind_over_time))
   return fulld

end



f_list = []

#generate data from sims
Threads.@threads for (i,d) in enumerate(dicts)
	@show d
	my_d = deepcopy(d)
	f = make_sim(my_d)
	push!(f_list,f)
end

##

#save sims
Threads.@threads for (i,f) ∈ enumerate(f_list)
	println("Running on thread $(Threads.threadid)")
    f["constraint"] = String(Symbol(f["constraint"]))
    my_savename = datadir("rasters","io_test",savename(f,"jld2")) 
    @show my_savename
    wsave(my_savename,f)
    #wsave(my_savename,f)
    my_savename_pdf = replace(my_savename,"jld2" => "svg")
    best_ind_i = f["best_ind_over_time"][length(f["best_ind_over_time"])]
    #my_plot = plot_population_with_tess(best_ind_i,fitness_function,geo_info,title = "$(f["constraint"]): $(round(f["const_val"])), beta: $(f["beta_val"]), alpha: $(f["alpha_val"])")
    #savefig(my_plot,my_savename_pdf)
end

##


#best_ind_i = Point2.(last(f["best_ind_over_time"]))
#best_ind_i = f["best_ind_over_time"][length(f["best_ind_over_time"])]
 #my_plot = plot_genome(best_ind_i,geo_info,title = "$(f["constraint"]): $(f["const_val"])")



# constraint = perimeter_constraint
# const_val = 20
# fitness_over_time,best_ind_over_time,info = sim_anneal(get_fitness,generate_genome,x -> constraint(x,const_val),(x,y,z) -> mutate(x,y,num_inds_to_change = z,remove_fac_prob = 0.1,add_fac_prob = 0.1),geo_info,total_generations = 500,delta_t = 5,genome_length = 50)
# best_genome = best_ind_over_time[length(best_ind_over_time)]
# my_plot =  plot_genome(best_genome,geo_info,title = "$(constraint): $(const_val)")
# my_savename_pdf = datadir("constraints","layout_$(constraint)_$(const_val).pdf")
# savefig(my_plot,my_savename_pdf)
