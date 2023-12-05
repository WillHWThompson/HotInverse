using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise
include(srcdir("HotInverse.jl"))#this line will make all the code available


us_data_dir = srcdir("data","rect_data")
geo_info = init_geoinfo(us_data_dir,pop_dir = "/pop_points/pop_points.csv")

area_const_range = 1:20.0:20
perimiter_const_range = 1:20.0:20.0

#concatenate the two lists
# const_val_list = vcat(area_const_range,perimiter_const_range)
# const_function_list = vcat(area_const_list,perimiter_const_list)


# best_ind_by_const = []
# fitness_over_time_by_const = []
# genome_length_over_time_by_const = []


# all_params_dict = Dict(
#     "constraint_func" => const_function_list,
#     "const_val_list" => const_val_list
# )


#constraint_funcs = [area_constraint,perimeter_constraint]
constraint_funcs = [area_constraint,perimeter_constraint,number_constraint]
constraint_range = LinRange(1,40,4) |> collect 
beta_value_range = [0,0.5,1]
alpha_value_range = [-1,0.5,1]


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

    mean_distance_exp = x -> mean_distance(x,beta_val)#closure for mean fitness function with correct exponent
    constaint_closure = x -> constraint_func(x,const_val,alpha_val)#create closure with constraint value and exponent

    fitness_over_time,best_ind_over_time,info = sim_anneal((x,y) -> get_fitness(x,y,mean_distance_exp),generate_genome,constaint_closure,(x,y,z) -> mutate(x,y,num_inds_to_change = z,remove_fac_prob = 0.1,add_fac_prob = 0.1),geo_info,total_generations = 500,delta_t = 5,genome_length = 1)
    println("finished sim")
    fulld = copy(d)
    #println("copied directory")
    fulld["fitness_over_time"] = fitness_over_time
    #fulld["best_ind_over_time"] = convert_point_to_vec.(best_ind_over_time)
    fulld["best_ind_over_time"] = best_ind_over_time
    #println("saving data")
   # push!(best_ind_by_const,best_ind_over_time[length(best_ind_over_time)])
   # push!(fitness_over_time_by_const,fitness_over_time)
   # push!(genome_length_over_time_by_const,length.(best_ind_over_time))
   return fulld

end



f_list = []

#generate data from sims
for (i,d) in enumerate(dicts)
my_d = deepcopy(d)
f = make_sim(my_d)
push!(f_list,f)
end


#save sims
for (i,f) ∈ enumerate(f_list)
    f["constraint"] = String(Symbol(f["constraint"]))
    my_savename = datadir("constraints","constraint_test2",savename(f,"jld2")) 
    @show my_savename
    wsave(my_savename,f)
    #wsave(my_savename,f)
    my_savename_pdf = replace(my_savename,"jld2" => "svg")
    best_ind_i = f["best_ind_over_time"][length(f["best_ind_over_time"])]
    my_plot = plot_genome(best_ind_i,geo_info,title = "$(f["constraint"]): $(round(f["const_val"])), beta: $(f["beta_val"]), alpha: $(f["alpha_val"])")
    savefig(my_plot,my_savename_pdf)
end


f_list[1]


f = f_list[1]
#best_ind_i = Point2.(last(f["best_ind_over_time"]))
best_ind_i = f["best_ind_over_time"][length(f["best_ind_over_time"])]
 my_plot = plot_genome(best_ind_i,geo_info,title = "$(f["constraint"]): $(f["const_val"])")








# constraint = perimeter_constraint
# const_val = 20
# fitness_over_time,best_ind_over_time,info = sim_anneal(get_fitness,generate_genome,x -> constraint(x,const_val),(x,y,z) -> mutate(x,y,num_inds_to_change = z,remove_fac_prob = 0.1,add_fac_prob = 0.1),geo_info,total_generations = 500,delta_t = 5,genome_length = 50)
# best_genome = best_ind_over_time[length(best_ind_over_time)]
# my_plot =  plot_genome(best_genome,geo_info,title = "$(constraint): $(const_val)")
# my_savename_pdf = datadir("constraints","layout_$(constraint)_$(const_val).pdf")
# savefig(my_plot,my_savename_pdf)
