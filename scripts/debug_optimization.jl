using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise
include(srcdir("HotInverse.jl"))#this line will make all the code available


us_data_dir = srcdir("data","rect_data")
geo_info = init_geoinfo(us_data_dir,pop_dir = "/pop_points/pop_points.csv")
constraint_func = area_constraint

const_val = 10.
const_exp = -1.#the exponent on the area constraint

my_beta = 0.5
mean_distance_exp = x -> mean_distance(x,my_beta)

fitness_over_time,best_ind_over_time,info = sim_anneal((x,y) -> get_fitness(x,y,mean_distance_exp),generate_genome,x -> constraint_func(x,const_val,const_exp),(x,y,z) -> mutate(x,y,num_inds_to_change = z,remove_fac_prob = 0.1,add_fac_prob = 0.1),geo_info,total_generations = 500,delta_t = 5,genome_length = 1)

best_genome = best_ind_over_time[length(best_ind_over_time)]
my_plot =  plot_genome(best_genome,geo_info,title = "$(constraint_func): $(const_val)")

length.(best_ind_over_time)