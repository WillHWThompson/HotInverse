
 function run_evo_model(get_fitness_func,generate_genome_func,geo_info_struct,param_dict)
    """
    run_evo_model(): a wrapper function that runs the evolutionary algorithim <num_runs> number of times and returns the results
    """
     num_runs = pop!(param_dict,:num_runs)
     total_generations = param_dict[:total_generations]
     #algo = pop!(param_dict,:algorithim)
     algo = param_dict[:algorithim]
     println("num_runs: $(num_runs)")
     #println(param_dict)
     #initialize record keeping
     best_ind_over_time_list = [] 
     fitness_over_time_matrix = []
     population_stats_over_time_matrix = []
     #fitness_over_time_matrix = zeros(num_runs,total_generations) 
     #fitness_over_time_matrix =[[(0.0,0.0) for i in 1:total_generations] for i in num_runs] 
     #for each run...
     for i in 1:num_runs
         println("Run: $i")
         #run the evolutionary algorithim
         #algo = param_dict[:algorithim]
         # algo = param_dict[:algorithim]
         # println(algo)

         println("algorithim: $algo")
         if algo == "evolutionary"
             fitness_over_time,best_ind_over_time,population_stats = evolutionary_algorithm(get_fitness_func,generate_genome_func,geo_info_struct ; param_dict...)
         elseif algo == "sim-anneal"
             fitness_over_time,best_ind_over_time,population_stats = sim_anneal(get_fitness_func,generate_genome_func,geo_info_struct ; param_dict...)
         else
             println("algorithim: $algo is not valid. Please choose 'evolutionary' or 'sim-anneal'")
         end

         #keep tack for results
         @show fitness_over_time 
         @show length(fitness_over_time)

         #fitness_over_time_matrix[i,:] = fitness_over_time
         push!(fitness_over_time_matrix,fitness_over_time)
         push!(best_ind_over_time_list,best_ind_over_time)
         push!(population_stats_over_time_matrix,population_stats)
     end
     return fitness_over_time_matrix,best_ind_over_time_list,population_stats_over_time_matrix
 end
