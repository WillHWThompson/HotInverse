#imports
using StatsBase
import Base.Threads.@threads

#setup
rng = Random.MersenneTwister() #unseeded

#function genome_length(i::Individual)
#    length(i.genome)
#end

#function evolutionary_algorithm(fitness_function::Function,generate_genome_function::Function,geo_info;
#                                num_elements_to_mutate=10,total_generations=100, 
#                                num_parents=10, num_children=10, crossover_points=1,genome_length = 500,
#                                tournament_size=4, num_tournament_winners::Integer=3,elitisim=false,
#                                num_perturbations = 0,perturbation_size=1,perturbation_distribution = Pareto(),kwargs...)
#
#
function evolutionary_algorithm(fitness_function::Function,generate_genome_function::Function,geo_info;
                                selection_function::Function=rank_assignment,
                                num_elements_to_mutate=10,total_generations=100, 
                                num_parents=10, num_children=10, crossover_points=1,genome_length = 500,elitisim=false,
                                num_perturbations = 0,perturbation_size=1,perturbation_distribution = Pareto(),kwargs...)

    """ Evolutinary Algorithm

    parameters: 
    fitness_function: function that returns the fitness of a genome given the genome as an input parameter
    gen_individual_function: function returns a new individual in the population
    geo_info<struct> contains geo_info.border which specifies the border in which to place points and 
                        geo_info.population which specifies the population as points
    total_generations: (int) number of total iterations for stopping condition
    num_parents: (int) the number of parents we downselect to at each generation (mu)
    num_children: (int) the number of children (note: parents not included in this count) that we baloon to each generation (lambda)
    genome_length: (int) length of the genome to be evoloved
    num_elements_to_mutate: (int) number of alleles to modify during mutation (0 = no mutation)
    perturbation_size: (int) size of perturbation
    num_perturbations: (int) number of perturbations applied per generation
    perturbation_distribution: (distribution) to use in perturbation gen_catastrophe_pos, default pareto

    returns:
    fitness_over_time: (numpy array) track record of the top fitness value at each generation
    best_ind_over_time: highest fitness individual
    """
    # initialize record keeping
    solution = nothing # best genome so far
    solution_fitness = -999 # fitness of best genome so far
    solution_fitness_diff = -999  # fitness of best genome so far
    #solution_generation = 0 # time (generations) when solution was found
    #solutions_over_time = zeros((total_generations, genome_length),dtype=int)
    #diversity_over_time = zeros(total_generations)
    #fitness_over_time = zeros(total_generations)
    #@infiltrate

    fitness_over_time = [(0.0,0.0) for i in 1:total_generations]
    best_ind_over_time = [] 
    best_ind = nothing 
    population_stats = []

    # if cmp(perturbation_distribution, "pareto")==0
    #     perturbation_distribution = Pareto()
    # if cmp(perturbation_distribution, "uniform")==0
    #     perturbation_distribution = Uniform()
    if perturbation_distribution == "uniform"
        perturbation_distribution = Distributions.Uniform()
    elseif perturbation_distribution == "pareto"
        perturbation_distribution = Distributions.Pareto()
    elseif perturbation_distribution == "normal"
        perturbation_distribution = Distributions.Normal()
    elseif perturbation_distribution == "arcsine"
        perturbation_distribution = Arcsine()
    elseif perturbation_distribution == "exponential"
        perturbation_distribution = Exponential()
    else
        println("perturbation dist does not match list of: 'pareto','normal','arcsine','exponential')")
        perturbation_distribution = "normal"
    end

    population = []
    for i in 0:num_parents # only create parents for initialization (the mu in mu+lambda)
        new_genome = generate_genome_function(geo_info.border,genome_length)
        new_fitness = fitness_function(new_genome,geo_info.pop_points)
        new_ind = Individual(new_fitness,-999,new_genome)
        push!(population, new_ind)
    end

    best_ind = population[rand(1:end)].genome

    for generation in 1:total_generations
     #for generation in 1:total_generations

       new_children = []
       for _ in 1:Int(round(num_children/2))#multithreaded fitness evaluation
            # inheritance
            parent1 = population[rand(1:end)] #random.choice(population, size=2) # pick 2 random parents
            parent2 = population[rand(1:end)] #random.choice(population, size=2) # pick 2 random parents
            
            #do crossover            
            if crossover_points > 0
                #perform the crossover on the two parent genomes
                new_genome1,new_genome2 = crossover(parent1.genome,parent2.genome,crossover_points) 
            else
               #if crossover if false, set the new genomes equal to the parent genomes
                new_genome1,new_genome2 = parent1.genome,parent2.genome
            end

            if num_elements_to_mutate > 0  
                new_genome1 = mutate(new_genome1,geo_info.border,num_elements_to_mutate)
                new_genome2 = mutate(new_genome2,geo_info.border,num_elements_to_mutate)
            end


            fit_diff = 0 #difference it fitness before and after

            #evaluate the original fitness
            fitness_orig_1 = fitness_function(new_genome1,geo_info.pop_points)
            fitness_orig_2 = fitness_function(new_genome2,geo_info.pop_points)


            #@show length(new_genome1)
            

            #create copy of genoem to perturb
            new_genome_perturbed1 = new_genome1
            new_genome_perturbed2 = new_genome2

            if num_perturbations > 0
            # do perturbation
                for i in 1:num_perturbations
                    new_genome_perturbed1 = perturb(new_genome_perturbed1, geo_info, perturbation_size, perturbation_distribution)
                    new_genome_perturbed2 = perturb(new_genome_perturbed2, geo_info, perturbation_size, perturbation_distribution)
                end
            else
                new_genome_perturbed1 = new_genome1
                new_genome_perturbed2 = new_genome2
            end

            #evaluate the fitness after the perturbation
            
            fitness_perturb_1 = fitness_function(new_genome_perturbed1,geo_info.pop_points)
            fitness_perturb_2 = fitness_function(new_genome_perturbed2,geo_info.pop_points)
            #calculate the drop in fitness, negative so that it can increases are good
            fit_diff_1 = fitness_perturb_1 - fitness_orig_1
            fit_diff_2 = fitness_perturb_2 - fitness_orig_2
            

            #create new individual structs with the data
            new_individual1 = Individual(fitness_orig_1,fit_diff_1,new_genome1)
            new_individual2 = Individual(fitness_orig_2,fit_diff_2,new_genome2)
            #add these individuals to the list of new children
            append!(new_children,[new_individual1,new_individual2]) 
        end
        
        # the assessement procedure
        # Add ALL children to population
        append!(population, new_children) # append is for extending lists with lists' individual elements
    
        #println("type of best ind $(typeof(best_ind))")

        elitisim
         if elitisim
             push!(population,Individual(solution_fitness,solution_fitness_diff,best_ind))
         end

         #println("population size $(length(population))")
         #println("population: $(population[1])")
        #println(population[1].fitness)
        #pop_type = [typeof(i) for i in population]
        sort!(population, by = p -> p.fitness, rev=true)
        
        #start with tournament selection
        #new_population = []
        #push!(new_population, population[1])
        #while length(new_population) < num_parents
        #    tournament = StatsBase.sample(rng, population, tournament_size, replace=false)
        #    tournament = sort!(tournament, by= indiv -> indiv.fitness, rev=true)
        #    winners = num_tournament_winners
        #    append!(new_population, tournament[1:winners]) # for some reason using num_tournament_winners doesn't work here
        #end
        new_population = selection_function(population,num_parents; kwargs...)
        population = new_population

        #take the individual which is in the top n% of fitness with the highest robustness
        #sort!(population, by = p -> p.fitness, rev=true)
       #n_th_percentile = length(population) ÷ (10)#find the index value corresponding to n% of the list
        #population_top_nth_percent = sort(population[1:n_th_percentile], by =  p->p.fitness_diff,rev = true)
        #solution = population_top_nth_percent[1]
        
        sort!(population, by = p -> p.fitness_diff, rev=true)
        solution = population[1]
        #@infiltrate(generation == total_generations-1)

        solution_fitness = solution.fitness
        solution_fitness_diff = solution.fitness_diff
        best_ind = solution.genome

       # if population[1].fitness > solution_fitness # if the new parent is the best found so far
       #     solution = population[1].genome                 # update best solution records
       #     solution_fitness = population[1].fitness
       #     solution_fitness_diff = population[1].fitness_diff
       #     best_ind = population[1].genome
       #     # solution_generation = generation_num
       # end

       #save all the fitness values for every individual in the population
       ind_fit_values = map((x) -> (x.fitness,x.fitness_diff),population)
       push!(population_stats,ind_fit_values)

       push!(best_ind_over_time,best_ind)
       fitness_over_time[generation] = (solution_fitness,solution_fitness_diff) # record the fitness of the current best over evolutionary time
       #solutions_over_time[generation] = solution
    end 
      return fitness_over_time,best_ind_over_time,population_stats
end

