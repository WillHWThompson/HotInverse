#imports
using StatsBase
import Base.Threads.@threads
const RAND_INTERVAL = 0.001

#setup
rng = Random.MersenneTwister() #unseeded



function sim_anneal(fitness_function::Function,generate_genome_function::Function,constraint_function::Function,mutation_function::Function;
                                num_elements_to_mutate=1,total_generations=100,genome_length=500,
                                num_parents=10,num_perturbations = 0,perturbation_size=1,perturbation_distribution = Pareto(),
                                init_t = 1000,delta_t = 1,kwargs... )

    """simulated annealing 

    parameters: 
    fitness_function: function that returns the fitness of a genome given the genome as an input parameter
    gen_individual_function: function returns a new individual in the population
    geo_info<struct> contains geo_info.border which specifies the border in which to place points and 
                        geo_info.population which specifies the population as points
    total_generations: (int) number of total iterations for stopping67 condition
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
    solution_fitness = -99999 # fitness of best genome so far
    #solution_generation = 0 # time (generations) when solution was found
    #solutions_over_time = zeros((total_generations, genome_length))
    solutions_over_time = []
    #diversity_over_time = zeros(total_generations)
    
    fitness_over_time = []
    best_ind_over_time = [] 
    best_ind = nothing 

    #create a initial random solution
    #genome = generate_genome_function()
    #ind = make_voronoi_individual(genome,fitness_function,geo_info.border)
    ind = generate_ind_function()


    solution = ind
    solution_fitness = ind.fitness
    solution_fitness_diff = -999

    best_ind = ind

    t = init_t
    generation = 0
    
    population_stats = [] 

    #set delta t to number after the decimal so 99 => 0.99
    delta_t = delta_t/(10^(ceil(log10(delta_t))))

    #a/(10^(ceil(log10(a)))
    println("Starting Annealing Process...")
    for t in 1:total_generations
        if t % 50 == 0
           println("\tgeneration $t")
        end
        new_individual1 = mutation_function(ind)

        if constraint_function(new_individual1)
            fit_diff = new_individual1.fitness-ind.fitness#the difference in fitness between this generation and the next
            quality_exp = exp((fit_diff)/(init_t*(delta_t^t)))

            if fit_diff > 0#if the fitness of the new individual is higher
                #println("accepting new ind")
                ind = new_individual1 elseif rand() < quality_exp#or with a certian probabilty set by the schedule println("accepting new ind with lower fitness")
                ind = new_individual1
            end
        else
           a= constraint_function(new_individual1)
           #println("individual violated constraint. value: $value, constraint: $my_constraint")
        end
        solution = ind
        solution_fitness = solution.fitness
        #@show solution_fitness

        #best_ind = solution.genome
        best_ind = solution

        push!(best_ind_over_time,best_ind)
        push!(fitness_over_time,(solution_fitness)) # record the fitness of the current best over evolutionary time
      
        ind_fit_values = ind.fitness
        push!(population_stats,ind_fit_values)
        #@infiltrate
        #solutions_over_time[generation] = solution
        push!(solutions_over_time,solution)

    end
     

    return fitness_over_time,best_ind_over_time,population_stats
end










function sim_anneal(fitness_function::Function,generate_genome_function::Function,geo_info;
                                num_elements_to_mutate=1,total_generations=100,genome_length=500,
                                num_parents=10,
                                init_t = 1000,delta_t = 1,kwargs... )

    """simulated annealing 

    parameters: 
    fitness_function: function that returns the fitness of a genome given the genome as an input parameter
    gen_individual_function: function returns a new individual in the population
    geo_info<struct> contains geo_info.border which specifies the border in which to place points and 
                        geo_info.population which specifies the population as points
    total_generations: (int) number of total iterations for stopping67 condition
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
    solution_fitness = -99999 # fitness of best genome so far
    #solution_generation = 0 # time (generations) when solution was found
    solutions_over_time = zeros((total_generations, genome_length),dtype=int)
    #diversity_over_time = zeros(total_generations)
    
    fitness_over_time = []
    best_ind_over_time = [] 
    best_ind = nothing 


    #create a initial random solution
    new_genome = generate_genome_function(geo_info.border,genome_length)
    new_fitness = fitness_function(new_genome,geo_info.pop_points)
    ind = Individual(new_fitness,-999,new_genome)

    solution = ind
    solution_fitness = new_fitness
    solution_fitness_diff = -999

    best_ind = ind

    t = init_t
    generation = 0
    
    population_stats = [] 

    #set delta t to number after the decimal so 99 => 0.99
    delta_t = delta_t/(10^(ceil(log10(delta_t))))

    #a/(10^(ceil(log10(a)))
    for t in 1:total_generations
        #println("generation $t")
        
        new_genome1 = mutate(ind.genome,geo_info.border,num_elements_to_mutate)
        fitness_orig_1 = fitness_function(new_genome1,geo_info.pop_points)

        new_genome_perturbed1 = new_genome1

            if num_perturbations > 0
            # do perturbation
                for i in 1:num_perturbations
                    new_genome_perturbed1 = perturb(new_genome_perturbed1, geo_info, perturbation_size, perturbation_distribution)
                end
            else
                new_genome_perturbed1 = new_genome1
            end

            fitness_perturb_1 = fitness_function(new_genome_perturbed1,geo_info.pop_points)
            fit_diff_1 = fitness_perturb_1 - fitness_orig_1


            #create new individual structs with the data
            #new_individual1 = Individual(fitness_function(new_genome_perturbed1,-999,geo_info.pop_points),new_genome1)
            new_individual1 = Individual(fitness_orig_1,fit_diff_1,new_genome1)

            #decrease the temperature 
            #calcualte the difference in fitness

            #println("ind fitness: $(ind.fitness), new ind fitness: $(new_individual1.fitness)") 

            fit_diff = new_individual1.fitness-ind.fitness#the difference in fitness between this generation and the next
            quality_exp = exp((fit_diff)/(init_t*(delta_t^t)))

            if constraint_function(new_individual1)
                if fit_diff > 0#if the fitness of the new individual is higher
                    #println("accepting new ind")
                    ind = new_individual1
                elseif rand() < quality_exp#or with a certian probabilty set by the schedule
                    #println("accepting new ind with lower fitness")
                    ind = new_individual1
                end
            end


        solution = ind
        solution_fitness = solution.fitness
        solution_fitness_diff = solution.fitness_diff
        #@show solution_fitness_diff
        best_ind = solution.genome

        push!(best_ind_over_time,best_ind)
        push!(fitness_over_time,(solution_fitness,solution_fitness_diff)) # record the fitness of the current best over evolutionary time
       
        ind_fit_values = [(ind.fitness,ind.fitness_diff)]
        push!(population_stats,ind_fit_values)
        #@infiltrate
        #solutions_over_time[generation] = solution
    end
     

    return fitness_over_time,best_ind_over_time,population_stats
end
#end
