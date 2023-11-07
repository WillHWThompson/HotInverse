using DrWatson

"""
SELECTION METHODS
"""

"""
tournament_selection()
takes in a list of <Individuals> structs and returns the tournament selection
"""

function tournament_selection(my_population,num_parents;tournament_size=30,num_tournament_winners=1,kwargs...)
#start with tournament selection
    new_my_population = []
    push!(new_my_population, my_population[1])
    while length(new_my_population) < num_parents
        tournament = StatsBase.sample(rng, my_population, tournament_size, replace=false)
        tournament = sort!(tournament, by= indiv -> indiv.fitness, rev=true)
        winners = num_tournament_winners
        append!(new_my_population, tournament[1:winners]) # for some reason using num_tournament_winners doesn't work here
    end
    return new_my_population
end


"""
NON-DOMINATED SORTING

""" 
function pareto_domination(individual_1::Individual, individual_2::Individual;objectives::Vector= [:fitness,:fitness_diff]) 
    a = false
    for objective in objectives
        individual_1_objective = getfield(individual_1,objective)
        individual_2_objective = getfield(individual_2,objective)
        if individual_1_objective > individual_2_objective
            a = true
        elseif individual_2_objective > individual_1_objective
            return false
        end
    end
    return a
end
    

function pareto_nondominated_front(population;objectives = [:fitness,:fitness_diff],kwargs...)
    front = [Individual(-9999,-9999,[Point2(0.,0.)])] # initialize an empty individual]
    for individual in population
        front = union(front, [individual])
        for front_individual in setdiff(front, [individual]) ##IS G gonna be removed?
            if pareto_domination(front_individual, individual,objectives = objectives)
                front = setdiff(front, [individual])
                break
            elseif pareto_domination(individual, front_individual)
                front = setdiff(front, [front_individual])
            end
        end
    end
    return front
end

function rank_assignment(population,num_parents;objectives = [:fitness,:fitness_diff],dominated_by = 2,kwargs...)
    population_buffer = deepcopy(population)
    population_final = []
    rank = 1
    fronts = []
    while length(population_buffer) > 0
        front = pareto_nondominated_front(population_buffer;kwargs...)
        sort!(front, by = p -> getfield(p,objectives[2]),rev=true)
        front = shuffle(front)#shuffle the front instead of selecting 
        #sort!(front, by = p -> getfield(p,objectives[2]),rev=true)
        push!(fronts,front)#add each front to the list
        population_buffer = setdiff(population_buffer, front)
        rank+=1

    end
    fronts_flat = reduce(vcat,fronts)#flatten to vector of individuals
    fronts_truncated = fronts_flat[1:num_parents]
    return fronts_truncated
end
