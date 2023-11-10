
"""
-------------
Mutation Functions
-------------
"""

function mean_distance(nn_df,beta = 1)
"""
given a dataframe of nearest neigbors and distances to each point, return the mean distance to every point
"""
    return -mean(nn_df[!,:dist])
end


function mean_distance_by_cell(nn_df,beta = 1)
"""
given a dataframe of nearest neigbors and distances to each point, return the mean distance to every point
"""
    grouped = groupby(nn_df,:idxs)
    fac_mean_dist = combine(grouped,:dist=> mean)#average the mean distance to the facility for each facility
    mean_fac_mean_dist = combine(fac_mean_dist,:dist_mean=>mean).dist_mean_mean[1]#average the mean facility desnity for each facility
    return mean_fac_mean_dist
end


function get_fitness(genome,pop_points::Vector{Point2{Float64}},loss_function::Function = mean_distance)
"""
get_fitness:given an individual and a loss function, compute the fitness of the individual
inputs: 
    genome::SVector - A static vector containing the locations of the facilities 
    pop_points::Vector{Point2{Float64}} - A vector containing the locations of citizens
    loss_function::Function - a the function to calculate the fitness with - takes an a set of arrays and takes out the value of the fitness
output:
    return fitness::Float64 - the fitness value for the individual
"""
    idxs,dist = get_nearest_neighbors(genome,pop_points)#perform neartest neighbors search
    #convert result to dataframe
    df = DataFrame([idxs,dist],["idxs","dist"])
    fitness = loss_function(df)
    return fitness
end


function get_fitness(idxs::Vector{Int64},dist::Vector{Float64};loss_function = mean_distance)
"""
get_fitness:given the results of a nearest neighbors search, return the fitness of an individual
inputs: 
    idxs::Vector{Int64} - the index of the nearest neighbor facility for each inividual in the population
    dist::Vector{Float64} - the distance to the nearest neighbor facility for each inividual in the population
output:
    return fitness::Float64 - the fitness value for the individual
"""
    df = DataFrame([idxs,dist],["idxs","dist"])
    fitness = loss_function(df)
    return fitness
end

function get_nearest_neighbors(genome,pop_points)
   """
   Given a genome, a set of facilites and a set of population points, determine the nearest facility to each pop point
   inputs: 
    genome::Genome - a genome for a facility layout
    pop_points:: a list of Point2 locations for each pop
   """
    genome_mat = transpose(mapreduce(permutedims, vcat, genome))#concat array of lon lat points into 2xn matrix
    kdtree = KDTree(genome_mat)
    @infiltrate
    idxs,dist = nn(kdtree,pop_points)#return the nearest facility and the distacne to the nearest facility
    return idxs,dist
end


function mutate(genome, boundary; num_inds_to_change = 1,add_fac_prob = 0.0,remove_fac_prob = 0.0)
    """
    mutate: mutate the ganome of an individual, randomly relocating a number of facilities
    """
    #generate indices to mutate
    indices_to_change = sample(1:length(genome),num_inds_to_change)

    new_genome = Vector(deepcopy(genome))
    for index in indices_to_change #for each indicie sampled, move the facility to a random location
        #editing the array indices of an immutable struct..
        new_genome[index] = gen_fac_pos(boundary)
    end

    if rand()<add_fac_prob
        push!(new_genome, gen_fac_pos(boundary))
    end
    if (rand()< remove_fac_prob) & (length(new_genome) > 1)
        index_to_remove = sample(1:length(genome))
        splice!(new_genome,index_to_remove)
    end





    return SVector{(length(new_genome)),Point2{Float64}}(new_genome)
end



#function mutate(genome, boundary, num_inds_to_change = 1)
#    """
#    mutate: mutate the ganome of an individual, randomly relocating a number of facilities
#    """
#    #generate indices to mutate
#    indices_to_change = sample(1:length(genome),num_inds_to_change)
#    new_genome = Vector(deepcopy(genome))
#    for index in indices_to_change #for each indicie sampled, move the facility to a random location
#        #editing the array indices of an immutable struct..
#        new_genome[index] = gen_fac_pos(boundary)
#    end
#
#    return SVector{(length(new_genome)),Point2{Float64}}(new_genome)
#end






function split_by(to_split,idx)
"""
split_by(): takes in an array and a list of indicies, splits the array into sub arrays at the indicies
    input: 
        to_split<Vector{Any}>: a vector you want to split
        idx<Vector{Int64}>: a vector of indicies to split on
    returns: 
        <Vector{Vector{Any}}>: returns <to_split> split into sub arrays at indicies
"""
    ranges = [(:)(i==1 ? 1 : idx[i-1]+1, idx[i]) for i in eachindex(idx)]
    push!(ranges,idx[end]+1:length(to_split))
    map(x->to_split[x],ranges) 
end

function crossover(genome1,genome2,crossover_points=2)
"""
crossover(): take two genomes and perform an n_point crossover
input: 
    genome1<Vector{Point2}>: parent 1's genome
    genome2<Vector{Point2}>: parent 2's genome
    crossover_points<Int64>: the number of points at which to perform the crossover
output:
    genome1<Vector{Point2}>: the first crossover over genome
    genome2<Vector{Point2}>: the second crossover over genome
"""
    new_genome1 = []
    new_genome2 = []
    #generate random indicies, crossover will happen at these points
    gen_length = length(genome1)
    inds = sort(sample(1:gen_length,crossover_points))
    #split the genomes at the crossover points
    split1 = split_by(genome1,inds)
    split2 = split_by(genome2,inds)
    #swap out alternating parts of the genome
    for i in eachindex(split1)
        if i%2 == 0
            append!(new_genome1,split1[i])
            append!(new_genome2,split2[i])
        else
            append!(new_genome1,split2[i])
            append!(new_genome2,split1[i])
        end
    #println("length of new_genome1: $(length(new_genome1))")
    #println("legth split[$i]: $(length(split[i]))")
    end
    return new_genome1,new_genome2
end

