#using Plots                                                                                                                                                                                                                                      
using Plots
using VoronoiCells
using GeometryBasics
using Random
using Shapefile
using GeoInterface
using NearestNeighbors
using StaticArrays
using DataFrames
using Statistics
using StatsBase
using Distributions
using Rasters   



abstract type AbstractIndividual end

struct Individual <: AbstractIndividual
    fitness::Float64
    fitness_diff::Float64
    genome::Vector{Point2{Float64}}
end

struct VoronoiIndividual <: AbstractIndividual
    fitness::Float64
    genome::Vector{Point2{Float64}}
    perimeters::Vector{Float64}
    areas::Vector{Float64}
    populations::Vector{Int64}

end



"""
------------------
Constructors
------------------
"""

"""
generate a test facility location
"""
function gen_fac_pos(my_border)
    is_in_bounds = false
    while is_in_bounds == false
        test_pos = Point2(rand(my_border.MBR.left:RAND_INTERVAL:my_border.MBR.right),rand(my_border.MBR.bottom:RAND_INTERVAL:my_border.MBR.top))
        #println("testpost: $test_pos")
        is_in_bounds = in_bounds(test_pos,my_border)
        if is_in_bounds
            return test_pos
        end 
    end
end

"""
determine if a candidate facility is within bounds
"""
function in_bounds(pos,my_border)
    return inpolygon(pos,my_border) ? true : false
end

"""
returns a single individual object with a randomized genome
"""
function generate_genome(my_border::Shapefile.Polygon,n_fac = 500)
    #fac_points = SVector{n_fac,Tuple{Float64,Float64}} #Vector{Tuple{Float64,Float64}}
    #build an array, check if each individual is in bounds
    fac_points = []
    
     while length(fac_points) <= n_fac
         fac_pos = gen_fac_pos(my_border)
             push!(fac_points,fac_pos)
         end
    static_fac_points = SVector{length(fac_points)}(fac_points)
    return static_fac_points
end




"""
-------------
Mutation Functions
-------------
"""


"""
given a dataframe of nearest neigbors and distances to each point, return the mean distance to every point
"""
function mean_distance(nn_df,beta = 1)
    grouped = groupby(nn_df,:idxs)
    fac_mean_dist = combine(grouped,:dist=> mean)#average the mean distance to the facility for each facility
    mean_fac_mean_dist = combine(fac_mean_dist,:dist_mean=>mean).dist_mean_mean[1]#average the mean facility desnity for each facility
    return mean_fac_mean_dist
end

function get_fitness(genome::SVector{Any},pop_points::Vector{Point2{Float64}},loss_function::Function = mean_distance)
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
    #@infiltrate
    kdtree = KDTree(genome_mat)
    idxs,dist = nn(kdtree,pop_points)#return the nearest facility and the distacne to the nearest facility
    return idxs,dist
end


function mutate(genome, boundary, num_inds_to_change = 1)
    """
    mutate: mutate the ganome of an individual, randomly relocating a number of facilities
    """
    #generate indices to mutate
    indices_to_change = sample(1:length(genome),num_inds_to_change)
    new_genome = deepcopy(genome)
    for index in indices_to_change #for each indicie sampled, move the facility to a random location
        #editing the array indices of an immutable struct..
        new_genome[index] = gen_fac_pos(boundary)
    end
    return new_genome 
end




"""
split_by(): takes in an array and a list of indicies, splits the array into sub arrays at the indicies
    input: 
        to_split<Vector{Any}>: a vector you want to split
        idx<Vector{Int64}>: a vector of indicies to split on
    returns: 
        <Vector{Vector{Any}}>: returns <to_split> split into sub arrays at indicies
"""
function split_by(to_split,idx)
    ranges = [(:)(i==1 ? 1 : idx[i-1]+1, idx[i]) for i in eachindex(idx)]
    push!(ranges,idx[end]+1:length(to_split))
    map(x->to_split[x],ranges) 
end

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
function crossover(genome1,genome2,crossover_points=2)
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


