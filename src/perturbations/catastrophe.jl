
"""
gen_catastrophe_pos
"""
function gen_catastrophe_pos(my_border, distribution)
    is_in_bounds = false
    while is_in_bounds == false
        y_uniform_distr = Uniform(my_border.MBR.bottom, my_border.MBR.top)
        y_uniform_vector = rand(y_uniform_distr, 10000) #1000 is an arbitrarily large number
        x_uniform_distr = Uniform(my_border.MBR.left, my_border.MBR.right)
        x_uniform_vector = rand(x_uniform_distr, 10000) #1000 is an arbitrarily large number
        sort!(x_uniform_vector)
        sort!(y_uniform_vector)
        # println(x_uniform_vector)
        weights = rand(distribution, 10000)
        sort!(weights)
        # print(size(weights))
        weights = AnalyticWeights(weights/maximum(broadcast(abs, weights)))
        # print(weights[2])
        x = wsample(x_uniform_vector, weights)
        sort!(weights, rev=true)
        y = wsample(y_uniform_vector,weights)
        # print(x)
        # print(y)
        # do wsample many times and create and plot a vector of Point2
        test_pos = Point2(x,y)
        #println("testpost: $test_pos")
        is_in_bounds = in_bounds(test_pos,my_border)
        if is_in_bounds
            # return x,y
            return test_pos
        end 
    end
end


function generate_genome_of_catastrophes(distribution, my_border,n_fac = 500)
    fac_points = [] #Vector{Tuple{Float64,Float64}}
    #build an array, check if each individual is in bounds
    x =[]
    y=[]
    distribution = (2.0,2.0)
    while length(fac_points) <= n_fac
        fac_pos = gen_catastrophe_pos(my_border, distribution)#gen_cata_point(my_border, distribution)#
            push!(fac_points,fac_pos)
    end
    # while length(x) <= n_fac
    #     xv, yv = gen_catastrophe_pos(my_border, distribution)
    #     append!(x, xv)
    #     append!(y, yv)
    # end
    scatter(x, y, label="data")
    savefig("scatter_cata")
    #static_fac_points = SVector{length(fac_points)}(fac_points)
    #return static_fac_points
    return fac_points
end


"""
simulate_catastrophe(): takes in a genome and removes the k-nearest facilities to a point selected from a distribution
input: 
    genome: genome to be modified
    size_of_catastrophe<Int64>: the catastrophe radius
    geo_info: the boundary of the zone we pick where the catastrophe will strike from, and pop_points
    distribution: distribution to use during gen_catastrophe_pos
output:
    perturbed_individual<Vector{Point2}>: the modified genome
"""
function perturb(gen, geo_info, size_of_catastrophe, perturbation_distribution) #my_border::Shapefile.Polygon, pop_points)
    genome = deepcopy(gen)

    #@infiltrate
    ##if size is -1, generate a remove a random facility
    if size_of_catastrophe == -1
        to_delete = rand(1:length(genome))
        genome = deleteat!(genome, [to_delete])
        return genome
    end

    genome_mat = transpose(mapreduce(permutedims, vcat, genome))#concat array of lon lat points into 2xn matrix
    kdtree = KDTree(genome_mat)
    point = gen_catastrophe_pos(geo_info.border, perturbation_distribution) #select a point within the border of the US from random distribution
    r = size_of_catastrophe
    idxs = inrange(kdtree, point, r, true) #returns all the facilities within a range r of point as an index vector
    count_removed = inrangecount(kdtree, point, r)
    
    genome = deleteat!(genome, idxs)

    if length(genome) > 0
        return genome
    else
        println("empty genome")
        return gen
    end
end


function percent_fac_retained(genome_old, genome_new)
    percent = length(genome_new)/length(genome_old) * 100
end

function robustness_retention_testing(best_ind_over_time_list,iterations, geo_info, perturbation_size, distribution)
    percents_list = []
    genome = csv_to_genome(best_ind_over_time_list)
    for i in range(1,1,iterations)
        #println(i)
        percent_retained = percent_fac_retained(genome, perturb(genome, geo_info, perturbation_size, distribution))
        push!(percents_list, percent_retained)
    end
    CSV.write("$iterations-percents_$best_ind_over_time_list.csv",  Tables.table(percents_list), writeheader=false)
end
