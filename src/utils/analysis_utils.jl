using DrWatson


"""
Misc. Timelines
"""

function has_repeated_vector(matrix)
    vector_sets = Set(matrix)
    return length(vector_sets) != size(matrix, 1)
end

"""
STATISTICAL SIGNIFICANE INTERVALS
"""

"""
clac_ci(): given a 1d slice of an array(a bunch of samples), reurn the <ci_percent> element of the sorted array 
input: 
    slice<Array{Any}>: a 1d array of samples for data
    ci_percent<Float64>: the confidence interval percentage to return 
output: 
    <Any>: the element at the <ci_percent> position in the sorted array

"""
calc_ci(slice,ci_percent) = sort(slice)[1+Int(floor(length(slice)*ci_percent))]#find the confidence interval value for a given array

"""
get_confidence_interval()
input:
    data_matrix<Matrix>: a MxN matrix where M is the number of runs of the fitness algorithim and M is the number of generations per run. 
    ci_low<Float64>: the lower bound percentage for the confidence interval
    ci_high<Float64>: the upper bound percentage for the confidence interval
output:
    ci_5<Array{Float64}: the 
"""
function get_confidence_interval(data_matrix,ci_low=.05,ci_high=.95)
    slices_array = [b for  (i,b) in enumerate(eachslice(data_matrix;dims=2))]#take an array and create a matrix of slices
    ci_low_arr = map(x->calc_ci(x,ci_low),slices_array)
    ci_high_arr = map(x->calc_ci(x,ci_high),slices_array)
    return ci_low_arr,ci_high_arr
end


"""
calc_bootstrap(): given a NxM matrix of evolutiaonary runs, where N is the number of runs and M is the number of generations in each one, calculate the mean value for each geneartion and the bootstrapped CI
input: 
 sample_matrix<NxM Matrix{Float64}>: a matrix of runs of evolutionary algorithim, N is the number of runs, M is the number of generations 
 ci_low<Float64>: the lower bound for the bootstrapped confidence interval
 ci_high<Float64>: the upper bound for the bootstrapped confidence interval
returns: 
 mean_ts: the mean value across all N runs for each generation 
 ci_low_ts: the lower bootstraped CI, specified by <ci_low>
 ci_high_ts: the upper bootstraped CI, specified by <ci_high>
"""
function calc_bootstrap(sample_matrix;ci_low_val = 0.05,ci_high_val = 0.95,n_bootstrap=1000)
    mean_ts = transpose(mean(sample_matrix,dims =1))
    #initialize bootstraping matrix
    bootstrap_matrix = zeros(n_bootstrap,size(sample_matrix)[2])
    #for each bootstrapping, sample all iterations of runs with replacement and calculate mean
    for i in 1:size(bootstrap_matrix,1)
        bootstrap_matrix[i,:] = mean(sample_matrix[sample(begin:end,num_runs,replace =true),:],dims =1)
    end
    #get confidence intercal from bootstrap
    ci_low_ts,ci_high_ts = get_confidence_interval(bootstrap_matrix,ci_low_val,ci_high_val)
    return mean_ts,ci_low_ts,ci_high_ts
end


function plot_genome(genome,geo_info;title = "")
    my_ind = make_voronoi_individual(genome,get_fitness,geo_info.border)
    idxs,dist = get_nearest_neighbors(genome,geo_info.pop_points)
    colors = [:red,:blue,:green,:purple,:yellow]
    #idx_color = Dict(i => i%length(colors) for i in collect(1:maximum(idxs)))
    idx_color = Dict(i => i%length(colors) for i in collect(1:maximum(idxs)))
    scatter(geo_info.pop_points,color =[idx_color[i] for i in idxs])
    plot!(my_ind.tess,color = :black)
    scatter!(Vector([i for i in genome]),color = :black,markersize = 5,label = "facilities",title = title)
    #annotate!([(genome[i][1] + 0.02, genome[i][2] + 0.03, Plots.text(i)) for i in eachindex(genome)])
    return current()
end



function plot_voronoi_raster(genome,fitness_function,geo_info;title = "")
    v_ind = genome
    #v_ind = make_voronoi_individual(genome,fitness_function,geo_info)
    mbr_rect = convertMBRtoRectangle(geo_info.MBR)
    my_tess = voronoicells(Vector([i for i in genome]),mbr_rect)
    tess_poly = GeometryBasics.Polygon.(my_tess.Cells)
    tess_mp = MultiPolygon(tess_poly)
    ras = geo_info.population.pop_points
    ras_tess = create_rasterized_tess(tess_mp,ras)
    ras_tess = ras_tess .|> i -> replace_missing(i,0.0)
    ras_tess_collected = collect_voronoi_raster(ras_tess)
    heatmap(ras_tess_collected)
    return current()
end

function plot_population_with_tess(genome,fitness_function,geo_info;title = "")
    #v_ind = make_voronoi_individual(genome,fitness_function,geo_info)
    v_ind  = genome
    mbr_rect = convertMBRtoRectangle(geo_info.MBR)
    my_tess = voronoicells(Vector([i for i in genome]),mbr_rect)
    tess_poly = GeometryBasics.Polygon.(my_tess.Cells)
    tess_mp = MultiPolygon(tess_poly)
    heatmap(geo_info.population.pop_points,legend = false)
    plot!(my_tess,color = :white,legend = false)
    scatter!(genome,color = :white,title = title,legend = false)
    return current()
end





