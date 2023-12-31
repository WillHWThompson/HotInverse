
abstract type AbstractIndividual end

struct Individual <: AbstractIndividual
    fitness::Float64
    fitness_diff::Float64
    genome::Vector{Point2{Float64}}
end

"""
------------------
Constructors
------------------
"""

"""
generate a test facility location
"""
function gen_fac_pos(pop_data::PopulationPoints,geo_info::GeoInfo;check_bounds::Bool = false)

    my_border = geo_info.border
    mbr = geo_info.MBR

    is_in_bounds = false
    while is_in_bounds == false
        test_pos = Point2(rand(mbr.left:RAND_INTERVAL:mbr.right),rand(mbr.bottom:RAND_INTERVAL:mbr.top))
        is_in_bounds =  check_bounds ? in_bounds(test_pos,my_border) : true #if check bounds is true, ensure point falls within polygon
        if is_in_bounds
            return test_pos
        end 
    end
end



"""
For Rasters
"""
function gen_fac_pos(pop_data::PopulationRaster,geo_info::GeoInfo;n_facs::Int = 10)
    ras = pop_data.pop_points
    my_size = size(ras)#get the size of the Raster in x and y
    range_tuples = map((i,j) -> (i,j),(1,1),my_size)#create ranes 1_size_x etc. 
    ranges = map(x-> UnitRange(x...),range_tuples)#convert to Range objects


    in_bounds_index_pairs = []
    #@generate facility locatons and ensure theyare in bounds, iterate until we have n_facs
    while length(in_bounds_index_pairs) < n_facs
        point_indicies = map(x-> sample(x,(1,n_facs)),ranges)#sample random points from these range objects for x and y
        indices_pairs = map((i,j) -> [i,j],point_indicies...)#get pairs of indicie values
        is_in_bounds = map(index_i -> ras[index_i...] != MISSING_VAL,indices_pairs)#is the proposed facility in bounds?
        append!(in_bounds_index_pairs, indices_pairs[is_in_bounds])
    end
    in_bounds_index_pairs  =   in_bounds_index_pairs[1:n_facs]#cut to exactly n facs

    in_bound_hcat = hcat(in_bounds_index_pairs...)#create a 2d matrix
    in_bound_seperated = [in_bound_hcat[i,:] for i in 1:size(in_bound_hcat,1)]#and break it into X and Y arrays
    point_coords = map(i -> ras.dims[i][vec(in_bound_seperated[i])],1:length(ras.dims))#get the coord values form raster index
    point_pairs = map((i,j) -> [i,j],point_coords...)
    fac_points = GeometryBasics.Point2.(point_pairs)#conver to Point2 structs
    return fac_points
end




"""
determine if a candidate facility is within bounds
"""
function in_bounds(pos,my_border::Shapefile.Polygon)
    return inpolygon(pos,my_border) ? true : false
end

function in_bounds(pos,my_border::GeometryBasics.Polygon)
    return in(pos,my_border) ? true : false
end


"""
returns a single individual object with a randomized genome
"""
function generate_genome(my_border::Shapefile.Polygon;n_fac = 500,kwargs...)
    #fac_points = SVector{n_fac,Tuple{Float64,Float64}} #Vector{Tuple{Float64,Float64}}
    #build an array, check if each individual is in bounds
    fac_points = []
    
     while length(fac_points) <= n_fac
         fac_pos = gen_fac_pos(my_border,check_bounds = check_bounds)
             push!(fac_points,fac_pos)
         end
    static_fac_points = SVector{length(fac_points)}(fac_points)
    return static_fac_points

end


function generate_genome(my_population::PopulationPoints,my_geo_info::GeoInfo;n_facs = 10,kwargs...)
    #fac_points = SVector{n_fac,Tuple{Float64,Float64}} #Vector{Tuple{Float64,Float64}}
    #build an array, check if each individual is in bounds
    fac_points = []
     while length(fac_points) <= n_facs
         fac_pos = gen_fac_pos(my_geo_info.border,my_geo_info.MBR,kwargs...)
             push!(fac_points,fac_pos)
         end
    static_fac_points = SVector{length(fac_points)}(fac_points)
    return static_fac_points
end

function generate_genome(my_population::PopulationRaster,my_geo_info::GeoInfo;n_facs = 10,kwargs...)
    fac_points = gen_fac_pos(my_population,geo_info)
    static_fac_points = SVector{length(fac_points)}(fac_points)
    return static_fac_points

end




