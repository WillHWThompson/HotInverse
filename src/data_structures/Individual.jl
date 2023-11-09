
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
function gen_fac_pos(my_border,mbr,check_bounds = false)
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


function generate_genome(my_geo_info::GeoInfo;n_fac = 500,kwargs...)
    #fac_points = SVector{n_fac,Tuple{Float64,Float64}} #Vector{Tuple{Float64,Float64}}
    #build an array, check if each individual is in bounds
    @show n_fac
    fac_points = []
     while length(fac_points) <= n_fac
         fac_pos = gen_fac_pos(my_geo_info.border,my_geo_info.MBR,kwargs...)
             push!(fac_points,fac_pos)
         end
    static_fac_points = SVector{length(fac_points)}(fac_points)
    return static_fac_points
end





