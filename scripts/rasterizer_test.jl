
using Revise
Revise.revise()
using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise
include(srcdir("HotInverse.jl"))#this line will make all the code available
#using LibGEOS
using LinearAlgebra


using GeometryBasics
using Rasters

abstract type Population end

struct PopulationRaster{P<:AbstractRaster} <: Population
    pop_points::P
end

struct PopulationPoints{P<:AbstractArray{Point2{Float64}}} <: Population
    pop_points::P
end


struct MBR{P}
    mbr::P
    left::Real
    right::Real
    top::Real
    bottom::Real
end


struct GeoInfo{ T <: AbstractGeometry}
    border::T
    MBR::MBR
    population::Population
end


"""
    mbr(polygon::AbstractGeometry)

Compute the minimum bounding rectangle (MBR) of a polygon.

# Arguments
- `polygon::AbstractGeometry`: A polygon represented as an `AbstractGeometry` object.

# Returns
- `MBR`: A minimum bounding rectangle object.

# Example
"""
function mbr(polygon::AbstractGeometry)
    my_mbr_rect = Rect(coordinates(polygon.exterior))
    mbr_extrema = extrema(coordinates(my_mbr_rect))
    return MBR(my_mbr_rect,mbr_extrema[1][1],mbr_extrema[2][1],mbr_extrema[1][2],mbr_extrema[2][2])
end

function mbr(shp_mbr::Shapefile.Rect)
    mbr_poly = convertShapefileMBRtoPoly(shp_mbr)
    return mbr(mbr_poly)
end



"""
population
"""
# function generate_pop_points(dist;min::Float64 = 0.,::Float64 = 1.0,inc::Float64 = 0.1)
#     x = min:inc:max
#     y = min:inc:max
#     #pop tuples
#     pop_tuples = [[i,j] for i in x, j in y]
#     pop_points = Point2.(pop_tuples)
#     pdf_values = map( x -> pdf(dist,x),pop_tuples)
#     return vec(pop_points),pdf_values
# end



# function generate_pop_points(dist;min::Float64 = 0.,max::Float64 = 1.0,inc::Float64 = 0.1)
#     x = min:inc:max
#     y = min:inc:max
#     #pop tuples
#     pop_tuples = [[i,j] for i in x, j in y]
#     pop_points = Point2.(pop_tuples)
#     pdf_values = map( x -> pdf(dist,x),pop_tuples)
#     return vec(pop_points),pdf_values
# end


# function generate_pop_points_mvGaussian(;μ = [0.5,0.5], Σ = [1.0 0.0;0.0 1.0],kwargs...)
#     dist  = MvNormal(μ,Σ)
#     return Population(generate_pop_points(dist;kwargs...)...)

# end


"""
Geometry Utils
"""
function convertRecttoPoly(rect::Rect)
    coordinates(rect) |> collect .|> Point2 |> Polygon
end


"""
Raster Generation 
"""

"""
    generate_pdf_raster(dist; my_width = 1, my_height = 1, my_res = 100)

Generate a raster of a probability density function `dist` over a rectangular region.

# Arguments
- `dist`: A probability density function.
- `my_width`: The width of the rectangular region.
- `my_height`: The height of the rectangular region.
- `my_res`: The resolution of the raster.

# Returns
A `Raster` object representing the probability density function over the rectangular region.
"""
function generate_pdf_raster(dist; my_width = 1, my_height = 1, my_res = 100)
    x_inc = my_width/my_res
    y_inc = my_height/my_res

    x_range =0:x_inc:my_width
    y_range= 0:y_inc:my_height

    raster = Raster(zeros(X(x_range),Y(y_range)))

    for (i,j) in Iterators.product(1:length(x_range),1:length(y_range))
        raster[i,j] = pdf(dist,[x_range[i],y_range[j]])
    end
    return raster
end


"""
    init_geoinfo(width::Real,height::Real,res::Real,distribution::Distribution)

Initialize the geometry information by generating a rectangle with the given width and height, 
a minimum bounding rectangle (mbr), a probability density function (pdf) raster, and a population raster. 

# Arguments
- `width::Real`: The width of the rectangle.
- `height::Real`: The height of the rectangle.
- `res::Real`: The resolution of the pdf raster.
- `distribution::Distribution`: The probability distribution function.

# Returns
- `GeoInfo`: A structure containing the geometry information.

# Example
"""
function init_geoinfo(width::Real,height::Real,res::Real,distribution::Distribution)
    #generate geomtery basics rectangle with widht and height
    my_border = Rect([Point2(0,0),Point2(width,height)]) |> convertRecttoPoly
    my_mbr = mbr(my_border)
    my_raster = generate_pdf_raster(distribution,my_width = width, my_height = height, my_res = res)
    println("raster generated")
    my_pop = PopulationRaster(my_raster)
    return GeoInfo(my_border,my_mbr,my_pop)    
    #return my_border,my_mbr,my_pop
end



dist = MvNormal([0.5,0.5],[1.0 0.0;0.0 1.0])
init_geoinfo(1,1,100,dist)

us_data_dir = srcdir("data","rect_data")

@which init_geoinfo(us_data_dir,pop_dir = "/pop_points/pop_points.csv")
init_geoinfo(us_data_dir,pop_dir = "/pop_points/pop_points.csv")


a.border |> convertShapefilePolytoPoly


genome = generate_genome(geo_info.border,10)
polygon = convertShapefileMBRtoRectangle(geo_info.border.MBR)
mbr = generate_MBR_from_points(geo_info.pop_points,delta = 5.0) 


#generate polygon
using GeometryBasics

# Define the vertices of the polygon
my_points = [Point2(0, 0), Point2(0, 1), Point2(1, 1), Point2(1, 0)]

x = X(0:0.1:10)
y = Y(0:0.1:10)
raster = Raster(rand(y,x))

# Create a polygon from the points
parts = [0]
polygon = Polygon(my_points, parts)

# Print the polygon

Polygon(Poinpoints)

plot(polygon)

rasterize(polygon,fill = 1,to =  raster)

# using Downloads
# # Download a borders shapefile
# shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
# shapefile_name = "country_borders.shp"
# isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# # Load the shapes for china
# china_border = Shapefile.Handle(shapefile_name).shapes[10]

# # Rasterize the border polygon
# china = rasterize(last, china_border; res=0.1, missingval=0, fill=1, boundary=:touches, progress=false)


# shp =  open(shapefile_name) |> read



using Shapefile, GeometryBasics

# Load the shapefile

# Get the first polygon in the shapefile
#poly = shp[1].shape

# Define a point to check
#pt = Point2(0.5, 0.5)

# Check if the point is inside the polygon
using GeometryTypes
poly = china_border

#poly = Polygon(Point{2, Int}[(3, 1), (4, 4), (2, 4), (1, 2), (3, 1)])

poly = geo_info.border

gb_poly = Polygon([Point2(point.x, point.y) for point in poly.points[2:end]])

generate_MBR_from_Poly(gb_poly)


poly.points




is_inside = GeometryTypes.in(pt, gb_poly)

poly = Circle(Point((0.,0.)), 1.)

rect(gb_poly)

generate_MBR_from_points(gb_poly.points,delta = 5.0)

gb_poly.points


using GeometryBasics

# Define a polygon
poly = Polygon([Point2(0, 0), Point2(0, 1), Point2(1, 1), Point2(1, 0)])


GeometryTypes.Polygon([Point2(0, 0), Point2(0, 1), Point2(1, 1), Point2(1, 0)])
# Get the bounding box of the polygon



abstract type GeoInfo end

struct GeoInfoShapefile <: GeoInfo
    border::Shapefile.Polygon
    pop_points::Vector{Point2{Float64}}
end




# Print the bounding box
println("Bounding box: $bbox")

function generate_genome(my_border::GeometryBasics.Polygon,n_fac = 500)
    #fac_points = SVector{n_fac,Tuple{Float64,Float64}} #Vector{Tuple{Float64,Float64}}
    #build an array, check if each individual is in bounds
    fac_points = []
    mbr = generate_MBR_from_Poly(gb_poly)
     while length(fac_points) <= n_fac
         fac_pos = gen_fac_pos(my_border,mbr)
             push!(fac_points,fac_pos)
         end
    static_fac_points = SVector{length(fac_points)}(fac_points)
    return static_fac_points
end






extrema(gb_poly)

generate_genome(gb_poly)


extrema


# Print the result
