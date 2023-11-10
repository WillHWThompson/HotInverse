
using Shapefile

const MISSING_VAL = -1


#abstract type AbstractGeoInfo end
#
#struct GeographicInfo
#    border::Shapefile.Polygon
#    pop_points::Vector{Point2{Float64}}
#end
#
#struct GeographicInfo
#    border::Shapefile.Polygon
#    pop_points::Vector{Point2{Float64}}
#end



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


struct GeoInfo{T}
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
    return MBR(my_mbr_rect,mbr_extrema[1][1],mbr_extrema[2][1],mbr_extrema[2][2],mbr_extrema[1][2])
end

function mbr(shp_mbr::Shapefile.Rect)
    mbr_poly = convertShapefileMBRtoPoly(shp_mbr)
    return mbr(mbr_poly)
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

    raster = Raster(zeros(X(x_range;sampling = Rasters.Intervals()),Y(y_range;sampling = Rasters.Intervals())),missingval = MISSING_VAL)

    for (i,j) in Iterators.product(1:length(x_range),1:length(y_range))
        raster[i,j] = pdf(dist,[x_range[i],y_range[j]])
    end
    return raster
end

"""
Constructors
"""

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

function init_geoinfo(dir::String,;pop_dir::String = "pop_points/pop_points.csv",boundary_dir::String = "boundary/boundary.shp")
"""
get_geo_info(): returns a GeoGraphicInfo struct with all the points in the population as well as the shape of the boundary
inputs: 
   dir::String - the location of the parent dir with both the pop_points and boundarys as shape files
   pop_dir::String: the relative path from <dir> to the csv file with population in it 
   boundary_dir::String - relative path to the .shp file containing all the boundary info NOTE: same directory must also contain all ancilary files for a .shp file 
returns: 
    geo_info::GeoGraphicInfo -  a datastructure that contains both the border and pop_points
"""

    """
    LOAD POPULATION DATA
    """
    println("loading population data")
    pop_path = "$dir/$pop_dir"
    pop_points = load_pop_data(pop_path)
    pop_struct = PopulationPoints(pop_points)

    """
    GET BOUNDARY SHAPE FILE
    """
     BOUNDARY_LOC = "$dir/$boundary_dir"
     us_shp = Shapefile.Handle(BOUNDARY_LOC)
     us_border = Shapefile.shapes(us_shp)[1]
    us_border_poly = convertShapefilePolytoPoly(us_border)
    border_mbr = mbr(us_border_poly)

     geo_info = GeoInfo(us_border_poly,border_mbr,pop_struct)

     return geo_info

    
end


function convertShapefileMBRtoRectangle(ShapefileMBR)
"""
convertShapefileMBRtoRectangle
input: a Shapefile Minimum Bounding Ractangle to a GeometryBasics Rectangle
output: A GeometryBasics Rectangle of the same shape
"""

    return Polygon(Point2(ShapefileMBR.left,ShapefileMBR.bottom),Point2(ShapefileMBR.right,ShapefileMBR.top))
end


function convertShapefileMBRtoPoly(ShapefileMBR)
"""
convertShapefileMBRtoRectangle
input: a Shapefile Minimum Bounding Ractangle to a GeometryBasics Rectangle
output: A GeometryBasics Rectangle of the same shape
"""
    Polygon([Point2(ShapefileMBR.left,ShapefileMBR.bottom),Point2(ShapefileMBR.left,ShapefileMBR.top),Point2(ShapefileMBR.right,ShapefileMBR.top),Point2(ShapefileMBR.right,ShapefileMBR.bottom)])
end


function convertShapefilePolytoPoly(ShapefilePoly)
    """
    Given a Shapefile.Poly object convert to a GeometryBasics.Polygon
    inputs: 
        ShapefilePoly::Shapefile.Poly - a polygon loaded from a shapefile
    outputs: 
        polygon::GeometryBasics.Polygon - the converted polygon
    """
    polygon = GeometryBasics.Polygon(map(x -> Point2(x.x,x.y),ShapefilePoly.points))
    print("hello")
    return polygon
end

function generate_MBR_from_points(my_points;delta::Float64= 0.01)
    """
    Given a set a points, constuct a minimum bounding rectangle for Vornoi Tesselation
    inputs: 
        my_points::Vector{Any}, must contain a vector a interables where the iterable has two elements for x and y
        delta::Float64, the padding added to the minimum bounding rectangle
    outputs: 
        my_mbr: GeometryBasics.Rectangle - the returned mbr
    """

    x_points = map(x->x[1], my_points)
    y_points = map(x->x[2], my_points)

    x_max = maximum(x_points)+delta
    x_min = minimum(x_points)-delta

    y_max = maximum(y_points)+delta
    y_min = minimum(y_points)-delta

    my_mbr = Rectangle(Point2(x_min,y_min),Point2(x_max,y_max))
    return my_mbr
end
