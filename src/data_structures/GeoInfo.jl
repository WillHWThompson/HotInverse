
using Shapefile

struct GeographicInfo
    border::Shapefile.Polygon
    pop_points::Vector{Point2{Float64}}
end

struct GeographicInfo
    border::Shapefile.Polygon
    pop_points::Vector{Point2{Float64}}
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

    """
    GET BOUNDARY SHAPE FILE
    """
     BOUNDARY_LOC = "$dir/$boundary_dir"
     us_shp = Shapefile.Handle(BOUNDARY_LOC)
     us_border = Shapefile.shapes(us_shp)[1]
     geo_info = GeographicInfo(us_border,pop_points)
    
     return geo_info
end


function convertShapefileMBRtoRectangle(ShapefileMBR)
"""
convertShapefileMBRtoRectangle
input: a Shapefile Minimum Bounding Ractangle to a GeometryBasics Rectangle
output: A GeometryBasics Rectangle of the same shape
"""

    return Rectangle(Point2(ShapefileMBR.left,ShapefileMBR.bottom),Point2(ShapefileMBR.right,ShapefileMBR.top))
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
