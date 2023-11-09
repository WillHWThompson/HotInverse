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
## %
us_data_dir = srcdir("data","rect_data")
dist = MvNormal([0.5,0.5],[1.0 0.0;0.0 1.0])
raster_geo_info = init_geoinfo(1,1,100,dist)
points_geo_info = init_geoinfo(us_data_dir,pop_dir = "/pop_points/pop_points.csv")
geo_info = raster_geo_info
##%
genome = generate_genome(geo_info,n_fac = 10)
#genome = generate_genome(geo_info)
#make_voronoi_individual(genome,get_fitness,geo_info)
## %
mbr_rect = convertMBRtoRectangle(geo_info.MBR)
my_tess = voronoicells(Vector([i for i in genome]),mbr_rect)
##%
#

#function create_rasterized_tess(tess::Tessellation{Point2{Float64}},ras::AbstractRaster)
#    tess_poly = GeometryBasics.Polygon.(tess.Cells)#convert voronoi cells to polygons
#    A = copy(ras) .= 0
#    #rasterized_poly =  map((i,poly_i) -> rasterize!(A,poly_i, fill=i, progress=false),enumerate(tess_poly))#rasterize voronoi polygons
#    [rasterize!(A,poly_i, fill=i, progress=false) for (i,poly_i) in enumerate(tess_poly)]#rasterize voronoi polygons
#    return A
#end
#
#function create_rasterized_tess(ind::VoronoiIndividual)
#    tess_poly = Polygon.(ind.tess.Cells)#convert voronoi cells to polygons
#    return create_rasterized_tess(tess_poly)
#end
##%
ras = geo_info.population.pop_points
A = create_rasterized_tess(my_tess,ras)
plot(A)
## %


n_facs = 10 
my_size = size(ras)#get the size of the Raster in x and y
range_tuples = map((i,j) -> (i,j),(1,1),my_size)#create ranes 1_size_x etc. 
ranges = map(x-> UnitRange(x...),range_tuples)#convert to Range objects

point_indicies = map(x-> sample(x,(1,n_facs)),ranges)#sample random points from these range objects for x and y
point_coords = map(i -> ras.dims[i][vec(point_indicies[i])],1:length(ras.dims))#get the coord values form raster index
fac_points = GeometryBasics.Point2.(map((i,j) -> (i,j),point_coords...))#conver to Point2 structs

in_boundary  = x -> in(x,geo_info.border)

in_bounds(fac_place)


## %
ras.dims[1][point_tuples[1]]




## %








#function create_rasterized_tess(tess::Tessellation{Point2{Float64}},ras::AbstractRaster)
#    tess_poly = Polygon.(tess.Cells)#convert voronoi cells to polygons
#    A = copy(ras) .= 0
#    rasterized_poly =  map(poly_i -> rasterize!(A,poly_i, fill=1, progress=false),tess_poly)#rasterize voronoi polygons
#    return A
#end
#
#create_rasterized_tess(my_tess,ras.dims)
##create_rasterized_tess(my_tess)
#
#function create_rasterized_tess(ind::VoronoiIndividual)
#    tess_poly = Polygon.(ind.tess.Cells)#convert voronoi cells to polygons
#    return create_rasterized_tess(tess_poly)
#end
#     rasterized_poly =  map(poly_i -> rasterize!(A,poly_i, fill=1, progress=false),tess_poly)#rasterize voronoi polygons
#     A = zeros(UInt32, dimz; missingval=UInt32(0))
# [rasterize!(A,poly_i, fill=i, progress=false) for (i,(raster_i,poly_i)) in enumerate(zip_a)]
