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
dist = MvNormal([0.5,0.5],[0.1 0.0;0.0 0.1])
raster_geo_info = init_geoinfo(1,1,100,dist)
points_geo_info = init_geoinfo(us_data_dir,pop_dir = "/pop_points/pop_points.csv")
geo_info = raster_geo_info
##%
raster_geo_info = init_geoinfo(1,1,100,dist)

##%
genome = generate_genome(geo_info.population,geo_info,n_facs = 10)
#make_voronoi_individual(genome,get_fitness,geo_info)
## %
mbr_rect = convertMBRtoRectangle(geo_info.MBR)
my_tess = voronoicells(Vector([i for i in genome]),mbr_rect)
tess_poly = GeometryBasics.Polygon.(my_tess.Cells)
tess_mp = MultiPolygon(tess_poly)
##%
voronoi_population(genome,geo_info.population,my_tess)
##%
my_beta = 1
dist_metric_cloj = (x,y) -> euclidean_exp(x,y;beta = my_beta) 
fitness_function = (x,y,z) -> mean_dist(x,y,z,dist_metric_cloj)
## %
fitness_function(geo_info.population.pop_points,tess_mp,genome)
voronoiperimeters(my_tess)
voronoiarea(my_tess)
voronoi_population(genome,geo_info.population,my_tess)
## %
v_ind = make_voronoi_individual(genome,fitness_function,geo_info)

ras = geo_info.population.pop_points
ras_tess = create_rasterized_tess(tess_mp,ras)

ras_tess = ras_tess .|> i -> replace_missing(i,0.0)

## % 
function collect_voronoi_raster(rass_tess)
    ras_tess_changed = map(i -> ras_tess[i]*i,1:length(ras_tess))
    voronoi_tess_ras = mosaic(first,ras_tess_changed)
return voronoi_tess_ras

##%
plot(voronoi_tess_ras)

##%
my_dist_mat = dist_mat(ras,tess_mp,genome,euclidean_exp)
## % 
#    my_size = size(ras)#get the size of the Raster in x and y
#    range_tuples = map((i,j) -> (i,j),(1,1),my_size)#create ranes 1_size_x etc. 
#    ranges = map(x-> UnitRange(x...),range_tuples)#convert to Range objects
#
#    in_bounds_index_pairs = []
#    #@generate facility locatons and ensure theyare in bounds, iterate until we have n_facs
#    while length(in_bounds_index_pairs) < n_facs
#        point_indicies = map(x-> sample(x,(1,n_facs)),ranges)#sample random points from these range objects for x and y
#        indices_pairs = map((i,j) -> [i,j],point_indicies...)#get pairs of indicie values
#        is_in_bounds = map(index_i -> ras[index_i...] != MISSING_VAL,indices_pairs)#is the proposed facility in bounds?
#        append!(in_bounds_index_pairs, indices_pairs[is_in_bounds])
#    end
#    in_bounds_index_pairs  =   in_bounds_index_pairs[:n_facs]#cut to exactly n facs
#
#    in_bound_hcat = hcat(in_bounds_index_pairs...)#create a 2d matrix
#    in_bound_seperated = [in_bound_hcat[i,:] for i in 1:size(in_bound_hcat,1)]#and break it into X and Y arrays
#    point_coords = map(i -> ras.dims[i][vec(in_bound_seperated[i])],1:length(ras.dims))#get the coord values form raster index
#    point_pairs = map((i,j) -> [i,j],point_coords...)
#    fac_points = GeometryBasics.Point2.(point_pairs)#conver to Point2 structs
#    return fac_points
#end

#function voronoi_population(ras::AbstractRaster,tess::Tessellation)
#    B = create_rasterized_tess(tess,ras)
#    masked_pop = mask_raster_vec(ras,B)
#    cell_pops = sum.(masked_pop .|> i -> replace_missing(i,0.0))#replace the missing values with 0 and sum reach array
#    return cell_pops
#end

## %


#
#function euclidean_dist(p1,p2)
#    return √(sum((p1.-p2).^2))
#end
#
#function euclidean_exp(p1,p2;beta = 1)
#    return euclidean_dist(p1,p2)^beta
#end
#
#tess_mp = MultiPolygon(Polygon.(tess.Cells))
#dist_closure(x) = euclidean_exp(fac_loc,x)
#
#function dist_mat(ras::AbstractRaster,my_genome::AbstractArray,dist_func::Function)
#    points_ras = points(ras)|> collect |> Raster#get Raster with coordinates of each point
#    masked_cells = map(cell_i -> mask(points_ras;with =cell_i),tess_mp)#create a list of Raster masks for each voronoi cell
#    dist_rasters = [(x -> dist_func(genome[i],x)).(masked_cell_i) for (i,masked_cell_i) in enumerate(masked_cells)]#calculate the distance between each pixel in the masked region on the corresponding facility, save this as the pixel value in a new raster
#    dist_series = RasterSeries(dist_rasters,Ti(1:length(dist_rasters)))##combine these into a RasterSeries
#    my_mosaic = mosaic(first,dist_series)#strich them together into a single raster
#    return my_mosaic
#end


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
#
## % 
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
#
#
#end
#function mean_dist(ras::AbstractRaster,tess_mp::AbstractArray,genome::AbstractArray,dist_func::Function)
#    dist_ras = dist_mat(ras,tess_mp,genome,dist_func)#calculate distance to nearest facility for each pixel 
#    pop_dist_ras = pop_ras .* dist_ras#weight by population
#    mean_distance = sum(zonal(sum,pop_dist_ras,of = tess_poly))#calculate for each zone and sum
#    return mean_weighted_distance
