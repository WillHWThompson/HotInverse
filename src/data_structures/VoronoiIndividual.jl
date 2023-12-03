
struct VoronoiIndividual <: AbstractIndividual
    fitness::Float64
    genome::Vector{Point2{Float64}}
    tess::Tessellation{Point2{Float64}}
    perimeters::Vector{Float64}
    areas::Vector{Float64}
    populations::Vector{Float64}

end


"""
-------------
Distance Metrics
-------------
"""
function euclidean_dist(p1,p2)
    return √(sum((p1.-p2).^2))
end

function euclidean_exp(p1,p2;beta = 1)
    return euclidean_dist(p1,p2)^beta
end

"""
------------------
Raster Methods
------------------
"""

function create_rasterized_tess(poly_list::AbstractArray,ras::AbstractRaster)
    x,y = ras.dims
    rasters_list = [MISSING_VAL*Raster(ones(x,y),missingval =MISSING_VAL ) for i in 1:length(poly_list)]#create empty raster for each polygon
    [rasterize!(rasters_list[i],poly_i, fill=1, progress=false) for (i,poly_i) in enumerate(poly_list)]#rasterize voronoi polygons
    return rasters_list
end


function create_rasterized_tess(tess::Tessellation{Point2{Float64}},ras::AbstractRaster)
    tess_poly = GeometryBasics.Polygon.(tess.Cells)#convert voronoi cells to polygons
    return create_rasterized_tess(tess_poly,ras)
end


function create_rasterized_tess(ind::VoronoiIndividual)
    tess_poly = Polygon.(ind.tess.Cells)#convert voronoi cells to polygons
    return create_rasterized_tess(tess_poly)
end



function mask_raster_vec(ras::AbstractRaster,B::AbstractArray)
    masked_cells = boolmask.(B)#create boolean mask from rasterized polygon
    masked_arr=map(mask_i -> mask(ras;with = mask_i),masked_cells)#mask arrulation raster with rastreized voronoi cells
    return masked_arr
end


"""
------------------
Voronoi Methods
------------------
"""

#function voronoi_population(genome::AbstractArray,population::PopulationRaster,tess::Tessellation)
#    ras = population.pop_points
#    B = create_rasterized_tess(my_tess,ras)
#    masked_pop = mask_raster_vec(ras,B)
#    cell_pops = sum.(masked_pop .|> i -> replace_missing(i,0.0))#replace the hmissing values with 0 and sum reach array
#    return cell_pops
#end

function voronoi_population(genome::AbstractArray,population::PopulationRaster,tess::Tessellation)
    tess_poly = GeometryBasics.Polygon.(tess.Cells)
    cell_pops = zonal(sum,population.pop_points,of = tess_poly)
    return cell_pops
end


function voronoi_population(genome::AbstractArray,population::PopulationPoints,tess::Tessellation)
    idx_list,dist = get_nearest_neighbors(genome,population.pop_points)
    counts =Int.(zeros(maximum(idx_list)))
    map(x ->  counts[x] = get(counts,x,0)+1,idx_list)
    return counts
end



function dist_mat(my_genome::AbstractArray,geo_info::GeoInfo,tess::Tessellation{Point2{Float64}},dist_func::Function)
    ras = geo_info.population.pop_points
    tess_mp = tess.Cells .|> Polygon |> MultiPolygon
    return dist_mat(ras,tess_mp,my_genome,dist_func)
end

function dist_mat(ras::AbstractRaster,tess_mp::AbstractArray,my_genome::AbstractArray,dist_func::Function)
    points_ras = points(ras)|> collect |> Raster#get Raster with coordinates of each point
    masked_cells = map(cell_i -> mask(points_ras;with =cell_i),tess_mp)#create a list of Raster masks for each voronoi cell
    #@infiltrate
    dist_rasters = [(x -> dist_func(my_genome[i],x)).(masked_cell_i) for (i,masked_cell_i) in enumerate(masked_cells)]#calculate the distance between each pixel in the masked region on the corresponding facility, save this as the pixel value in a new raster
    dist_series = RasterSeries(dist_rasters,Ti(1:length(dist_rasters)))##combine these into a RasterSeries
    my_mosaic = mosaic(first,dist_series)#strich them together into a single raster
    return my_mosaic
end


function mean_dist(ras::AbstractRaster,tess_mp::AbstractArray,genome::AbstractArray,dist_func::Function)
    dist_ras = dist_mat(ras,tess_mp,genome,dist_func)#calculate distance to nearest facility for each pixel 
    pop_dist_ras = ras .* dist_ras#weight by population
    mean_distance = sum(zonal(sum,pop_dist_ras,of = tess_mp))#calculate for each zone and sum
    return mean_distance
end



function voronoi_population(idx_list::Vector{Int64})
    """
    calc_population: given a list of the nearest facility for each citizen in the population, returen the total population of each Voronoi Cell 
    input: 
        idx_list::Vector{Int64}: a list of the nearest facility for each individual in the population
    returns: 
        counts::Vector{Int64}: the total number of individuals in each index
    """
    counts =Int.(zeros(maximum(idx_list)))
    map(x ->  counts[x] = get(counts,x,0)+1,idx_list)
    return counts
end




function calc_perimeter(vertices)
    """
    calc_perimeter: take in a list of veritices of a polygon and returns a the perimiter of the shape
    """
    perimeter = 0.0
    # Loop through each pair of adjacent vertices and add the distance between them to the perimeter
    n = length(vertices)
    @inbounds for i in eachindex(vertices)
        j = mod(i, n) + 1
        dx = vertices[j][1] - vertices[i][1]
        dy = vertices[j][2] - vertices[i][2]
        perimeter += sqrt(dx^2 + dy^2)
    end
    return perimeter
end


function voronoiperimeters(tess::Tessellation{Point2{Float64}})
    """
    get_cell_perimeters: take in a Tesselation object from VoronoiCells.jl and returns the perimiters of each cell in a list
    input: 
        tess::Tesselation, the Voronoi Tesselation object returned by vornoitess()
    """
    map(x ->calc_perimeter(x),tess.Cells)
end


function make_voronoi_individual(genome,fitness_function::Function,border::Shapefile.Polygon)

    #@infiltrate
    polygon = convertShapefileMBRtoRectangle(border.MBR)
    my_tess = voronoicells(Vector([i for i in genome]),polygon)
    idxs,dist = get_nearest_neighbors(genome,geo_info.pop_points)

    VoronoiIndividual(
        get_fitness(idxs,dist),
        genome,
        my_tess,
        voronoiperimeters(my_tess),
        voronoiarea(my_tess),
        voronoipopulation(idxs)
    )
end


function make_voronoi_individual(genome::AbstractArray,fitness_function::Function,geo_info::GeoInfo)
    #@infiltrate
    mbr_rect = convertMBRtoRectangle(geo_info.MBR)
    tess_poly = 0

    try
        tess_poly  = voronoicells(Vector([i for i in genome]),mbr_rect)
    catch e
        println("error")
        @infiltrate
        println("error for real!")
    end


    tess_mp = MultiPolygon(Polygon.(tess_poly.Cells))

    VoronoiIndividual(
        fitness_function(geo_info.population.pop_points,tess_mp,genome),
        genome,
        tess_poly,
        voronoiperimeters(tess_poly),
        voronoiarea(tess_poly),
        voronoi_population(genome,geo_info.population,tess_poly)
    )
end


function make_voronoi_individual(generate_genome_function::Function,fitness_function::Function,geo_info::GeoInfo)
    #@infiltrate
    genome = generate_genome_function()
    mbr_rect = convertMBRtoRectangle(geo_info.MBR)
    tess_poly  = voronoicells(Vector([i for i in genome]),mbr_rect)
    tess_mp = MultiPolygon(Polygon.(tess_poly.Cells))

    VoronoiIndividual(
        fitness_function(geo_info.population.pop_points,tess_mp,genome),
        genome,
        tess_poly,
        voronoiperimeters(tess_poly),
        voronoiarea(tess_poly),
        voronoi_population(genome,geo_info.population,tess_poly)
    )
end




