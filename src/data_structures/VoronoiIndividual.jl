
struct VoronoiIndividual <: AbstractIndividual
    fitness::Float64
    genome::Vector{Point2{Float64}}
    tess::Tessellation{Point2{Float64}}
    perimeters::Vector{Float64}
    areas::Vector{Float64}
    populations::Vector{Int64}

end




"""
------------------
Raster Methods
------------------
"""

function create_rasterized_tess(tess::Tessellation{Point2{Float64}},ras::AbstractRaster)
    tess_poly = GeometryBasics.Polygon.(tess.Cells)#convert voronoi cells to polygons
    A = copy(ras) .= 0
    #rasterized_poly =  map((i,poly_i) -> rasterize!(A,poly_i, fill=i, progress=false),enumerate(tess_poly))#rasterize voronoi polygons
    [rasterize!(A,poly_i, fill=i, progress=false) for (i,poly_i) in enumerate(tess_poly)]#rasterize voronoi polygons
    return A
end

function create_rasterized_tess(ind::VoronoiIndividual)
    tess_poly = Polygon.(ind.tess.Cells)#convert voronoi cells to polygons
    return create_rasterized_tess(tess_poly)
end




"""
------------------
Voronoi Methods
------------------
"""

function voronoipopulation(idx_list::Vector{Int64})
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


function make_voronoi_individual(genome,fitness_function::Function,geo_info::GeoInfo)

    #@infiltrate
    mbr_rect = convertMBRtoRectangle(geo_info.MBR)
    my_tess = voronoicells(Vector([i for i in genome]),mbr_rect)
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

