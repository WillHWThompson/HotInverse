using DrWatson
using Plots
using Reexport
@reexport using VoronoiCells, GeometryBasics,Random,Shapefile,GeoInterface
@reexport using GeometryTypes: in as GT_in
@reexport using NearestNeighbors, StaticArrays, DataFrames,Statistics,StatsBase,Distributions,Rasters   
@reexport using Revise, Infiltrator
@reexport using Logging

#make sure you include Shapefile seperatley

const RAND_INTERVAL = 0.001

@quickactivate "HOTInverse"

#include(srcdir("genome.jl"))
include(srcdir("data_structures","GeoInfo.jl"))
include(srcdir("data_structures","Individual.jl"))
include(srcdir("data_structures","VoronoiIndividual.jl"))

include(srcdir("evolutionary_operators","selection.jl"))
include(srcdir("evolutionary_operators","mutation.jl"))
include(srcdir("evolutionary_operators","constraints.jl"))

include(srcdir("perturbations","catastrophe.jl"))

include(srcdir("utils/","io_utils.jl"))
include(srcdir("utils/","analysis_utils.jl"))
include(srcdir("utils/","run_evo_alg.jl"))


include(srcdir("algorithms/","evo_algorithm.jl"))
include(srcdir("algorithms/","sim_anneal.jl"))



