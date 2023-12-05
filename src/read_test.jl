using DrWatson
@quickactivate "HotInverse"
using Shapefile
using Revise
include(srcdir("HotInverse.jl"))#this line will make all the code available
using DataFrames

df = collect_results(datadir("rasters","io_test"))

##
