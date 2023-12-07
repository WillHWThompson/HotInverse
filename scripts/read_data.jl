using Distributed

using DrWatson
@quickactivate "HotInverse"



@everywhere begin 
    #using Shapefile
    using Revise
    include(srcdir("HotInverse.jl"))#this line will make all the code available
end
#
##
try 
    global num_cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"]) 
catch 
    global num_cores = Threads.nthreads() 
end 
addprocs(num_cores) 
 
println("Number of cores: ", nprocs()) 
println("Number of workers: ", nworkers()) 



Threads.@threads for i = 1:10
    println(Threads.threadid())
end

##
df = collect_results(datadir("rasters","io_test"))
##
ind = df[1,:best_ind_over_time]
##
