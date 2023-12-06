using Distributed

@everywhere begin
	using DrWatson
	using Pkg 
	Pkg.activate("/gpfs1/home/w/t/wthomps3/CSDS/research/HotInverse")
	#@quickactivate "HotInverse"
	using Revise
	include(srcdir("HotInverse.jl"))#this line will make all the code available
end

try 
	
    global num_cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"]) 
	println("using slurm cores")
catch 
    global num_cores = Threads.nthreads() 
end 

addprocs(num_cores) 
 
println("Number of cores: ", nprocs()) 
println("Number of workers: ", nworkers()) 






a = zeros(10)

Threads.@threads for i = 1:10
	a[i] = Threads.threadid()
end
@show a
