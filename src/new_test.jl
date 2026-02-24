using Distributed


try 
    global num_cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"]) 
	println("using slurm cores")
catch 
    global num_cores = Threads.nthreads() 


@everywhere begin
	using DrWatson
	@quickactivate "HotInverse"
	using Revise
	include(srcdir("HotInverse.jl"))#this line will make all the code available
end




##
