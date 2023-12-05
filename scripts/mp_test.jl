using Distributed


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
