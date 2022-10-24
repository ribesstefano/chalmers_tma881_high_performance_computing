#!/usr/bin/julia
include("/home/hpc2022/diffusion_opencl/diffusion.jl")

if length(ARGS) != 2
  println("Two or three arguments expected")
end

diffusion_const = nothing
num_iter = nothing
for arg in ARGS
  if arg[1] == '-'
    if arg[2] == 'd'
        global diffusion_const = parse(Float64, arg[3:end])
    elseif arg[2] == 'n'
      global num_iter = parse(Int64, arg[3:end])
    else
      println("Unknown argument $arg")
      exit(1)
    end
  end
end

println(diffusion_const, " ", num_iter)

parse_coordinate(s::String) :: Vector{Float32} = [parse(Float32, ss) for ss in split(s)]

initial_values = nothing
open("init") do file
    width, height = split(readline(file)) # Get dimensions from first line
    initial_values = fill(zero(Float64), height, width)
    println(width, " ", height)
end


# heat_diffusion(
#     width::Int,
#     height::Int,
#     initial_value::F,
#     diffusion_constant::F,
#     nmbiterations::Int
#    )