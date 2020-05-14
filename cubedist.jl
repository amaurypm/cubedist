#!/usr/bin/env julia

using ArgParse
using Images
using Printf
using Distances

function parse_commandline()
    s = ArgParseSettings()
    s.description = "Calculate the pairwise distances (rmsd) between a set of electrostatic maps, contained in Gaussian cube files, and report these distances in CSV and a Mega-compatible formats."
    s.version = "1.0"
    s.add_version = true

    @add_arg_table! s begin
        "--output", "-o"
            help = "output files base name"
            default = "map_distances_rmsd"

        "map"
            help = "electrostatic potential maps"
            required = true
            nargs = '*'

    end

    return parse_args(s)
end

struct GaussianCube
    header::String
    size::Tuple{Int, Int, Int}
    voxels::Array{Float64, 3}
end

"""
    parse_cube(filename::String)

Parse a Gaussian cube phimap file, `filename`, containg a electrostatic potential map and return a proper GaussianCube object.
"""
function parse_cube(filename::String)
    counter = 0
    header = ""
    (x, y, z) = (0, 0, 0)
    natoms = 1
    voxels = Float64[]
    open(filename, "r") do input_file
        for line in eachline(input_file)
            counter += 1
            if counter <= 6 + natoms
                header *= (line * "\n")
                if counter == 3
                    natoms = parse(Int, split(line)[1])
                elseif counter == 4
                    x = parse(Int, split(line)[1])
                elseif counter == 5
                    y = parse(Int, split(line)[1])
                elseif counter == 6
                    z = parse(Int, split(line)[1])
                end

            else
                line = strip(line)
                for field in split(line)
                    push!(voxels, parse(Float64, field))
                end
            end
        end
    end
    if x*y*z != length(voxels)
        write(stderr, "ERROR: Number of voxels read from $filename does not correspond with map dimensions.")
        exit(1)
    end
    voxels = reshape(voxels, x, y, z)
    GaussianCube(header, (x,y,z), voxels)
end

function distance(map1::GaussianCube, map2::GaussianCube)
    rmsd1 = rmsd2 = 0
    if map1.size != map2.size
        rmsd1 = rmsd(map1.voxels, imresize(map2.voxels, map1.size))
        rmsd2 = rmsd(map2.voxels, imresize(map1.voxels, map2.size))
    else
        rmsd1 = rmsd2 = rmsd(map1.voxels, map2.voxels)
    end
    (rmsd1 + rmsd2)/2
end

function write_csv(filename::String, names::Array{String, 1}, mat::Array{Float64,2})
    open(filename, "w") do output_file
        write(output_file, "maps")
        map_names::Array{String, 1} = []
        for map_name in names
            map_name = splitext(basename(map_name))[1]
            push!(map_names, map_name)
            write(output_file, ",$map_name")
        end
        write(output_file, "\n")
        for i in 1:size(mat)[1]
            for j in 1:size(mat)[2]
                if j == 1
                    write(output_file, "$(map_names[i])")
                end
                if i > j
                    write(output_file, ",$(mat[i,j])")
                else
                    write(output_file, ",")
                end
            end
            write(output_file, "\n")
        end
    end
end

function write_meg(filename::String, names::Array{String,1}, mat::Array{Float64,2})
    open(filename, "w") do output_file
        write(output_file, "#mega\n")
        write(output_file, "!Title: Electrostatic map distances;\n")
        write(output_file, "!Format DataType=Distance DataFormat=LowerLeft NTaxa=$(length(names));\n")
        write(output_file, "!Description\n")
        write(output_file, "\tDistance (rmsd) bewteen electrostatic maps calculated with cubdedist\n;\n\n")

        map_names::Array{String, 1} = []
        for map_name in names
            map_name = splitext(basename(map_name))[1]
            push!(map_names, map_name)
        end

        for (i, map_name) in enumerate(map_names)
            write(output_file, "[$i] #$map_name\n")
        end

        write(output_file, "\n[     ")
        for i in 1:length(map_names)
            @printf(output_file, "%9d", i)
        end
        write(output_file, "  ]\n")

        for i in 1:size(mat)[1]
            for j in 1:size(mat)[2]
                if j == 1
                    write(output_file, "[$i]   ")
                end
                if i > j
                    @printf(output_file, "%9.3g", mat[i,j])
                else
                    write(output_file, repeat(" ", 9))
                end
            end
            write(output_file, "\n")
        end
    end
end


function main()
    parsed_args = parse_commandline()
    unique_files::Array{String, 1} = unique(parsed_args["map"])
    n = length(unique_files)
    if n < 2
        write(stderr, "ERROR: At least two distinct files are required.")
        exit(1)
    end

    rmsd_mat = zeros(n, n)

    cube1 = nothing
    cube2 = nothing

    # These cicle has to transverse the lower left corner of the matrix, excluding the diagonal.
    for i in 1:n
        for j in 1:n
            if i > j
                if isfile(unique_files[i]) && isfile(unique_files[j])
                    try
                        cube1 = parse_cube(unique_files[i])
                        cube2 = parse_cube(unique_files[j])
                    catch
                        write(stderr, "ERROR: Couldn't parse provided map files, please check they are proper Gaussian cube files.")
                        exit(1)
                    end

                    rmsd_mat[i, j] = distance(cube1, cube2)
                else
                    write(stderr, "ERROR: At least one of the provided file names does not correspond with a real file.\n")
                    exit(1)
                end
            end
        end
    end

    write_csv(parsed_args["output"]*".csv", unique_files, rmsd_mat)
    write_meg(parsed_args["output"]*".meg", unique_files, rmsd_mat)
end

main()
