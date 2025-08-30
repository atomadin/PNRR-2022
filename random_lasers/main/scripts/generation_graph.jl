# !
# !    Copyright (C) 2025 Camillo Tassi
# !
# !    This program is free software: you can redistribute it and/or modify
# !    it under the terms of the GNU General Public License as published by
# !    the Free Software Foundation, either version 3 of the License, or
# !    (at your option) any later version.
# !
# !    This program is distributed in the hope that it will be useful,
# !    but WITHOUT ANY WARRANTY; without even the implied warranty of
# !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# !    GNU General Public License for more details.
# !
# !    You should have received a copy of the GNU General Public License
# !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# !



using Random
using Plots
using LinearAlgebra
using CSV
using SymPy
using DataFrames
using StatsBase  
#using FLoops
using Printf  # To format directory names with leading zeros
#using ProgressMeter # For progress bar


function gen_graph(target_dir)

    ###### Data folder
    # Define the path of the directory
    path = joinpath(@__DIR__, target_dir,"InputData")
    path_0 = joinpath(@__DIR__, target_dir)

    # Check if the directory exists
    if !isdir(path_0)
        mkdir(path_0)  # Create the directory if it does not exist
    end

    # Check if the directory exists
    if !isdir(path)
        mkdir(path)  # Create the directory if it does not exist
    end


    # Parameters
    RadSpot = 40.0   # radius of the pump laser spot in μm
    NLin = 18        # number of nanofibers
    tLengLines = 500 # length of the lines
    nanodiam = 0.500 # diameter of the nanofiber in nm

    # Generate NLin random points inside a circle of radius RadSpot
    rho = RadSpot .* sqrt.(rand(NLin)) 
    theta = 2 * π * rand(NLin)
    PointsDisk0X = rho .* cos.(theta)
    PointsDisk0Y = rho .* sin.(theta)
    PointsDisk0 = hcat(PointsDisk0X, PointsDisk0Y)


    ##################################################

    # Generate random angles for the segments
    theta_segments = -π/2 .+ π * rand(NLin)

    # Define the parametric equation for the fiber lines
    function FiberLines(t, index)
        x_direction = cos(theta_segments[index])
        y_direction = sin(theta_segments[index])
        return PointsDisk0[index, :] .+ t .* [x_direction, y_direction]
    end



    ##################################################

    # Data Analysis - Vertices
    # Intersection Points and knots. The point resolution is nanodiam (fiber diameter)

    # Define symbolic variables t and s
    t, s = symbols("t s")

    FibIFibJtsPInters = [] #[[[0, 0],[0.0,0.0],[0.0,0.0]] for _ in 1:(NLin*(NLin-1)/2)] [[line1,line2], [t,s], point]

    # Loop to calculate FibIFibJtsPInters
    for i in 1:(NLin - 1)
        for j in (i + 1):NLin
            # Solve FiberLines(t, i) == FiberLines(s, j)
            sol = solve(FiberLines(t, i) .- FiberLines(s, j), [t, s])  # Call FiberLines as a function
            sol_list = Float64.([sol[t], sol[s]])
            point_int = FiberLines(sol_list[1], i)
            # Append the result to FibIFibJtsPInters
            push!(FibIFibJtsPInters, ([i, j], [sol_list[1], sol_list[2]], point_int))

        end
    end


    # Function to round a point to the nearest multiple of nanodiam
    function round_to_nanodiam(point, nanodiam)
        return round.(point ./ nanodiam) .* nanodiam
    end

    # Iterate over FibIFibJtsPInters and round the intersection points to multiples of nanodiam
    FibIFibJtsPInters = [
        (
            FibIFibJtsPInters[i][1],  # [i, j] - fiber pair
            FibIFibJtsPInters[i][2],  # [t, s] - parameter values
            round_to_nanodiam(FibIFibJtsPInters[i][3], nanodiam)  # Round (x, y) to nearest nanodiam multiple
        )
        for i in 1:length(FibIFibJtsPInters)
    ]

    FibIFibJtsPInters[1]

    # Extract points from FibIFibJtsPInters (third element is the point [x, y])
    PInters = [FibIFibJtsPInters[i][3] for i in 1:length(FibIFibJtsPInters)]  # Points

    # Select only points inside the spot (within RadSpot - nanodiam)
    PIntersSpot = filter(p -> norm(p) <= (RadSpot - nanodiam), PInters)  # spot internal points

    # Remove duplicate points
    Knots = unique(PIntersSpot)  # vertices

    # Convert Knots (Vector of points) to a DataFrame with two columns for x and y
    Knots_df = DataFrame(x = [k[1] for k in Knots], y = [k[2] for k in Knots])

    # Export the vertices to a CSV file in the relative directory
    #CSV.write(joinpath(@__DIR__, "InputData/vertIn.csv"), Knots_df, header=false)
    CSV.write(joinpath(path, "vertIn.csv"), Knots_df, header=false)

    # Number of vertices
    NKnots = length(Knots)  # number of vertices


    #####################################

    # Define symbolic variable for t
    t = symbols("t")

    # Initialize the result list to store the external points (line, t, point)
    FibItPOut = []

    # Iterate over all fibers to find intersections with the circle of radius RadSpot
    for i in 1:NLin
        # Parametric equations for fiber i
        fiber = FiberLines(t, i)

        # Set up the equation: x^2 + y^2 == RadSpot^2
        x_eq = fiber[1]^2  # x coordinate
        y_eq = fiber[2]^2  # y coordinate
        eqn = Eq(x_eq + y_eq, RadSpot^2)  # Circle equation

        # Solve the equation for t
        solutions = solve(eqn, t)

        # Check if we have valid solutions
        if !isempty(solutions)
            for sol in solutions
                t_val = sol  # Extract the value of t
                point = Float64.(FiberLines(t_val, i))  # Calculate the point on the fiber
                push!(FibItPOut, (i, t_val, point))  # Store {i, t, (x, y)}
            end
        else
            println("No valid solution for fiber $i")
        end
    end


    FibItPOut


    # Extract the third element from each tuple in FibItPOut
    POut = [point[3] for point in FibItPOut]

    # Remove duplicates from POut
    KnotsOut = unique(POut)

    # Convert Knots (Vector of points) to a DataFrame with two columns for x and y
    KnotsOut_df = DataFrame(x = [k[1] for k in KnotsOut], y = [k[2] for k in KnotsOut])

    # Export the vertices to a CSV file in the relative directory
    #CSV.write(joinpath(@__DIR__, "InputData/vertOut.csv"), KnotsOut_df, header=false)
    CSV.write(joinpath(path, "vertOut.csv"), KnotsOut_df, header=false)


    #################################

    # Select the elements where the norm of [x,y] is less than or equal to (RadSpot - nanodiam)
    FibIFibJtsPIntersSpot_1 = filter(x -> norm(x[3]) <= (RadSpot - nanodiam), FibIFibJtsPInters)
    # Initialize FibIFibJtsPIntersSpot as an array of tuples with the desired structure
    FibIFibJtsPIntersSpot = [ 
        (Vector{Int64}(), Vector{Float64}(), Vector{Int64}()) for _ in eachindex(FibIFibJtsPIntersSpot_1)
    ]

    for i in eachindex(FibIFibJtsPIntersSpot_1)
        # Find the index in Knots
        index = findfirst(==(FibIFibJtsPIntersSpot_1[i][3]), Knots)
        
        # Replace (x, y) pair with a singleton vector containing the index in Knots
        FibIFibJtsPIntersSpot[i] = (
            Vector{Int64}(FibIFibJtsPIntersSpot_1[i][1]),  # Convert to Vector{Int64}
            Vector{Float64}(FibIFibJtsPIntersSpot_1[i][2]),  # Convert to Vector{Float64}
            index === nothing ? Vector{Int64}() : [Int64(index)]  # Ensure Vector{Int64}
        )
    end

    FibIFibJtsPIntersSpotAdj = Vector{Tuple{Vector{Int}, Vector{Float64}, Vector{Int}}}()

    for i in eachindex(FibIFibJtsPIntersSpot)
        original = FibIFibJtsPIntersSpot[i]
        
        # Append original and reciprocal terms
        push!(FibIFibJtsPIntersSpotAdj, (
            [original[1][2], original[1][1]],  # Swap (i, j) to (j, i) and convert to Vector
            [original[2][2], original[2][1]],  # Swap (t, s) to (s, t) and convert to Vector
            original[3]                        # Keep {k} as Vector
        ))
    end


    LinIntSpot = []

    for i in 1:NLin
        # Filter elements of FibIFibJtsPIntersSpot where the first component equals i
        filtered_fib_ifib = filter(x -> x[1][1] == i, FibIFibJtsPIntersSpot)
        
        # Filter elements of FibIFibJtsPIntersSpotAdj where the first component equals i
        filtered_fib_ifib_adj = filter(x -> x[1][1] == i, FibIFibJtsPIntersSpotAdj)
        
        # Concatenate the results and add them to LinIntSpot
        push!(LinIntSpot, vcat(filtered_fib_ifib, filtered_fib_ifib_adj))
    end

    LinIntSpotOrd = [
        # Sort each element of LinIntSpot by the first component of the second element
        sort(LinIntSpot[i], by = x -> x[2][1]) for i in 1:length(LinIntSpot)
    ]

    FibersIn = []

    for j in 1:length(LinIntSpotOrd)
        for i in 1:(length(LinIntSpotOrd[j]) - 1)
            # Retrieve the minimum and maximum values of connected knots
            knot1 = LinIntSpotOrd[j][i][3][1]
            knot2 = LinIntSpotOrd[j][i + 1][3][1]
            
            # Add the ordered pair {Min, Max} to FibersIn
            push!(FibersIn, (min(knot1, knot2), max(knot1, knot2)))
        end
    end

    # Sort FibersIn by the values of the pairs
    FibersIn = sort(FibersIn, by = x -> (x[1], x[2]))

    # Remove duplicate pairs from FibersIn
    FibersIn = unique(FibersIn)

    # Keep only the pairs where the two elements are different
    FibersIn_filtered = [pair for pair in FibersIn if pair[1] != pair[2]]

    # Create the DataFrame with two columns "Node1" and "Node2"
    df = DataFrame(Node1 = [x[1] for x in FibersIn_filtered], Node2 = [x[2] for x in FibersIn_filtered])

    # Define the relative path for output file
    #output_path = joinpath(@__DIR__, "InputData", "fibersIn.csv")
    output_path = joinpath(path, "fibersIn.csv")

    # Write the DataFrame to a CSV file without headers
    CSV.write(output_path, df, writeheader=false, header = false)


    ###############


    # Create a copy of FibItPOut to modify it
    FibItPOutSpot = deepcopy(FibItPOut)

    # Update FibItPOutSpot by finding the position in KnotsOut and adjusting by NKnots
    for i in 1:length(FibItPOutSpot)
        FibItPOutSpot[i] = (
            FibItPOutSpot[i][1], 
            FibItPOutSpot[i][2], 
            findfirst(==(FibItPOutSpot[i][3]), KnotsOut) + NKnots
        )
    end

    # Build LinOutSpot based on matching elements in FibItPOutSpot
    LinOutSpot = [filter(x -> x[1] == i, FibItPOutSpot) for i in 1:NLin]

    # Sort LinOutSpot by the second component in each sublist
    LinOutSpotOrd = [sort(LinOutSpot[i], by = x -> x[2]) for i in 1:length(LinIntSpot)]

    # Initialize FibersOut and add ordered pairs
    FibersOut = []

    for i in 1:NLin
        if length(LinIntSpotOrd[i]) > 0
            # Add pairs based on the first and last elements of LinIntSpotOrd and LinOutSpotOrd
            push!(FibersOut, (LinIntSpotOrd[i][1][3][1], LinOutSpotOrd[i][1][3][1]))
            push!(FibersOut, (LinIntSpotOrd[i][end][3][1], LinOutSpotOrd[i][2][3][1]))
        end
    end

    # Sort FibersOut by the pairs and remove duplicates
    FibersOut = unique(sort(FibersOut, by = x -> (x[1], x[2])))

    # Keep only the pairs where the two elements are different
    FibersOut_filtered = [pair for pair in FibersOut if pair[1] != pair[2]]

    # Export FibersOut to a CSV file without headers
    output_path = joinpath(path, "fibersOut.csv")#joinpath(@__DIR__, "InputData", "fibersOut.csv")
    CSV.write(output_path, DataFrame(FibersOut_filtered), writeheader=false)



end



# Main script execution
if length(ARGS) < 1
    println("Usage: julia script.jl <target_dir>")
    exit(1)
end

target_dir = ARGS[1]
println("Generating graph data in directory: $target_dir")
gen_graph(target_dir)







