########  Packages ##############
#using Pkg
#Pkg.add("Symbolics")
using LinearAlgebra, SparseArrays
using Arpack # ARnoldi PACKage: library in FORTRAN 77 for solving large scale eigenvalue problems
using Plots
using RootsAndPoles
using CSV, DataFrames
using Base.Threads
using QuadGK
using FLoops, Base.Iterators
using Distributions, Random
using JLD2
using Zygote, DualNumbers
using FiniteDiff
#using  Symbolics
using ForwardDiff
using Printf 
using SymPy
using DelimitedFiles



#########    Graph Definition   #############
# Vertices in the spot
#const vertIn = [[0.0, 1.0], [-1.2, 0.0], [0.2, -1.2], [1.5, 0.1]];
const vertIn = [row |> Vector for row in 
    eachrow(CSV.read(joinpath(dirname(@__FILE__), "InputData/vertIn.csv"), DataFrame, header=false))];
# Vertices out of the spot
#const vertOut = [[2, .2],[-1.5,0],[-1.6,-.5],[0.0,1.3],[0.1,1.4]];
const vertOut = [row |> Vector for row in 
    eachrow(CSV.read(joinpath(dirname(@__FILE__), "InputData/vertOut.csv"), DataFrame, header=false))];
# Internal Edges
#const fibersIn = [[1, 2], [1, 4], [2, 3], [2, 4], [3, 4]];
const fibersIn = [row |> Vector for row in 
    eachrow(CSV.read(joinpath(dirname(@__FILE__), "InputData/fibersIn.csv"), DataFrame, header=false))];
# External Edges
#const fibersOut = [[1,8],[1,9],[2,6],[2,7],[4,5]];
const fibersOut = [row |> Vector for row in 
    eachrow(CSV.read(joinpath(dirname(@__FILE__), "InputData/fibersOut.csv"), DataFrame, header=false))];
# check order
function check_order(fibersIn)
    for k in 1:length(fibersIn)-1
        if fibersIn[k][1] > fibersIn[k+1][1] || 
            (fibersIn[k][1] == fibersIn[k+1][1] && fibersIn[k][2] >= fibersIn[k+1][2])
            error("The array is not correctly ordered: element $k ($(fibersIn[k])) precedes element $(k+1) ($(fibersIn[k+1]))")
        end
    end
end
check_order(fibersIn)
check_order(fibersOut)

####### Derived Parameters
function derived_params(vertIn,vertOut,fibersIn,fibersOut)

    # Vertices
    NVIn = length(vertIn);
    NVOut = length(vertOut);
    vert = vcat(vertIn, vertOut);
    NV = NVIn + NVOut;

    # Edges
    NEI = length(fibersIn);
    NEO = length(fibersOut);
    fibers = vcat(fibersIn, fibersOut);
    NE = length(fibers);
    fibersLength = [norm(vert[fibers[i][1]] - vert[fibers[i][2]]) for i in 1:NE];  
    fibersPoints = [[vert[fibers[h][1]], vert[fibers[h][2]]] for h in 1:NE]; # couple  of  points

    return  NVIn,  NVOut, vert, NV, NEI, NEO, fibers, NE, fibersLength, fibersPoints

end
const NVIn,  NVOut, vert, NV, NEI, NEO, fibers, NE, fibersLength, fibersPoints = 
    derived_params(vertIn,vertOut,fibersIn,fibersOut);



###### Refractive Indices ##############
const nij = fill(1.5, NE) ; # refractive indices
const delta_pump = fill(1.0, NE)

###### Matrix #####
function LMat_sym(k, vertIn = vertIn, vertOut = vertOut,fibersIn = fibersIn, fibersOut = fibersOut,
    NVIn = NVIn, NVOut = NVOut, vert = vert, NV = NV,NEI = NEI, NEO = NEO, fibers = fibers, 
    NE = NE, fibersLength = fibersLength, fibersPoints = fibersPoints)

    
    #= PosVertFirst[i] refers to the i-node. It is the list of positions 
    of the form (i,j) in fibersIn (internal fibers) =#
    PosVertFirst = [findall(x -> x == i, [pair[1] for pair in fibersIn]) for i in 1:NVIn];
    #= PosVertSecond[i] refers to the i-node . It is the list of positions
    of the form (j, i) in fibersIn (internal fibers) =#
    PosVertSecond = [findall(x -> x == i, [pair[2] for pair in fibersIn]) for i in 1:NVIn];
    PosVert = [vcat(PosVertFirst[i], PosVertSecond[i]) for i in 1:NVIn];
    NPosVertFirst = [length(PosVertFirst[i]) for i in 1:NVIn];
    NPosVertSecond = [length(PosVertSecond[i]) for i in 1:NVIn];
    #= PosVertFiberOut[[i]] refers to the i-node . It is the list of 
    positions of the form (i,j) where j is an external node =#
    PosVertFiberOut = [findall(x -> x == i, [pair[1] for pair in fibersOut]) for i in 1:NVIn];
    NPosVertFiberOut = [length(PosVertFiberOut[i]) for i in 1:NVIn];
    
    # TO TEST USE LMat = Array{Any}(undef, NE + NEI, NE + NEI)
    #LMat = spzeros(Complex{Float64}, NE + NEI, NE + NEI)
    #LMat = spzeros(Complex{typeof(k)}, NE + NEI, NE + NEI)
    LMat = zeros(Sym, NE + NEI, NE + NEI)
     
    for i in eachindex(vertIn)
        # contiunuity eqs. internal fibers
        for j in 1:(NPosVertFirst[i] + NPosVertSecond[i] - 1)
            index1 = sum([NPosVertFirst[l] + NPosVertSecond[l] + NPosVertFiberOut[l] for l in 1:i-1]; init=0) + j;
            if ( 0 < index1 <= NE + NEI)
                index2 = PosVert[i][1];
                if (0 < index2 <= NE + NEI )
                    LMat[index1, index2] = exp(im * k * nij[PosVert[i][1]] * fibersLength[PosVert[i][1]] *
                    Int(i == fibersIn[PosVert[i][1]][2]))
                end
                index2 = PosVert[i][1] + NEI + NEO;
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = exp(im * k * nij[PosVert[i][1]] * fibersLength[PosVert[i][1]] * 
                    Int(i == fibersIn[PosVert[i][1]][1]))
                end
                index2 = PosVert[i][j + 1];
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -exp(im * k * nij[PosVert[i][j + 1]] * fibersLength[PosVert[i][j + 1]] *
                    Int(i == fibersIn[PosVert[i][j + 1]][2]))
                end
                index2 = PosVert[i][j + 1] + NEI + NEO;
                if ( 0 < index1 <= NE + NEI && 0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -exp(im * k * nij[PosVert[i][j + 1]] * fibersLength[PosVert[i][j + 1]] * 
                    Int(i == fibersIn[PosVert[i][j + 1]][1]))
                end
            end
        end
        # contiunuity eqs. external fibers
        for j in 1:NPosVertFiberOut[i]
            index1 = sum(NPosVertFirst[l] + NPosVertSecond[l] + 
            NPosVertFiberOut[l] for l in 1:i-1; init=0) + NPosVertFirst[i] + NPosVertSecond[i] - 1 + j
            if ( 0 < index1 <= NE + NEI)
                index2 = PosVert[i][1]
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = exp(im * k * nij[PosVert[i][1]] *
                    fibersLength[PosVert[i][1]] * Int(i == fibersIn[PosVert[i][1]][2]))
                end
                index2 = PosVert[i][1] + NEI + NEO
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = exp(im * k * nij[PosVert[i][1]] *
                    fibersLength[PosVert[i][1]] * Int(i == fibersIn[PosVert[i][1]][1]))
                end
                index2 = PosVertFiberOut[i][j] + NEI
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -1.0
                end
            end
        end
        # charge conservation
        for j = 1:NPosVertFirst[i]
            index1 = sum(NPosVertFirst[l] + NPosVertSecond[l] +
            NPosVertFiberOut[l] for l in 1:i-1; init=0) + NPosVertFirst[i] + NPosVertSecond[i] + NPosVertFiberOut[i];
            if 0 < index1 <= NE + NEI
                index2 = PosVertFirst[i][j];
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = nij[PosVertFirst[i][j]]
                end
                index2 = PosVertFirst[i][j] + NEI + NEO
                if ( 0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -nij[PosVertFirst[i][j]] *
                    exp(im * k * nij[PosVertFirst[i][j]] * fibersLength[PosVertFirst[i][j]]);
                end
            end
        end
        for j in 1:NPosVertSecond[i]
            index1 = sum(NPosVertFirst[l] + NPosVertSecond[l] +
             NPosVertFiberOut[l] for l in 1:i-1; init=0) + NPosVertFirst[i] + NPosVertSecond[i] + NPosVertFiberOut[i];
            if ( 0 < index1 <= NE + NEI )
                index2 = PosVertSecond[i][j]
                if ( 0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -nij[PosVertSecond[i][j]] * 
                    exp(im * k * nij[PosVertSecond[i][j]] * fibersLength[PosVertSecond[i][j]])
                end
                index2 = PosVertSecond[i][j] + NEI + NEO
                if ( 0 < index2 <= NE + NEI )
                    LMat[index1, index2] = nij[PosVertSecond[i][j]]
                end
            end
        end
        for j in 1:NPosVertFiberOut[i]
            index1 = sum(NPosVertFirst[l] + NPosVertSecond[l] + NPosVertFiberOut[l] for l in 1:i-1; init=0) +
            NPosVertFirst[i] + NPosVertSecond[i] + NPosVertFiberOut[i];
            index2 = NEI + PosVertFiberOut[i][j];
            if ( 0 < index1 <= NE + NEI && 0 < index2 <= NE + NEI )
                LMat[index1, index2] = nij[NEI + PosVertFiberOut[i][j]]
            end
        end
    end

    return LMat

end


#############################################



# Symbolic matrix
k = symbols("k")


LMat_sym_expr = 1.0.*LMat_sym(k)  # Generate symbolic matrix based on `k`




# # First step: replace I with im
# matrix_string_im = [replace(string(elem), "I" => "im") for elem in LMat_sym_expr]

# # Second step: replace exactly "0" with "0.0", only for elements that are exactly zero
# matrix_string = [elem == "0" ? replace(elem, "0" => "0.0") : elem for elem in matrix_string_im]




# # Open a file to save the function
# open(joinpath(dirname(@__FILE__), "LMat.jl"), "w") do file
#     # Write the function definition to the file with the correct name `LMat_sym`
#     write(file, "function LMat(k::Complex{Float64})::Matrix{Complex{Float64}}\n")
#     write(file, "    return [\n")
    
#     # Write each element of the matrix in Julia notation
#     for i in 1:size(LMat_sym_expr, 1)
#         # Join the elements of the row with spaces (Julia's matrix notation)
#         row = join([matrix_string[i, j] for j in 1:size(LMat_sym_expr, 2)], " ")
        
#         # Write the row, ensuring no semicolon is added at the end of the last row
#         if i < size(LMat_sym_expr, 1)
#             write(file, "        $row\n")  # Separate rows with semicolon, except the last row
#         else
#             write(file, "        $row\n")  # Last row should not have a semicolon
#         end
#     end
    
#     # Close the function definition
#     write(file, "    ]\n")
#     write(file, "end\n")
# end

#####################################


# Replace I with im in the string representation and generate Julia code for the function
matrix_string = [replace(string(elem), "I" => "im") for elem in LMat_sym_expr]


function extract_sparse_data(matrix_string)
    # Dimension of the matrix
    n = size(matrix_string, 1)
    
    # Initialize arrays to store row indices (I), column indices (J), and values (V)
    I = Int[]  # Row indices
    J = Int[]  # Column indices
    V = Any[]  # Non-zero values (numeric or symbolic)
    
    # Iterate over the matrix to find non-zero elements
    for i in 1:n
        for j in 1:n
            elem = matrix_string[i, j]
            
            # Skip zero elements
            if elem == "0" || elem == "0.0"
                continue
            end
            
            # Append the row index, column index, and value to the arrays
            push!(I, i)
            push!(J, j)
            push!(V, elem)
        end
    end
    
    return I, J, V, n
end





I, J, V, n = extract_sparse_data(matrix_string)

# write sparse matrix
function write_LMat_sparse(I, J, V, n, filename)
    # Open a file to save the function
    open(filename, "w") do file
        # Write the function definition to the file with the correct name `LMat_sparse`
        write(file, "function LMat_sparse(k)\n")
        
        # Preallocating I, J, and V arrays and inserting elements efficiently
        write(file, "    I = Int[$(join(I, ','))]  # Row indices\n")
        write(file, "    J = Int[$(join(J, ','))]  # Column indices\n")
        write(file, "    V = Complex{Float64}[$(join(V, ','))]  # Non-zero values\n")
        
        #  Create a sparse matrix with n rows and columns
        write(file, "    return sparse(I, J, V, $n, $n)\n")
        
        # End func
        write(file, "end\n")
    end
end

write_LMat_sparse(I, J, V, n, joinpath(dirname(@__FILE__), "LMat_sparse.jl"))

# write dense matrix
function write_LMat_dense(I, J, V, n, filename)
    # Open a file to save the function
    open(filename, "w") do file
        # Write the function definition to the file with the correct name `LMat_dense`
        write(file, "function LMat_dense(k)\n")
        
        # Write the matrix initialization (dense matrix with complex zeros)
        write(file, "    # Initialize a $n x $n matrix of complex zeros\n")
        write(file, "    LMat = Matrix{Complex{Float64}}(zeros(Complex{Float64}, $n, $n))\n")
        
        # Write the indices and values arrays
        write(file, "    # Define row indices, column indices, and values\n")
        write(file, "    I = Int[$(join(I, ','))]  # Row indices\n")
        write(file, "    J = Int[$(join(J, ','))]  # Column indices\n")
        write(file, "    V = Complex{Float64}[$(join(V, ','))]  # Non-zero values\n")
        
        # Vectorized assignment to populate the dense matrix
        write(file, "    # Use vectorized assignment to populate the matrix\n")
        write(file, "    LMat[I .+ (J .- 1) * size(LMat, 1)] .= V\n")
        
        # Return the dense matrix
        write(file, "    return LMat\n")
        
        # Close the function
        write(file, "end\n")
    end
end

write_LMat_dense(I, J, V, n, joinpath(dirname(@__FILE__), "LMat_dense.jl"))



#################################


# # Open a file to save the function
# open(joinpath(dirname(@__FILE__), "LMat_sparse_2.jl"), "w") do file
#     # Write the function definition to the file with the correct name `LMat_sym`
#     write(file, "function LMat_sparse_2(k::Complex{Float64})::Matrix{Complex{Float64}}\n")
#     write(file, "    MM = [\n")
    
#     # Write each element of the matrix in Julia notation
#     for i in 1:size(LMat_sym_expr, 1)
#         # Join the elements of the row with spaces (Julia's matrix notation)
#         row = join([matrix_string[i, j] for j in 1:size(LMat_sym_expr, 2)], " ")
        
#         # Write the row, ensuring no semicolon is added at the end of the last row
#         if i < size(LMat_sym_expr, 1)
#             write(file, "        $row\n")  # Separate rows with semicolon, except the last row
#         else
#             write(file, "        $row\n")  # Last row should not have a semicolon
#         end
#     end
    
#     # Close the function definition
#     write(file, "    ]\n")
#     write(file, "    return sparse(MM)\n")
#     write(file, "end\n")
# end