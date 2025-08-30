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
#using  Symbolics, ForwardDiff
using Printf 

using DelimitedFiles
 using SymPy


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

########### Gain medium parameters #########

const k_a = 10.68  # μm^(-1)D0 = mean(all_D0s_th_0) # Let consider D0 = mean value of non-int threshold
const gamma_perp = 0.5 # μm^(-1) (c=1 units)

function nij_g(k,D0)
    n = length(nij)
    nij_g = zeros(Sym, n)
    for i in eachindex(nij)
        nij_g[i] = sqrt(nij[i]^2 + delta_pump[i] * D0*gamma_perp/(k - k_a + im*gamma_perp))
    end
    return nij_g
end
nij_g(2.7+im,0.0)


###### Matrix with gain medium #####
function LMatGain_sym(k,D0)

    
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
     
    for i in 1:NVIn
        # contiunuity eqs. internal fibers
        for j in 1:(NPosVertFirst[i] + NPosVertSecond[i] - 1)
            index1 = sum([NPosVertFirst[l] + NPosVertSecond[l] + NPosVertFiberOut[l] for l in 1:i-1]; init=0) + j;
            if ( 0 < index1 <= NE + NEI)
                index2 = PosVert[i][1];
                if (0 < index2 <= NE + NEI )
                    LMat[index1, index2] = exp(im * k * nij_g(k,D0)[PosVert[i][1]] * fibersLength[PosVert[i][1]] *
                    Int(i == fibersIn[PosVert[i][1]][2]))
                end
                index2 = PosVert[i][1] + NEI + NEO;
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = exp(im * k * nij_g(k,D0)[PosVert[i][1]] * fibersLength[PosVert[i][1]] * 
                    Int(i == fibersIn[PosVert[i][1]][1]))
                end
                index2 = PosVert[i][j + 1];
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -exp(im * k * nij_g(k,D0)[PosVert[i][j + 1]] * fibersLength[PosVert[i][j + 1]] *
                    Int(i == fibersIn[PosVert[i][j + 1]][2]))
                end
                index2 = PosVert[i][j + 1] + NEI + NEO;
                if ( 0 < index1 <= NE + NEI && 0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -exp(im * k * nij_g(k,D0)[PosVert[i][j + 1]] * fibersLength[PosVert[i][j + 1]] * 
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
                    LMat[index1, index2] = exp(im * k * nij_g(k,D0)[PosVert[i][1]] *
                    fibersLength[PosVert[i][1]] * Int(i == fibersIn[PosVert[i][1]][2]))
                end
                index2 = PosVert[i][1] + NEI + NEO
                if (  0 < index2 <= NE + NEI )
                    LMat[index1, index2] = exp(im * k * nij_g(k,D0)[PosVert[i][1]] *
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
                    LMat[index1, index2] = nij_g(k,D0)[PosVertFirst[i][j]]
                end
                index2 = PosVertFirst[i][j] + NEI + NEO
                if ( 0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -nij_g(k,D0)[PosVertFirst[i][j]] *
                    exp(im * k * nij_g(k,D0)[PosVertFirst[i][j]] * fibersLength[PosVertFirst[i][j]]);
                end
            end
        end
        for j in 1:NPosVertSecond[i]
            index1 = sum(NPosVertFirst[l] + NPosVertSecond[l] +
             NPosVertFiberOut[l] for l in 1:i-1; init=0) + NPosVertFirst[i] + NPosVertSecond[i] + NPosVertFiberOut[i];
            if ( 0 < index1 <= NE + NEI )
                index2 = PosVertSecond[i][j]
                if ( 0 < index2 <= NE + NEI )
                    LMat[index1, index2] = -nij_g(k,D0)[PosVertSecond[i][j]] * 
                    exp(im * k * nij_g(k,D0)[PosVertSecond[i][j]] * fibersLength[PosVertSecond[i][j]])
                end
                index2 = PosVertSecond[i][j] + NEI + NEO
                if ( 0 < index2 <= NE + NEI )
                    LMat[index1, index2] = nij_g(k,D0)[PosVertSecond[i][j]]
                end
            end
        end
        for j in 1:NPosVertFiberOut[i]
            index1 = sum(NPosVertFirst[l] + NPosVertSecond[l] + NPosVertFiberOut[l] for l in 1:i-1; init=0) +
            NPosVertFirst[i] + NPosVertSecond[i] + NPosVertFiberOut[i];
            index2 = NEI + PosVertFiberOut[i][j];
            if ( 0 < index1 <= NE + NEI && 0 < index2 <= NE + NEI )
                LMat[index1, index2] = nij_g(k,D0)[NEI + PosVertFiberOut[i][j]]
            end
        end
    end


    return LMat

end

# Symbolic matrix
k = symbols("k")
D0 = symbols("D0")


LMat_sym_expr = 1.0.*LMatGain_sym(k, D0)  # Generate symbolic matrix based on `k`

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


# write dense matrix
function write_LMat_dense(I, J, V, n, filename)
    # Open a file to save the function
    open(filename, "w") do file
        # Write the function definition to the file with the correct name `LMat_dense`
        write(file, "function LMatGain_dense(k, D0)\n")
        
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

write_LMat_dense(I, J, V, n, joinpath(dirname(@__FILE__), "LMatGain_dense.jl"))
