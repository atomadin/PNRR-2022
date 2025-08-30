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
using  Symbolics, ForwardDiff
using Printf 
#using ProgressMeter
using DelimitedFiles


##### Parameters #########

# Function to read a parameter from a file
function read_parameter(filename, param_name)
    # Read all lines from the file
    lines = readlines(filename)
    
    # Loop through the lines to find the one containing the parameter name
    for line in lines
        if startswith(line, param_name)
            # Split the string based on the '=' symbol to extract the value
            _, value_str = split(line, "=")
            # Trim whitespace and convert the value to a float
            return parse(Float64, strip(value_str))
        end
    end
    # Throw an error if the parameter is not found
    error("Parameter not found")
end

# Example usage to get the value of int_root_atol
#int_root_atol = read_parameter(joinpath(dirname(@__FILE__), "InputData/input_par.txt"), "int_root_atol")
int_root_atol = 0.001;  #0.01
int_root_check_tol = 0.1;  #0.4



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
check_order(fibersIn);
check_order(fibersOut);

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

# Write to the file
file_path = joinpath(@__DIR__, "out_numb.txt");
open(file_path, "w") do file
    write(file, "$NVIn\n");
    write(file, "$NE\n");
end


# Number of bins using Sturges' rule
bins_sturges = ceil(Int, log2(length(fibersLength)) + 1);
# Plot histogram
histogram(
    fibersLength,
    bins=bins_sturges,
    title="Fiber Lengths, #bins = $bins_sturges, NE = $NE, NV = $NV",
    xlabel="Length",
    ylabel="Frequency",
    legend=false
);
# Save the plot
output_path = joinpath(dirname(@__FILE__),"fiber_lengths_histogram.pdf");  # Change this to your preferred file name and path
savefig(output_path);

###### Refractive Indices ##############
const nij = fill(1.5, NE) ; # refractive indices
const delta_pump = fill(1.0, NE);

##### Plotting Graph ##################
function plotFiberGraph(fibersPoints, vertIn, vertOut)
    
    # Create the plot
    plt = plot()

    # Add lines representing the fiber segments
    connected_points = Set()  # Set to store connected vertices
    for points in fibersPoints
        # Add the line segment
        plot!([points[1][1], points[2][1]], [points[1][2], points[2][2]], 
            color=:black, linewidth=1, label="")
        
        # Add points to the connected set
        push!(connected_points, points[1])
        push!(connected_points, points[2])
    end

    # Add red points (internal vertices)
    scatter!([point[1] for point in vertIn], [point[2] for point in vertIn], 
        color=:red, markersize=3, label="Int vert")

    # Filter vertOut: only include points connected to fibers
    filtered_vertOut = [point for point in vertOut if point in connected_points]

    # Add green points (external vertices connected to fibers)
    scatter!([point[1] for point in filtered_vertOut], [point[2] for point in filtered_vertOut], 
        color=:green, markersize=3, label="Ext vert")

    # Show the plot
    plot!(legend=:topleft, aspect_ratio=:equal)

    # Save and display the plot
    savefig(plt, joinpath(@__DIR__, "fibers.pdf"))
    display(plt);
end
plotFiberGraph(fibersPoints, vertIn, vertOut);


###### Matrix #####
function LMat_num(k, vertIn = vertIn, vertOut = vertOut,fibersIn = fibersIn, fibersOut = fibersOut,
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
    LMat = spzeros(Complex{Float64}, NE + NEI, NE + NEI)
    #LMat = spzeros(Complex{typeof(k)}, NE + NEI, NE + NEI)
     
    for i in 1:NVIn #eachindex(vertIn)
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

# determinant
function detLMat_num(k)
    Mat = LMat_num(k);

    return det(Mat)
end



################## Simbolics

# Include the generated function file
include(joinpath(dirname(@__FILE__), "LMat_dense.jl"));
# Include the generated function file
include(joinpath(dirname(@__FILE__),"LMat_sparse.jl"));





LMat = LMat_num;
detLMat = detLMat_num;


##### Determing Im(k) region in finding roots #####

# complex domain over which we would like to search roots of k -> detLMat(k)
# Re_k_min = 5.0 + 1.36;  # real part begin
# Re_k_max = 16.36 - 1.36;  # real part end
include(joinpath(dirname(@__FILE__),"InputData/input_data.txt"));
ext_min_roots_distance = 1.0/(maximum(real.(nij).*fibersLength)*2*NV);
# mesh
mesh_factor = 1; #20
r_mesh = ext_min_roots_distance/mesh_factor;

open(joinpath(dirname(@__FILE__), "out_data.txt"), "w") do file
    write(file, "Initial mesh: $r_mesh");
end



# Numerical/Automatic Differentiation
function num_auto_diff(func, x; h=1e-3, meth=5)
    if meth == 1
        return (func(x + h) - func(x - h)) / (2*h)
    elseif meth == 2
        _, back = Zygote.pullback(func, x)
        du, dv = back(1)[1], back(im)[1]
        return (du' + im*dv')/2
    elseif meth == 3
        return dualpart(func(Dual(x,1)))
    elseif meth == 4
        return (-func(x+2*h)+8*func(x+h)-8*func(x-h)+func(x-2*h))/(12*h)
    elseif meth==5
        return FiniteDiff.finite_difference_derivative(func, x)
    end
end





# # Newton-Raphson method with numerical derivative and Im(k)<0
function Newton_Raphson_mod(f, initialGuess, real_min, real_max, r_mesh; tol=1e-6, maxIter=10^3, exp_prec=1, r_mesh_fac_tol=500)
    
    x = initialGuess
    fx = f(x)
    lamdax = minimum(abs.(eigvals(Matrix(LMat(x)))))
    iter=0
    while lamdax > tol && iter < maxIter
        dfx = num_auto_diff(f, x, h=r_mesh*10.0^(-3))
        iter += 1
        if dfx == 0
            x = x - rand()*r_mesh*10.0^(exp_prec)*im
            fx = f(x)
            lamdax = minimum(abs.(eigvals(Matrix(LMat(x)))))
        elseif imag(x  - fx/dfx)>0 && real(x  - fx/dfx) < real_min - r_mesh_fac_tol*r_mesh
            x = (real_min + rand()*r_mesh*10.0^(exp_prec)) - rand()*r_mesh*10.0^(exp_prec)*im
            fx = f(x)
            lamdax = minimum(abs.(eigvals(Matrix(LMat(x)))))
        elseif imag(x  - fx/dfx)>0 && real(x  - fx/dfx) > real_max + r_mesh_fac_tol*r_mesh
            x = (real_max - rand()*r_mesh*10.0^(exp_prec)) - rand()*r_mesh*10.0^(exp_prec)*im
            fx = f(x)
            lamdax = minimum(abs.(eigvals(Matrix(LMat(x)))))
        elseif imag(x  - fx/dfx) > 0
            x = real(x  - fx/dfx) - rand()*r_mesh*10.0^(exp_prec)*im
            fx = f(x)
            lamdax = minimum(abs.(eigvals(Matrix(LMat(x)))))
        elseif real(x  - fx/dfx) < real_min - r_mesh_fac_tol*r_mesh
            x = (real_min + rand()*r_mesh*10.0^(exp_prec)) + imag(x  - fx/dfx)
            fx = f(x)
            lamdax = minimum(abs.(eigvals(Matrix(LMat(x)))))
        elseif real(x  - fx/dfx) > real_max + r_mesh_fac_tol*r_mesh
            x = (real_max - rand()*r_mesh*10.0^(exp_prec)) + imag(x  - fx/dfx)
            fx = f(x)
            lamdax = minimum(abs.(eigvals(Matrix(LMat(x)))))
        else
            x = x  - fx/dfx
            fx = f(x)
            lamdax = minimum(abs.(eigvals(Matrix(LMat(x)))))
        end
    end
    if lamdax < tol #abs(fx) < tol
        return x
    else
        return 
    end
end
# wrong side (test) Newton-Raphson method with numerical derivative and Im(k)<0



# Test to compile
Newton_Raphson_mod(detLMat, 9.6-0.00000000001im, Re_k_min, Re_k_max, r_mesh; tol=10^90, maxIter=1);




# #Determining the region Im(k) to consider


function guess_root(detLMat, Re_k_min, Re_k_max, r_mesh; samp = 100)
    initialGuesses = range(Re_k_min, stop=Re_k_max, length=samp)
    root=Array{Any}(nothing, length(initialGuesses))

    @floop for i in eachindex(initialGuesses)
        root[i] = Newton_Raphson_mod(detLMat, initialGuesses[i]-r_mesh/10*im, Re_k_min, Re_k_max, r_mesh; tol=10^-6, maxIter=50);
        #Newton_Raphson_mod(detLMat, initialGuesses[i], Re_k_min, Re_k_max, r_mesh; tol=1e-1*r_mesh, maxIter=10^3);
    end

    return root
end

##### SAVED DATA ####
# ~ 3 min
root = guess_root(detLMat, Re_k_min, Re_k_max, r_mesh, samp = 10^3);
# Exporting
s_root = something.(root, missing);
s_root_r_i =  Array{Union{Float64, Missing}, 2}(undef, length(s_root), 2);
for i in eachindex(s_root)
    if ismissing(s_root[i])
        s_root_r_i[i, 1] = missing;
        s_root_r_i[i, 2] = missing;
    else
        s_root_r_i[i, 1] = real(s_root[i]);
        s_root_r_i[i, 2] = imag(s_root[i]);
    end
end
s_root_r_i;
CSV.write(joinpath(dirname(@__FILE__), "out_init_guess_roots.csv"), Tables.table(s_root_r_i), writeheader=false);
# # Importing
# df = DataFrame(CSV.File(joinpath(dirname(@__FILE__), "out_init_guess_roots.csv"), header=false))
# function create_complex(real_part, imag_part)
#     if ismissing(real_part) || ismissing(imag_part)
#         return nothing
#     else
#         return real_part + imag_part * im
#     end
# end
# root = [create_complex(df.Column1[i], df.Column2[i]) for i in 1:size(df, 1)]
#### END SAVED DATA ####




if count(x -> x === nothing, root) > 100 #length(initialGuesses)/2
    @warn "Too many fails in Newton_Raphson_mod! Possible error in Im(k) size!";
    
    # Path to the output file
    file_path = joinpath(@__DIR__, "out_data.txt");
    
    # Create the file if it does not exist, and append the warning along with the code line
    open(file_path, "a") do io
        write(io, "Too many fails in Newton_Raphson_mod! Possible error in Im(k) size!\n");
        write(io, "Line: if count(x -> x === nothing, root) > 100 #length(initialGuesses)/2\n");
    end
end


root = filter(x -> x !== nothing, root);

# Plot
rr = real.(root);
ir = imag.(root);
scatter(rr,ir);

function remove_close_values(root, r_mesh)
    eps = r_mesh * 1e-0  # Soglia di prossimità
    cleaned = Complex[]  # Array pulito inizialmente vuoto
    
    for x in root
        if all(abs(x - y) >= eps for y in cleaned)
            push!(cleaned, x)
        end
    end
    
    return cleaned
end

root_cleaned = remove_close_values(root, r_mesh);

# Sort by imaginary part in descending order (from zero to more negative values)
sorted_roots = sort(root_cleaned, by = imag, rev = true);

# Select the first 35 elements, or fewer if the array contains less than 35
top_150_roots = sorted_roots[1:min(221, end)];
top_149_roots = sorted_roots[1:min(220, end)];

min_im_root = minimum(imag.(root_cleaned));
min_im_root_m_1 = minimum(imag.(top_149_roots));
max_im_root = maximum(imag.(top_150_roots));

# min_im_root = minimum(imag.(root_cleaned))
# max_im_root = maximum(imag.(root_cleaned))

# Im(k) region
Im_k_min = mean([min_im_root , min_im_root_m_1]); #minimum([min_im_root, max_im_root*100]) #ext_min_im_k  # imag part begin
#Im_k_min = minimum([min_im_root, max_im_root*100])
Im_k_max = 0.0;#ext_max_im_k  # imag part end



##### Counting Roots ####

# integrand functions
function gamma_rect(t,re_k_min,re_k_max,im_k_min,im_k_max)

        if 0 <= t < 1
            gamma = im*im_k_min + re_k_min + t*(re_k_max-re_k_min)
        elseif 1 <= t < 2
            gamma = re_k_max + im*(im_k_min + (t-1)*(im_k_max - im_k_min))
        elseif 2 <= t < 3
            gamma = im*im_k_max + re_k_max + (t-2)*(re_k_min - re_k_max)
        elseif 3 <= t <= 4
            gamma = re_k_min + im*(im_k_max + (t-3)*(im_k_min - im_k_max))
        else
            gamma = 0
        end
            
        return gamma
     
end
    
function gamma_rect_der(t,re_k_min,re_k_max,im_k_min,im_k_max)
    
        if 0 <= t < 1
            gamma =  (re_k_max-re_k_min)
        elseif 1 <= t < 2
            gamma =  im*( (im_k_max - im_k_min))
        elseif 2 <= t < 3
            gamma = (re_k_min - re_k_max)
        elseif 3 <= t <= 4
            gamma =  im*( (im_k_min - im_k_max))
        else
            gamma = 0
        end
            
        return gamma
     
end

fun_der = k -> num_auto_diff(detLMat, k, h=r_mesh*1e-3);

integrand_fun(t, min_re, max_re, min_im, max_im) = 1/(2*pi*im) * fun_der(gamma_rect(t, min_re, max_re, min_im, max_im)) / 
    detLMat(gamma_rect(t, min_re, max_re, min_im, max_im)) * gamma_rect_der(t, min_re, max_re, min_im, max_im);

############################################################# Remove (?)
#######################################################################




# definition of the partition to avoid that integrals are too close to the roots
function redef_part(Re_k_min, Re_k_max, Im_k_min, Im_k_max, numb_re_int, numb_im_int, max_iter, en_fac_mesh; new_r_mesh=r_mesh)
    
    # uniform definition of a partition of the complex region
    Delta_k_re = (Re_k_max - Re_k_min) / numb_re_int
    Delta_k_im = (Im_k_max - Im_k_min) / numb_im_int
    re_k_ext_or = [Re_k_min + Delta_k_re*i for i in 0:numb_re_int]
    im_k_ext_or = [Im_k_min + Delta_k_im*i for i in 0:numb_im_int]
    
    new_r_mesh=minimum([Delta_k_re, Delta_k_im])/en_fac_mesh # mesh

    #redifinition    re_k_ext_or = [Re_k_min + Delta_k_re*i for i in 0:numb_re_int]
    im_k_ext_or = [Im_k_min + Delta_k_im*i for i in 0:numb_im_int]
    re_k_ext = deepcopy(re_k_ext_or)
    im_k_ext = deepcopy(im_k_ext_or)
    
    #integrand functions
    fun_der_line = k -> num_auto_diff(detLMat, k, h=new_r_mesh*1e-3)
    function gamma_line(t,re_i,im_i,re_f,im_f)
        if 0 <= t <= 1
            gamma = re_i + im*im_i + t * (re_f + im*im_f - re_i - im*im_i)
        else
            gamma = 0
        end
        
        return gamma
    end
    
    function gamma_line_der(t,re_i,im_i,re_f,im_f)
        if 0 <= t <= 1
            gamma = (re_f + im*im_f - re_i - im*im_i)
        else
            gamma = 0
        end
    
        return gamma
    end
    
    integrand_fun_line(t, re_i,im_i,re_f,im_f) = 1/(2*pi*im) * fun_der_line(gamma_line(t, re_i,im_i,re_f,im_f)) / 
    detLMat(gamma_line(t, re_i,im_i,re_f,im_f)) * gamma_line_der(t, re_i,im_i,re_f,im_f)
    
    iter = 1

    # displacing real intervals too close to roots
    
    err_integral_lines_re = Array{Float64}(undef, numb_re_int + 1)
    @floop for i in 1:(numb_re_int + 1)
        _, err_integral_lines_re[i] = quadgk(t -> integrand_fun_line(t, re_k_ext[i], Im_k_min, re_k_ext[i], Im_k_max), 0, 1, atol=int_root_atol)
    end
    
    too_large_err_ind =  findall(x -> x > int_root_check_tol, err_integral_lines_re)
    
    if !isempty(too_large_err_ind)
        for i in too_large_err_ind
            while iter <= max_iter && err_integral_lines_re[i] > int_root_check_tol
                iter = iter + 1;
                re_k_ext[i] = re_k_ext_or[i] + rand(Uniform(-new_r_mesh/en_fac_mesh, new_r_mesh/en_fac_mesh)) #rand(Uniform(-Delta_k_re/2,Delta_k_re/2))
                _, err_integral_lines_re[i] = quadgk(t -> integrand_fun_line(t, re_k_ext[i], Im_k_min, re_k_ext[i], Im_k_max), 0, 1, atol=int_root_atol)
            end
        end
        if err_integral_lines_re[i] > int_root_check_tol
            error("Too many iterations in finding integral lines not too close to roots")
        end
    end
                        
     # displacing im intervals too close to roots
     
     err_integral_lines_im = Array{Float64}(undef, numb_im_int + 1)
     @floop for i in 1:(numb_im_int + 1)
        _, err_integral_lines_im[i] = quadgk(t -> integrand_fun_line(t, Re_k_min ,im_k_ext[i], Re_k_max, im_k_ext[i]), 0, 1, atol=int_root_atol)
    end
    
    too_large_err_ind =  findall(x -> x > int_root_check_tol, err_integral_lines_im)
    
    if !isempty(too_large_err_ind)
        for i in too_large_err_ind
            while iter <= max_iter && err_integral_lines_im[i] > int_root_check_tol
                iter = iter + 1;
                im_k_ext[i] = im_k_ext_or[i] + rand(Uniform(-new_r_mesh/en_fac_mesh, new_r_mesh/en_fac_mesh)) #rand(Uniform(-Delta_k_im/2,Delta_k_im/2))
                _, err_integral_lines_im[i] = quadgk(t -> integrand_fun_line(t, Re_k_min ,im_k_ext[i], Re_k_max, im_k_ext[i]), 0, 1, atol=int_root_atol)
            end
            if err_integral_lines_im[i] > int_root_check_tol
                error("Too many iterations in finding integral lines not too close to roots")
            end
        end
    end

    return re_k_ext, im_k_ext
end


###############################################################
###############################################################


function reg_with_few_sol_old(Re_k_min, Re_k_max, Im_k_min, Im_k_max, numb_roots; max_iter=300, en_fac_mesh=100, ref_fac=1, new_r_mesh=r_mesh, redef=false, parall=true)
    # numb of subregions
    numb_re_int = Int(round(sqrt(numb_roots * ref_fac)))
    numb_im_int = Int(round(sqrt(numb_roots * ref_fac)))

    # mesh factor
    Delta_k_re = (Re_k_max - Re_k_min) / numb_re_int
    Delta_k_im = (Im_k_max - Im_k_min) / numb_im_int


    # redefining regions
    if redef
        re_k_ext_or, im_k_ext_or = redef_part(Re_k_min, Re_k_max, Im_k_min, Im_k_max, numb_re_int, numb_im_int, max_iter, en_fac_mesh)
    else
        re_k_ext_or = Float64[Re_k_min + Delta_k_re * i for i in 0:numb_re_int]
        im_k_ext_or = Float64[Im_k_min + Delta_k_im * j for j in 0:numb_im_int]
    end
    re_k_ext = deepcopy(re_k_ext_or)
    im_k_ext = deepcopy(im_k_ext_or)

    # integrands
    fun_der_int = k -> num_auto_diff(detLMat, k, h=new_r_mesh * 1e-3)
    integrand_fun_int(t, min_re, max_re, min_im, max_im) = 1/(2 * pi * im) * fun_der_int(gamma_rect(t, min_re, max_re, min_im, max_im)) / 
        detLMat(gamma_rect(t, min_re, max_re, min_im, max_im)) * gamma_rect_der(t, min_re, max_re, min_im, max_im)

    # number of sol per region    
    integrals = Array{Complex{Float64}}(undef, numb_re_int, numb_im_int)
    errors = Array{Float64}(undef, numb_re_int, numb_im_int)

    # # Initialize progress bar
    # n_iterations = numb_re_int * numb_im_int
    # progress = Progress(n_iterations, 1)

    if parall
        @floop for (i, j) in product(1:numb_re_int, 1:numb_im_int)
            integrals[i, j], errors[i, j] = quadgk(t -> integrand_fun_int(t, re_k_ext[i], re_k_ext[i+1], im_k_ext[j], im_k_ext[j+1]), 0, 4, atol=int_root_atol)
            # next!(progress)
        end
    else
        for i in 1:numb_re_int
            for j in 1:numb_im_int
                integrals[i, j], errors[i, j] = quadgk(t -> integrand_fun_int(t, re_k_ext[i], re_k_ext[i+1], im_k_ext[j], im_k_ext[j+1]), 0, 4, atol=int_root_atol)
                # next!(progress)
            end
        end
    end

    # consistency test
    for i = 1:numb_re_int
        for j = 1:numb_im_int
            if errors[i, j] > int_root_check_tol || abs(imag(integrals[i, j])) > int_root_check_tol
                @warn "Something wrong with the integrals on ($i,$j)! Try to change ref_fac"
            end
        end
    end

    # # consistency test
    numb_sol_matrix = Int.(round.(real.(integrals)))
    # numb_sol_found = sum(numb_sol_matrix)
    # if numb_sol_found != numb_roots
    #     @warn "Something wrong with the number of solutions: numb sol found $numb_sol_found, numb roots $numb_roots. Maybe the region has changed!"
    # end

    # regions with roots
    max_numb_sol_matrix = maximum(numb_sol_matrix)
    ind_numb_sol = [[] for _ in 1:max_numb_sol_matrix]
    regions_with_roots = [[] for _ in 1:max_numb_sol_matrix]
    for i in 1:max_numb_sol_matrix
        ind_numb_sol[i] = findall(x -> x == i, numb_sol_matrix)
        for j in eachindex(ind_numb_sol[i])
            push!(regions_with_roots[i], [Re_k_min + (ind_numb_sol[i][j][1] - 1) * Delta_k_re, Re_k_min + ind_numb_sol[i][j][1] * Delta_k_re,
                Im_k_min + (ind_numb_sol[i][j][2] - 1) * Delta_k_im, Im_k_min + ind_numb_sol[i][j][2] * Delta_k_im]) # re_min, re_max, im_min, im_max
        end
    end

    return regions_with_roots #, max_numb_sol_matrix
end





function reg_with_few_sol(Re_k_min, Re_k_max, Im_k_min, Im_k_max, numb_roots; max_iter=300, en_fac_mesh=100, ref_fac=1, new_r_mesh=r_mesh, redef=false, parall=true, numb_sol_bound = 500)
    
    # Number of subregions
    numb_re_int = Int(round(sqrt(numb_roots * ref_fac)))
    numb_im_int = Int(round(sqrt(numb_roots * ref_fac)))

    # Mesh factor
    Delta_k_re = (Re_k_max - Re_k_min) / numb_re_int
    Delta_k_im = (Im_k_max - Im_k_min) / numb_im_int


    # Redefining regions
    if redef
        re_k_ext_or, im_k_ext_or = redef_part(Re_k_min, Re_k_max, Im_k_min, Im_k_max, numb_re_int, numb_im_int, max_iter, en_fac_mesh)
    else
        re_k_ext_or = Float64[Re_k_min + Delta_k_re * i for i in 0:numb_re_int]
        im_k_ext_or = Float64[Im_k_min + Delta_k_im * j for j in 0:numb_im_int]
    end
    re_k_ext = deepcopy(re_k_ext_or)
    im_k_ext = reverse(deepcopy(im_k_ext_or))
    # Integrands
    fun_der_int = k -> num_auto_diff(detLMat, k, h=new_r_mesh * 1e-3)
    integrand_fun_int(t, min_re, max_re, min_im, max_im) = 1/(2 * pi * im) * fun_der_int(gamma_rect(t, min_re, max_re, min_im, max_im)) / 
        detLMat(gamma_rect(t, min_re, max_re, min_im, max_im)) * gamma_rect_der(t, min_re, max_re, min_im, max_im)

    # Number of solutions per region    
    integrals = fill(0.0 + 0.0im, numb_re_int, numb_im_int)
    errors = fill(0.0, numb_re_int, numb_im_int)


    # Define the output file path
    file_path_ = joinpath(dirname(@__FILE__), "out_data.txt")

    # Define the message
    message = "A good parallelization in the function reg_with_few_sol is achieved when length(im_k_ext) = $(length(im_k_ext)), which is a multiple of the number of threads. If not, try changing ref_fac."

    # Write the message to the file
    open(file_path_, "w") do file
        write(file, message)
    end

    roots_count = 0

    if parall
        for j in 1:numb_im_int
            @floop for i in 1:numb_re_int
                integrals[i, j], errors[i, j] = quadgk(t -> integrand_fun_int(t, re_k_ext[i], re_k_ext[i+1], im_k_ext[j+1], im_k_ext[j]), 0, 4, atol=int_root_atol)
            end
            roots_count += sum(Int(round(real(integrals[i, j]))) for i in 1:numb_re_int)
            if roots_count > numb_sol_bound
                break
            end
        end
    else
        for j in 1:numb_im_int
            for i in 1:numb_re_int
                integrals[i, j], errors[i, j] = quadgk(t -> integrand_fun_int(t, re_k_ext[i], re_k_ext[i+1], im_k_ext[j+1], im_k_ext[j]), 0, 4, atol=int_root_atol)
            end
            roots_count += sum(Int(round(real(integrals[i, j]))) for i in 1:numb_re_int)
            if roots_count > numb_sol_bound
                break
            end
        end
    end



    # consistency test
    for i = 1:numb_re_int
        for j = 1:numb_im_int
            if errors[i, j] > int_root_check_tol || abs(imag(integrals[i, j])) > int_root_check_tol
                @warn "Something wrong with the integrals on ($i,$j)! Try to change ref_fac"
            end
        end
    end

    # # consistency test
    numb_sol_matrix = Int.(round.(real.(integrals)))
    # numb_sol_found = sum(numb_sol_matrix)
    # if numb_sol_found != numb_roots
    #     @warn "Something wrong with the number of solutions: numb sol found $numb_sol_found, numb roots $numb_roots. Maybe the region has changed!"
    # end

    # regions with roots
    max_numb_sol_matrix = maximum(numb_sol_matrix)
    ind_numb_sol = [[] for _ in 1:max_numb_sol_matrix]
    regions_with_roots = [[] for _ in 1:max_numb_sol_matrix]
    for i in 1:max_numb_sol_matrix
        ind_numb_sol[i] = findall(x -> x == i, numb_sol_matrix)
        for j in eachindex(ind_numb_sol[i])
            push!(regions_with_roots[i], [Re_k_min + (ind_numb_sol[i][j][1] - 1) * Delta_k_re, Re_k_min + ind_numb_sol[i][j][1] * Delta_k_re,
                Im_k_max - ind_numb_sol[i][j][2] * Delta_k_im, Im_k_max - (ind_numb_sol[i][j][2] - 1) * Delta_k_im]) # re_min, re_max, im_min, im_max
        end
    end

    fused_regions = vcat(regions_with_roots...)
    new_Im_k_min = minimum([fused_regions[i][3] for i in eachindex(fused_regions)])

    
    return regions_with_roots, new_Im_k_min #, max_numb_sol_matrix
end




#### SAVED DATA ###
#= the function partitions the region into approximately ref_fac*numb_roots =  subregions. 
    reg_with_few_sol[i] is of the from [[min_re_1, max_re_1, min_im_1, max_im_1], [min_re_2, max_re_2, min_im_2, max_im_2],...]
    where the arrays [min_re, max_re, min_im, max_im] indicate the regions with i roots =#
begin #=  7 min =#
    regions_with_roots, Im_k_min = reg_with_few_sol(Re_k_min, Re_k_max, Im_k_min, Im_k_max,  length(root_cleaned); max_iter=300, ref_fac=3, redef=false);
end
@save joinpath(dirname(@__FILE__), "out_regions_with_roots.jld2") regions_with_roots;
@save joinpath(dirname(@__FILE__), "out_Im_k_min.jld2") Im_k_min;
# @load joinpath(dirname(@__FILE__), "out_regions_with_roots.jld2") regions_with_roots
# @load joinpath(dirname(@__FILE__), "out_Im_k_min.jld2") Im_k_min
#### END SAVED DATA ####


# Initialize the plot with y-axis limits
p = plot(title="Regions with Roots", xlabel="Re(k)", ylabel="Im(k)", legend=true, ylim=(Im_k_min, Im_k_max));

# Generate colors
colors = palette(:viridis, length(regions_with_roots));

# Draw the filled rectangles
for (num_solutions, regions) in enumerate(regions_with_roots)
    for (region_idx, region) in enumerate(regions)
        re_min, re_max, im_min, im_max = region;

        # Coordinates of the rectangle (closed polygon)
        x = [re_min, re_max, re_max, re_min, re_min];
        y = [im_min, im_min, im_max, im_max, im_min];

        # Add the filled rectangle
        plot!(p, x, y, seriestype=:shape, color=colors[num_solutions], label=(region_idx == 1 ? "$(num_solutions) sol" : ""));
    end
end

# Display the plot
display(p);



regions_with_roots;
regions_with_roots[1];
regions_with_roots[2];
regions_with_roots[end];



numb_roots = sum(i*length(regions_with_roots[i]) for i in eachindex(regions_with_roots));

open(joinpath(dirname(@__FILE__), "out_data.txt"), "a") do file
    write(file, "\n");  # newline
    write(file, "New number of solutions: $numb_roots");
    write(file, "\n");
    #write(file, "Refined mesh: $new_r_mesh_s")
end




#### REFINEMENT BASED ON GRADIENT 


function plot_region(detLMat, region)

    re_min, re_max, im_min, im_max = region
    
    # Create the heatmap of the region
    re_vals_heatmap = range(re_min, re_max; length=100)
    im_vals_heatmap = range(im_min, im_max; length=100) 
    z_vals = [Complex(re, im) for re in re_vals_heatmap, im in im_vals_heatmap]
    abs_vals = abs.(detLMat.(z_vals))
                
    # Create the plot
    p = heatmap(re_vals_heatmap, im_vals_heatmap, abs_vals', title="Heatmap",
                xlabel="Re(z)", ylabel="Im(z)", color=:viridis)
    
    # Display the plot
    display(p)

end

plot_region(detLMat, regions_with_roots[1][1]);
#plot_region(detLMat, [11.29926,11.29929,-0.00387, -0.003847096551632433 ])
# it is better to refine mostly in the direction of real axis: large gradient along real(z)


# Function to estimate the average gradient along the real and imaginary axes
function ratio_average_gradient(f, region; samples=100, n_r_mesh=r_mesh)

    re_min, re_max, im_min, im_max = region
    r_mesh =  n_r_mesh#minimum([re_max-re_min, im_max-im_min])
    re_grads = []
    im_grads = []
    
    for _ in 1:samples
        z = Complex(rand(Uniform(re_min, re_max)), rand(Uniform(im_min, im_max)))
        h = r_mesh*1e-3
        re_grad = (abs(f(z+h))-abs(f(z-h)))/(2*h)
        im_grad = (abs(f(z+h*im))-abs(f(z-h*im)))/(2*h)
        re_abs = abs(re_grad)/(re_max-re_min)
        im_abs = abs(im_grad)/(im_max-im_min)
        if re_abs != NaN && im_abs != NaN
            push!(re_grads, re_abs)
            push!(im_grads, im_abs)
        end
    end
    
    avg_re_grad = mean(re_grads)
    avg_im_grad = mean(im_grads)
    ratio = avg_im_grad / avg_re_grad

    return ratio
end


# check
# @timed rat_a_g = ratio_average_gradient(detLMat, regions_with_roots[1][1]; samples=100, n_r_mesh=r_mesh)
# maximum(rat_a_g)
# minimum(rat_a_g)
# mean(rat_a_g)

function estimate_min_abs(region; samples=15, meth = 2)
    # Define the real and imaginary range for the rectangle
    re_min, re_max, im_min, im_max = region
    
    # Initialize the minimum value and corresponding k
    min_abs_val = Inf
    min_k = Complex(NaN, NaN)  # Store the k value corresponding to the minimum
    
    # Sample the function at random points within the rectangle
    for _ in 1:samples
        # Generate a random point within the rectangle
        re_k = rand(Uniform(re_min, re_max))
        im_k = rand(Uniform(im_min, im_max))
        k = Complex(re_k, im_k)
        
        if meth == 1
            # Calculate the absolute value of detLMat at k
            Mat = Matrix(LMat(k))
            abs_val = minimum(abs.(eigvals(Mat)))
        elseif meth == 2
            abs_val = abs(detLMat(k))
        end
        
        # Update the minimum if a smaller value is found
        if abs_val < min_abs_val
            min_abs_val = abs_val
            min_k = k
        end
    end
    
    # Return the estimated minimum value and corresponding k
    return min_abs_val#, min_k
end

# check
#@timed estimate_min_abs(regions_with_roots[1][1],samples=15)



function ref_reg_grad(init_region, numb_sol, int_root_atol; min_reg = 3, max_reg = 35, redef=false, n_mesh = r_mesh, debug = true)
    sub_regions = [[] for _ in 1:numb_sol]
    Re_k_min, Re_k_max, Im_k_min, Im_k_max = init_region

    # Calculate gradient ratio for the region
    ratio_av_grad = ratio_average_gradient(detLMat, init_region, samples=100)
    if ratio_av_grad == Inf
        ratio_av_grad = max_reg
    end
    if ratio_av_grad < 1
        # Set number of intervals based on the ratio
        numb_re_int = min_reg
        numb_im_int = minimum([Int(round(min_reg/ratio_av_grad)), max_reg])
    else
        numb_im_int = min_reg
        numb_re_int = minimum([Int(round(min_reg*ratio_av_grad)), max_reg])
    end

    if numb_re_int == 1 && numb_im_int == 1
        numb_re_int = 2
        numb_im_int = 2
    end

    # Calculate the mesh size for real and imaginary components
    Delta_k_re = (Re_k_max - Re_k_min) / numb_re_int
    Delta_k_im = (Im_k_max - Im_k_min) / numb_im_int
    new_r_mesh = r_mesh  # or choose minimum([Delta_k_re, Delta_k_im])

    # Redefine the grid if required
    if redef
        re_k_ext_or, im_k_ext_or = redef_part(Re_k_min, Re_k_max, Im_k_min, Im_k_max, numb_re_int, numb_im_int, max_iter, en_fac_mesh)
    else
        re_k_ext_or = Float64[Re_k_min + Delta_k_re * i for i in 0:numb_re_int]
        im_k_ext_or = Float64[Im_k_min + Delta_k_im * i for i in 0:numb_im_int]
    end
    re_k_ext = deepcopy(re_k_ext_or)
    im_k_ext = deepcopy(im_k_ext_or)

    # Define the regions based on the refined grid
    regions = []
    for (i, j) in product(1:numb_re_int, 1:numb_im_int)
        push!(regions, [re_k_ext[i], re_k_ext[i+1], im_k_ext[j], im_k_ext[j+1]])
    end

    # Estimate minimum absolute values in each region and sort regions by these values
    min_vals = estimate_min_abs.(regions; samples=15)
    sorted_indices = sortperm(min_vals)
    sorted_regions = [regions[sorted_indices[i]] for i in eachindex(sorted_indices)]
    
    # Define integrand functions for analysis
    fun_der_int = k -> num_auto_diff(detLMat, k, h=new_r_mesh*1e-3)
    integrand_fun_int(t, min_re, max_re, min_im, max_im) = 1/(2*pi*im) * fun_der_int(gamma_rect(t, min_re, max_re, min_im, max_im)) / 
        detLMat(gamma_rect(t, min_re, max_re, min_im, max_im)) * gamma_rect_der(t, min_re, max_re, min_im, max_im)

    # Find regions that contain the required number of solutions
    found_reg = 0
    for l in eachindex(sorted_regions)
        re_k_min, re_k_max, im_k_min, im_k_max = sorted_regions[l]
        # Compute the integral to determine the number of solutions in the sub-region
        integr, _ = quadgk(t -> integrand_fun_int(t, re_k_min, re_k_max, im_k_min, im_k_max), 0, 4, atol=int_root_atol)
        numb_sol_subreg = Int(round(real(integr)))
        
        # If the sub-region contains solutions, add it to the list
        if numb_sol_subreg > 0
            push!(sub_regions[numb_sol_subreg], sorted_regions[l])
            found_reg = found_reg + numb_sol_subreg
            if !debug && found_reg == numb_sol
                break
            end
        end
    end

    return sub_regions
end


# test to compile regions_with_roots[1][1] #time 3 sol, 4 min
sub_regions = ref_reg_grad(regions_with_roots[1][1], 1, int_root_atol; min_reg = 1, max_reg = 1, redef=false, n_mesh = r_mesh);


function recursive_refinement(regions, max_recursion; current_recursion=1, ref_fact=1, debug=false)
    
    
    # Array to hold the refined regions
    refined_regions_with_roots = [[] for _ in 1:length(regions)]
    # Copy regions with a single solution (no refinement needed)
    # refined_regions_with_roots[1] = deepcopy(regions[1])
    

    # Parallelized loop through regions that contain more than one solution
    for i in eachindex(regions)#2:length(regions)  # Start from i=2, refine only regions with multiple solutions
        @floop for j in eachindex(regions[i])
            init_region = regions[i][j]
            
            # Apply the refinement function to split the region
            refined_subregions = ref_reg_grad(init_region, i, int_root_atol; min_reg = 2, 
                max_reg = 45, redef=false, n_mesh = r_mesh, debug=debug)

            # Append the refined subregions to the corresponding arrays
            for k in eachindex(refined_subregions)
                append!(refined_regions_with_roots[k], refined_subregions[k])
            end
        end
    end
    
    # If max recursion depth is reached, return regions without further refinement
    if current_recursion > max_recursion || all(isempty, refined_regions_with_roots[2:length(refined_regions_with_roots)])
        return refined_regions_with_roots, current_recursion
    end
    
    
    # Recursively call the function on the refined regions
    return recursive_refinement(refined_regions_with_roots, max_recursion; current_recursion=current_recursion+1)
end



# ### SAVED DATA ###

# non deb 5 min
research_int_with_deg, n_recursion = recursive_refinement(regions_with_roots, 7; debug=false);
@save joinpath(dirname(@__FILE__), "out_research_int_with_deg.jld2") research_int_with_deg;
# @load joinpath(dirname(@__FILE__), "out_research_int_with_deg.jld2") research_int_with_deg
##### END SAVED DATA  ####


#################### check DEBUG
# # #  15 min deb
# # @timed research_int_with_deg_DB, n_recursion_DB = recursive_refinement(regions_with_roots, 7; debug=true)
# # @save joinpath(dirname(@__FILE__), "out_research_int_with_deg_DB.jld2") research_int_with_deg_DB
# @load joinpath(dirname(@__FILE__), "out_research_int_with_deg_DB.jld2") research_int_with_deg_DB
# isequal(Set(research_int_with_deg_DB[1]), Set(research_int_with_deg[1]))
########### END CHECK ##################



research_int_with_deg[1];
research_int_with_deg[2];
research_int_with_deg[3];
research_int_with_deg[4];
research_int_with_deg[end];



plot_region(detLMat, research_int_with_deg[1][1]);
plot_region(detLMat, research_int_with_deg[1][89]);
plot_region(detLMat, research_int_with_deg[1][182]);


# Function to compute and display the heatmap
function plot_min_eigvals_heatmap(region; n_points =35)
    # Extract the region boundaries
    Re_k_min, Re_k_max, Im_k_min, Im_k_max = region

    # Create grid of k points in the specified complex plane region
    re_values = range(Re_k_min, Re_k_max, length=n_points)
    im_values = range(Im_k_min, Im_k_max, length=n_points)
    heatmap_data = zeros(Float64, n_points, n_points)

    # Compute minimum(abs.(eigvals(LMat(k)))) for each point in the grid using @floop for parallelization
    @floop for i in 1:n_points
        re_p = re_values[i]
        for j in 1:n_points
            im_p = im_values[j]
            k = re_p + im_p * im  # Construct complex number k
            eig_values = eigvals(Matrix(LMat(k)))
            heatmap_data[i, j] = minimum(abs.(eig_values))
        end
    end

    # Display heatmap
    heatmap(re_values, im_values, heatmap_data, xlabel="Re(k)", ylabel="Im(k)", title="Minimum of abs(eigvals(LMat(k)))")
end

# # check
# #plot_min_eigvals_heatmap(research_int_with_deg[1][1])
# plot_min_eigvals_heatmap(research_int_with_deg[1][89])
# plot_min_eigvals_heatmap(research_int_with_deg[1][182])



function find_highest_degeneration(arrays)
    # Start from the end and go backwards
    for i in length(arrays):-1:1
        if !isempty(arrays[i])  # If the array is not empty
            return i  # Return the index
        end
    end
    return 0  # If all arrays are empty, return 0 or another appropriate value
end
new_max_deg = find_highest_degeneration(research_int_with_deg);

degeneracy_table = [0 for i = 1:new_max_deg];
for i = 1:new_max_deg
    degeneracy_table[i] = length(research_int_with_deg[i]);
end

degeneracy_table;

if length(degeneracy_table) > 1
    @warn "Possible degenerations! Try to change ref_fac.";
    open(joinpath(dirname(@__FILE__), "out_data.txt"), "a") do file
        write(file, "\n");  # newline
        write(file, "Degeneracy table: $degeneracy_table");
    end
end

research_int = vcat(research_int_with_deg...);

if length(research_int) != numb_roots
    @warn "Problem with the numbers of roots! Try to check degeneracies.";

    # Path to the output file
    file_path = joinpath(@__DIR__, "out_data.txt");

    # Create the file if it does not exist, and append the warning along with the code line
    open(file_path, "a") do io
        write(io, "Problem with the numbers of roots! Try to check degeneracies.\n");
        write(io, "Line: if length(research_int) != numb_roots\n");
    end
end



#### ATTEMPT TO FIND ROOTS: NEWTON - RAPHSON ####

### NOTE: avoiding newton_raphson exiting the region does not work!


function filter_points_in_region(path, region)
    re_min, re_max, im_min, im_max = region
    min_candidates = [z for z in path if re_min <= real(z) <= re_max && im_min <= imag(z) <= im_max]
    return min_candidates
end

# memorize the paths
function newton_raphson_path(detLMat, region; max_iter_nr=100, grid_size=20, tolerance=1e-10, show_plot=true)
    re_min, re_max, im_min, im_max = region
    r_mesh = minimum([re_max - re_min, im_max - im_min])
    
    # Generate a grid of initial points within the region
    re_vals = LinRange(re_min, re_max, grid_size)
    im_vals = LinRange(im_min, im_max, grid_size)
    initial_points = [Complex(re, im) for re in re_vals, im in im_vals]
    #max_point_iter = grid_size^2

    # Define the function for the derivative
    f_prime(z) = num_auto_diff(detLMat, z, h=r_mesh * 10.0^(-3))
    
    all_paths = []  # Array to store all Newton-Raphson paths
    min_values = []  # Initialize the minimum value to infinity
    #min_candidates = [] # values in the region

    for z0 in initial_points
        # Initialize the trace for the Newton-Raphson path
        path = [z0]
        z = z0
        nr_iter = 0  # Track the number of Newton-Raphson iterations
        for i in 1:max_iter_nr
            fz = detLMat(z)
            Matz = Matrix(LMat(z))
            if any(x -> isinf(x) || isnan(x), Matz) || isnan(fz) || isinf(fz)
                break
            end
            eig_z = minimum(abs.(eigvals(Matz)))
            fz_prime = f_prime(z)
            if fz_prime == 0
                break
            end
            z = z - fz / fz_prime
            push!(path, z)
            nr_iter = i
            if eig_z < tolerance && re_min <= real(z) <= re_max && im_min <= imag(z) <= im_max #abs(fz) < tolerance 
                break  # Convergence when the value is sufficiently small
            end
        end

        # Store the path
        push!(all_paths, path)
        # store the minimum in the region
        min_candidates = filter_points_in_region(path, region)
        arg_cand = argmin([minimum(abs.(eigvals(Matrix(LMat(z))))) for z in min_candidates])
        push!(min_values, min_candidates[arg_cand])
        #push!(min_values,min_candidates[argmin(abs.(detLMat.(min_candidates)))])
        if minimum(abs.(eigvals(Matrix(LMat(min_candidates[arg_cand]))))) < tolerance
            break
        end
    end
    
    arg_best_path = argmin([minimum(abs.(eigvals(Matrix(LMat(z))))) for z in min_values])
    #arg_best_path = argmin(abs.(detLMat.(min_values)))

    if show_plot && length(all_paths[arg_best_path]) > 0
        # Create the heatmap of the region
        re_vals = LinRange(re_min, re_max, 300)
        im_vals = LinRange(im_min, im_max, 300)
        z_vals = [Complex(re, im) for re in re_vals, im in im_vals]
        abs_vals = abs.(detLMat.(z_vals))
        
        # Create the plot
        p = heatmap(re_vals, im_vals, abs_vals', title="Newton-Raphson Path + Heatmap",
                    xlabel="Re(z)", ylabel="Im(z)", color=:viridis)

        # Add the best Newton-Raphson path
        path_re = real.(all_paths[arg_best_path])
        path_im = imag.(all_paths[arg_best_path])
        plot!(path_re, path_im, label="Best Path", marker=:circle, linecolor=:red, linewidth=2)
        
        # Add the initial point in red
        scatter!(path_re[1:1], path_im[1:1], color=:red, label="Initial Point", markersize=6)
        
        # Add the intermediate points in gray
        if length(all_paths[arg_best_path]) > 2
            scatter!(path_re[2:end-1], path_im[2:end-1], color=:gray, label="Intermediate Points", markersize=4)
        end
        
        # Add the final point in green
        scatter!(path_re[end:end], path_im[end:end], color=:green, label="Final Point", markersize=6)

        # Display the plot
        display(p)
    end
    
    return min_values[arg_best_path]
end



# to compile 
newton_raphson_path(detLMat, research_int[1]; max_iter_nr=1, grid_size=2, tolerance=1e-0, show_plot=false);

### Checks
# @timed try88 = newton_raphson_path(detLMat, research_int[88]; max_iter_nr=50, grid_size=5, tolerance=1e-10, show_plot=false)
# minimum(abs.(eigvals(Matrix(LMat(try88)))))
# @timed try1 = newton_raphson_path(detLMat, research_int[1]; max_iter_nr=50, grid_size=5, tolerance=1e-10, show_plot=false)
# minimum(abs.(eigvals(Matrix(LMat(try1)))))


function roots_attempt_nr_func(research_int)
    roots_attempt_nr = [0.0+im*0.0 for i in eachindex(research_int)]

    @floop for i in eachindex(research_int)
        roots_attempt_nr[i] = newton_raphson_path(detLMat, research_int[i]; max_iter_nr=50, grid_size=5, tolerance=1e-10, show_plot=false)
    end

    return roots_attempt_nr
end

##### SAVED DATA ########
# 7 min
roots_attempt_nr = roots_attempt_nr_func(research_int);
@save joinpath(dirname(@__FILE__), "out_roots_attempt_nr.jld2") roots_attempt_nr;
# @load joinpath(dirname(@__FILE__), "out_roots_attempt_nr.jld2") roots_attempt_nr
#### END SAVED DATA ###


roots_attempt_nr_abs=[minimum(abs.(eigvals(Matrix(LMat(k))))) for k in roots_attempt_nr]; #abs.(detLMat.(roots_attempt_nr))
sort(roots_attempt_nr_abs);
histogram(roots_attempt_nr_abs, bins=100, xlabel="Values", ylabel="Freq.");
# Create the file if it does not exist, and append the maximum value
# Path to the output file
file_path_datatxt = joinpath(@__DIR__, "out_data.txt");

open(file_path_datatxt, "a") do io
    write(io, "Maximum of roots_attempt_nr_abs: $(maximum(roots_attempt_nr_abs))\n");
end

# Function to find close pairs of complex numbers
function find_close_pairs(kroots, min_dist)
    close_pairs = []  # List to store pairs of indices that are too close

    # Brute force approach: check the distance between each pair of complex numbers
    for i in 1:length(kroots)
        for j in i+1:length(kroots)
            # Calculate the Euclidean distance between the two complex numbers
            if norm(kroots[i] - kroots[j]) < min_dist
                # If the distance is less than the minimum threshold, store the pair
                push!(close_pairs, (i, j))
            end
        end
    end

    return close_pairs  # Return the list of close pairs
end


# Call the function to find close pairs
close_pairs = find_close_pairs(roots_attempt_nr, r_mesh);

close_pairs;

#### kroots #######


kroots = roots_attempt_nr;


function center_regions(research_int, kroots)
    centered_research_int = []

    for (region, root) in zip(research_int, kroots)
        re_min, re_max, im_min, im_max = region
        re_center, im_center = real(root), imag(root)
        
        # Calculate current center of the region
        current_re_center = (re_min + re_max) / 2
        current_im_center = (im_min + im_max) / 2

        # Calculate the shifts needed to move the root to the center
        re_shift = re_center - current_re_center
        im_shift = im_center - current_im_center

        # Shift the region bounds
        new_re_min = re_min + re_shift
        new_re_max = re_max + re_shift
        new_im_min = im_min + im_shift
        new_im_max = im_max + im_shift
        
        push!(centered_research_int, [new_re_min, new_re_max, new_im_min, new_im_max])
    end

    return centered_research_int
end

# Example usage
centered_research_int = center_regions(research_int, kroots);


function ordered_indices_by_proximity(kroots)
    n = length(kroots)
    # Calculate the distance matrix for all pairs
    dist_matrix = [abs(kroots[i] - kroots[j]) for i in 1:n, j in 1:n]
    inds = Vector{Int}(undef, n)
    
    # Start with the pair having the minimum distance
    _, min_idx = findmin(dist_matrix + I * maximum(dist_matrix)) # Avoid zero-distance with self by adding a large I
    i, j = Tuple(CartesianIndex(min_idx))
    inds[1] = i
    inds[2] = j
    used = Set([i, j])
    
    for k in 3:n
        # Find the next closest element to the last one added
        last_index = inds[k - 1]
        candidates = [(idx, dist_matrix[last_index, idx]) for idx in setdiff(1:n, used)]
        
        # Find the candidate with the minimum distance
        _, closest = findmin(candidates)
        next_index = candidates[closest][1]
        
        inds[k] = next_index
        push!(used, next_index)
    end
    
    return inds
end

# Proximity ordering
ordered_indices = ordered_indices_by_proximity(kroots);


plot_min_eigvals_heatmap(centered_research_int[ordered_indices[1]]);
plot_min_eigvals_heatmap(centered_research_int[ordered_indices[2]]);

function smallest_bounding_region(region1, region2)
    # Extract the bounds for each region
    re_min1, re_max1, im_min1, im_max1 = region1
    re_min2, re_max2, im_min2, im_max2 = region2
    
    # Calculate the bounds of the smallest bounding region
    re_min_bounding = min(re_min1, re_min2)
    re_max_bounding = max(re_max1, re_max2)
    im_min_bounding = min(im_min1, im_min2)
    im_max_bounding = max(im_max1, im_max2)
    
    # Return the smallest bounding region as a new array
    return [re_min_bounding, re_max_bounding, im_min_bounding, im_max_bounding]
end

#
bounding_region = smallest_bounding_region(centered_research_int[ordered_indices[1]], centered_research_int[ordered_indices[2]]);
plot_min_eigvals_heatmap(bounding_region);


function max_min_abs_eigvals_segment_with_point(k1, k2, LMat, n_points=100)
    # Define points k along the segment between k1 and k2
    k_values = [k1 + t * (k2 - k1) for t in range(0, stop=1, length=n_points)]
    
    # Compute the minimum of the absolute eigenvalues for each sampled k
    min_abs_eigvals = [minimum(abs.(eigvals(Matrix(LMat(k))))) for k in k_values]
    
    # Find the maximum of these minimum values and its corresponding point
    max_min_value = maximum(min_abs_eigvals)
    max_min_index = argmax(min_abs_eigvals)  # Index of the maximum
    max_min_point = k_values[max_min_index]  # Corresponding k point

    # Return the maximum value and its corresponding point
    return max_min_value, max_min_point
end

# Usage of the function
k1 = kroots[ordered_indices[1]];
k2 = kroots[ordered_indices[2]];
result, point = max_min_abs_eigvals_segment_with_point(k1, k2, LMat);

######### EIGENVALUES  ######


# checks: two smallest eigenvalues for each k in kroot
# LMat_array = LMat.(kroots)
# @timed eigs(LMat_array[100], nev=2, which=:SM, ncv=40)
# @timed eigvals(Array(LMat_array[100]))
# E = eigen(Array(LMat_array[100]))
# amin = argmin(abs.(E.values))
# E.vectors[:,amin]
# norm(LMat(kroots[100])*E.vectors[:,amin])
# typeof(E)
# size(LMat_array[1])[1]

# return the two smallest eigenvalues and the eigenvector corresponding to the smallest eigenvalue
function smallest_eigen(kroots; meth=1)
    if meth == 1
        LMat_array = Array.(LMat.(kroots))
        E = Vector{LinearAlgebra.Eigen}(undef, length(kroots))
        eigenvalues = Vector{Vector{ComplexF64}}(undef, length(kroots))
        sorted_eigenvalues = Vector{Vector{ComplexF64}}(undef, length(kroots))
        eigenvectors = Vector{Matrix{ComplexF64}}(undef, length(kroots))
        pairs_eigenvalues = [Vector{ComplexF64}(undef, 2) for _ in eachindex(kroots)]
        first_eigenvectors = Vector{Vector{ComplexF64}}(undef, length(kroots))
        
        @threads for i in eachindex(kroots)
            E[i] = eigen(LMat_array[i])
            eigenvalues[i] = E[i].values
            eigenvectors[i] = E[i].vectors
            sorted_eigenvalues[i] = sort(eigenvalues[i], by=abs)
            pairs_eigenvalues[i] = sorted_eigenvalues[i][1:2]
            first_eigenvectors[i] = eigenvectors[i][:, argmin(abs.(eigenvalues[i]))]
        end
    elseif meth == 2
        LMat_array = LMat.(kroots)
        results = Vector{Any}(undef, length(kroots))
        
        @threads for i in eachindex(kroots)
            results[i] = eigs(LMat_array[i], nev=2, which=:SM, ncv=40)
        end
        
        pairs_eigenvalues = [result[1] for result in results]
        first_eigenvectors = [result[2][:, 1] for result in results]
    end
    
    return pairs_eigenvalues, first_eigenvectors
end

########### SAVED DATA #####
begin # ~14 sec
sm_eigvls, sm_eigvec = smallest_eigen(kroots, meth=1);
end
@save joinpath(dirname(@__FILE__), "out_sm_eigvls.jld2") sm_eigvls;
@save joinpath(dirname(@__FILE__), "out_sm_eigvec.jld2") sm_eigvec;
# @load joinpath(dirname(@__FILE__), "out_sm_eigvls.jld2") sm_eigvls
# @load joinpath(dirname(@__FILE__), "out_sm_eigvec.jld2") sm_eigvec
########## END SAVED DATA


# checks
# sm_eigvls[100]==sort(E.values, by=abs)[1:2]
# sm_eigvec[100]==E.vectors[:,amin]

function check_degeneracy_debug(eigvls; sef=1000)
    # We set the second smallest eigenvalue divided by 1000 as a small value with respect to the matrix elements
    small_value_ref = [abs(eigvls[i][2])/sef  for i = 1:length(kroots)] ;


    # If small_value_ref is smaller than the smallest eigenvalue, it raises an ErrorException
    for i in eachindex(kroots)
        if small_value_ref[i] < abs(eigvls[i][1])
            error("Probable degeneracy problem! The program will be stopped! \U1F4A5")
        end
    end

    print("No degeneracy problem! \U1F44D")
end

check_degeneracy_debug(sm_eigvls; sef=2);

# It gives the largest norm of the vectors L(k)*eigevector(k), for k in kroots (it must be very small!)
function check_eigensys_debug(kroots,sm_eigvec)
    small_value = 0.0

    for i in eachindex(kroots)
        small_vec_norm = norm(LMat(kroots[i]) * sm_eigvec[i])
        if small_vec_norm > small_value
            small_value = small_vec_norm
        end
    end

    return small_value
end

check_eigensys_debug(kroots,sm_eigvec);

# check if there are proportional vectors

# compare two vectors
function are_proportional(arr1, arr2; tolerance=1e-2)
    if length(arr1) != length(arr2)
        return false
    end

    ratios = []
    for i in 1:length(arr1)
        if arr2[i] != 0
            push!(ratios, arr1[i] / arr2[i])
        elseif arr2[i] == 0 && abs(arr2[i])>tolerance
            return false
        end
    end

    ratio_mean = mean(ratios)

    for r in ratios
        if abs(r - ratio_mean) > tolerance
            return false
        end
    end

    return true
end

# compare all vectors
function check_proport(sm_eigvec)
    # Length of the array
    n = length(sm_eigvec)

    # Initialize an empty array to store the pairs (i, j); i<j
    indices = []
    # Iterate through all possible pairs (i, j) where i < j
    for i in 1:(n-1)  # Loop i from 1 to n-1
        for j in (i+1):n  # Loop j from i+1 to n
            push!(indices, [i, j])
        end
    end

    # Initialize an empty array to store the pairs of proportional vectors indices
    prop_vec_indices = []
    @threads for i in eachindex(indices)
        if are_proportional(sm_eigvec[indices[i][1]],sm_eigvec[indices[i][2]])
            push!(prop_vec_indices, indices[i])
        end
    end

    return prop_vec_indices
end

prop_vectors_ind = check_proport(sm_eigvec);
@save joinpath(dirname(@__FILE__), "out_prop_vectors_ind.jld2") prop_vectors_ind;
if !isempty(prop_vectors_ind)
    @warn "Proportional pairs detected: check if they correspond to close eigenvalues!";

    # Path to the output file
    file_path = joinpath(@__DIR__, "out_data.txt");

    # Write the warning and the content of prop_vectors_ind to the file
    open(file_path, "a") do io
        write(io, "Proportional pairs detected: check if they correspond to close eigenvalues!\n");
        write(io, "Line: if !isempty(prop_vectors_ind)\n");
        write(io, "Proportional vectors indices: $(prop_vectors_ind)\n");
    end
end


# Electric field

# Define the eta_mu function (passive modes)
function eta_mu(x, kroot, sm_eigvec, nij, fibersLength, N_int, N_ext)
  
    k = kroot
    # Initialize the result
    result = zeros(Complex{Float64}, N_int + N_ext)
    
    # Calculate the result for i = 1,...,N_int
    for i in 1:N_int
        lambda_plus = sm_eigvec[i]
        lambda_minus = sm_eigvec[N_int + N_ext + i]
        n_i = nij[i]
        l_i = fibersLength[i]
        result[i] = lambda_plus * exp(im * k * n_i * x) + lambda_minus * exp(im * k * n_i * (l_i - x))
    end
    
    # Calculate the result for i = N_int+1,...,N_int+N_ext
    for i in N_int+1:N_int+N_ext
        lambda_plus = sm_eigvec[i]
        n_i = nij[i]
        result[i] = lambda_plus * exp(im * k * n_i * x)
    end
    
    return result
end



####### GAIN MEDIUM  ####



function nij_g(k,D0)
    n = length(nij)
    nij_g = [0.0 + 0.0im for i in 1:n]
    for i in eachindex(nij)
        nij_g[i] = sqrt(nij[i]^2 + delta_pump[i] * D0*gamma_perp/(k - k_a + im*gamma_perp))
    end
    return nij_g
end
nij_g(2.7+im,0.0)


###### Matrix with gain medium #####
function LMatGain_num(k,D0)

    
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
    LMat = spzeros(Complex{Float64}, NE + NEI, NE + NEI)
    #LMat = spzeros(Complex{typeof(k)}, NE + NEI, NE + NEI)
     
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


# Include the generated function file
include(joinpath(dirname(@__FILE__), "LMatGain_dense.jl"));


###### CHECKS #####
# LMatGain_num(2.1-0.01im, 0.3)
# LMatGain_dense(2.1-0.01im, 0.3)
# @timed LMatGain_num(2.1-0.01im, 0.3)
# @timed LMatGain_dense(2.1-0.01im, 0.3)

# maximum(abs.(LMatGain_num(2.1-0.01im, 0.3) .- LMatGain_dense(2.1-0.01im, 0.3)))
# maximum(abs.(LMatGain_num(2.1-0.01im, 0.0) .- LMat(2.1-0.01im)))

# maximum(imag.(kroots))
# minimum(imag.(kroots))





LMatGain = LMatGain_dense;



# determinant
function detLMatGain(k,D0)
    Mat = LMatGain(k,D0);

    return det(Mat)
end

# determinant
function detLMatGain_num(k,D0)
    Mat = LMatGain_num(k,D0);

    return det(Mat)
end


############ CHECKS
# Initialize variables to store the maximum difference
# max_diff = 0.0
# max_k = 0 + 0im  # The k value that maximizes the difference
# # Loop over the grid of k (complex numbers) and D0
# for x in x_vals
#     for y in y_vals
#         k = Complex(x, y)  # Construct the complex k = x + yi
#             # Calculate the absolute difference for the current (k, D0)
#             diff = abs(detLMatGain_num(k, 0.0) - detLMat_num(k))
#             if diff > max_diff
#                 # Update the maximum difference and store the corresponding (k, D0)
#                 max_diff = diff
#                 max_k = k
#             end
#     end
# end

# Output the maximum difference and the corresponding (k, D0)
#println("Maximum difference: $max_diff found at k = $max_k")
################# END CHECKS


# checking that δD0 ~ 1/100 is a good choice (since max| n_NI_ref^2 | < 1 and n1 >= 1.0)
n_NI_ref_sqr(k) = gamma_perp/(k-k_a-im*gamma_perp);
if maximum(abs.(n_NI_ref_sqr.(kroots))) > 1
    @warn "n_NI_ref^2 is greater than 1"
end

# set δD0
l_max = maximum(fibersLength);
k_max = maximum(abs.(kroots));
sf_fac_gain = 15; # safety factor
delta_D0 = 1/(l_max*k_max)^2/sf_fac_gain;
max_pump_est = 1;  # max pump estimated D0
iter_D0 = round(Int, max_pump_est/delta_D0);  # max iteration est

######  GAIN MEDIUM: TLMs ######


function gamma_mu(k_mu, k_a, gamma_perp)
    gamma = gamma_perp/(k_mu - k_a + im*gamma_perp)
    return gamma
end
deltak = abs(maximum(imag.(kroots)))*10^(-4);
maximum(imag.(kroots));


# integrand(x) = eta_mu(x, kroots[1], sm_eigvec[1], nij, fibersLength, NEI, NEO)[4]^2
# quadgk(integrand, 0, fibersLength[4], atol=deltak)
all(x -> x == nij[1], nij) && all(x -> x == 1.0, delta_pump);

function f_mu(k_mu, sm_eigvec_mu, nij, fibersLength, NEI, NEO, deltak; delta_pump=delta_pump)
    
    if all(x -> x == nij[1], nij) && all(x -> x == 1.0, delta_pump)
        f=1/nij[1]^2
    else
        integrals1 = zeros(Complex{Float64}, NEO + NEI)
        errors1 = zeros(Float64, NEO + NEI)
        @floop for i in 1:NEI+NEO
            integrals1[i], errors1[i] = quadgk(x-> delta_pump[i] *
                eta_mu(x, k_mu, sm_eigvec_mu, nij, fibersLength, NEI, NEO)[i]^2,
                0, fibersLength[i], atol=deltak)
        end

        if maximum(errors1) > deltak*100
            @warn "Problem with eta_ij integration!"
        end

        integrals2 = zeros(Complex{Float64}, NEO + NEI)
        errors2 = zeros(Float64, NEO + NEI)
        @floop for i in 1:NEI+NEO
            integrals2[i], errors2[i] = quadgk(
                x -> eta_mu(x, k_mu, sm_eigvec_mu, nij, fibersLength, NEI, NEO)[i]^2*nij[i]^2,
                0, fibersLength[i], atol=deltak)
        end

        if maximum(errors2) > deltak*100
            @warn "Problem with eta_ij integration!"
        end

        f = sum(integrals1)/sum(integrals2)
    end
    
    return f

end

f_values = zeros(Complex{Float64}, length(kroots))
##### SAVED DATA
begin  
@floop for i in eachindex(kroots)
    f_values[i] = f_mu(kroots[i], sm_eigvec[i], nij, fibersLength, NEI, NEO, deltak);
end
end
@save joinpath(dirname(@__FILE__), "out_f_values.jld2") f_values;
# @load joinpath(dirname(@__FILE__), "out_f_values.jld2") f_values
####### END SAVED DATA




# Newton Raphson (with memorization of the path)
function newton_raphson_gain(detLMatGain, k_D0_in_mu, D0_in, delta_D0, sm_eigvec_mu, nij, fibersLength, NEI, NEO, f_values_mu; max_iter_nr=100, tolerance=1e-9, show_plot=false, meth=1)
    # typical length scale
    r_mesh = abs(imag(k_D0_in_mu)) 
    
    # define the function
    D0 = D0_in + delta_D0
    detLMat0(z) = detLMatGain(z, D0)
    #######
    LMatGain0(z) = LMatGain(z, D0)
    ########
    # define gamma_mu
    gamma_mu_val = gamma_mu(k_D0_in_mu, k_a, gamma_perp)

    # Define the derivative
    f_prime(z) = num_auto_diff(detLMat0, z, h=r_mesh * 10.0^(-3))
    
    if meth == 2
        z0 = k_D0_in_mu
    elseif meth == 1
        z0 = k_D0_in_mu*sqrt((1 + D0_in*gamma_mu_val*f_values_mu)/(1+D0*gamma_mu_val*f_values_mu))
    end
    
    # Initialize the trace for the Newton-Raphson path
    path = [z0]
    z = z0
    fz = detLMat0(z)
    ###########
    LMatz = LMatGain0(z)
    eigenvalz = sort(eigvals(LMatz), by=abs)[1]
    ##########
    fz_prime = f_prime(z)
    nr_iter = 0  # Track the number of Newton-Raphson iterations
    while nr_iter < max_iter_nr && abs(eigenvalz) > tolerance #abs(fz) > tolerance 
        while fz_prime == 0 && nr_iter < max_iter_nr
            nr_iter = nr_iter + 1
            z = z + rand(Uniform(-r_mesh*10.0^(-6)/2, r_mesh*10.0^(-6)/2)) +
                im*rand(Uniform(-r_mesh*10.0^(-6)/2, r_mesh*10.0^(-6)/2))
            fz = detLMat0(z)
            ##############
            LMatz = LMatGain0(z)
            eigenvalz = sort(eigvals(LMatz), by=abs)[1]
            ############
            fz_prime = f_prime(z)
        end
        z = z - fz / fz_prime
        push!(path, z)
        fz = detLMat0(z)
        ##############
        LMatz = LMatGain0(z)
        eigenvalz = sort(eigvals(LMatz), by=abs)[1]
        ############
        fz_prime = f_prime(z)
        nr_iter = nr_iter + 1
    end
    
    # define plotting region
    real_path = real.(path)
    imag_path = imag.(path)
    re_k_D0_in_mu = real(k_D0_in_mu)
    im_k_D0_in_mu = imag(k_D0_in_mu)
    re_min = minimum(vcat(real_path, re_k_D0_in_mu))
    re_max = maximum(vcat(real_path, re_k_D0_in_mu))
    im_min = minimum(vcat(real_path, im_k_D0_in_mu))
    im_max = maximum(vcat(imag_path, [0.0, im_k_D0_in_mu]))
    delta_re = re_max - re_min
    delta_im = im_max - im_min
    re_min = re_min - delta_re/50.0
    re_max = re_max + delta_re/50.0
    im_min = im_min - delta_im/50.0
    im_max = im_max + delta_im/50.0

    # best value for the root
    reverse_path = reverse(path)
    min_index_abs = argmin(abs.(detLMat0.(reverse_path)))

    if show_plot
        # Create the heatmap of the region
        re_vals = LinRange(re_min, re_max, 200)
        im_vals = LinRange(im_min, im_max, 200)
        z_vals = [Complex(re, im) for re in re_vals, im in im_vals]
        abs_vals = abs.(detLMat0.(z_vals))
        
        # Create the plot
        p = heatmap(re_vals, im_vals, abs_vals', title="Newton-Raphson Path + Heatmap",
                    xlabel="Re(z)", ylabel="Im(z)", color=:viridis)

        # Add the best Newton-Raphson path
        plot!(real_path, imag_path, label="NR Path", marker=:circle, linecolor=:red, linewidth=2)
        
        # Add the initial point in red
        scatter!(real_path[1:1], imag_path[1:1], color=:red, label="Initial Point", markersize=6)
        
        # Add the intermediate points in gray
        if length(path) > 2
            scatter!(real_path[2:end-1], imag_path[2:end-1], color=:gray, label="Intermediate Points", markersize=4)
        end
        
        # Add the final point in green
        scatter!(real_path[end:end], imag_path[end:end], color=:green, label="Final Point", markersize=6)

        # Add starting point
        scatter!([re_k_D0_in_mu], [im_k_D0_in_mu], color=:red, label="Starting Point", markersize=6, marker=:star)

        # Display the plot
        display(p)
        return reverse_path[min_index_abs], nr_iter
    end
    
    return reverse_path[min_index_abs], nr_iter
end

## trying NR
try1, _ = newton_raphson_gain(detLMatGain, kroots[1], 0, delta_D0, sm_eigvec[1], nij, fibersLength, NEI, NEO, f_values[1]; max_iter_nr=1000, tolerance=1e-9, show_plot=false, meth=2);
try1;
sort(eigvals(LMatGain(try1, delta_D0)), by=abs)[1];

# for the last point Im k ~ 0 is required
function NR_last_point_gain(k, D0, detLMatGain, delta_k, delta_D0; tol=1e-9, max_iter=1000)
    
    iter_indx = 0

    while iter_indx < max_iter

        iter_indx += 1
        # Current values of a and b
        a, b = real(k), imag(k)

        # Compute the value of (M_g)_{2,2}(k, D0)
        M_g = detLMatGain(k, D0)
        M_g_real = real(M_g)
        M_g_imag = imag(M_g)

        # Compute partial derivatives
        ∂Mg_∂a = num_auto_diff((x) -> detLMatGain(x + im * b, D0), a, h=delta_k * 10.0^(-3))
        ∂Mg_∂b = num_auto_diff((x) -> detLMatGain(a + im * x, D0), b, h=delta_k * 10.0^(-3))
        ∂Mg_∂D0 = num_auto_diff((x) -> detLMatGain(k, x), D0, h=delta_D0 * 10.0^(-3))

        ∂Mg_∂a_real = real(∂Mg_∂a)
        ∂Mg_∂b_real = real(∂Mg_∂b)
        ∂Mg_∂D0_real = real(∂Mg_∂D0)

        ∂Mg_∂a_im = imag(∂Mg_∂a)
        ∂Mg_∂b_im = imag(∂Mg_∂b)
        ∂Mg_∂D0_im = imag(∂Mg_∂D0)

        # Construct the Jacobian matrix
        J = [
            ∂Mg_∂a_real   ∂Mg_∂b_real   ∂Mg_∂D0_real
            ∂Mg_∂a_im     ∂Mg_∂b_im     ∂Mg_∂D0_im
            0.0           1.0           0.0
        ]

        # Right-hand side vector
        F = [M_g_real, M_g_imag, b]

        if det(J) == 0
            # If matrix is singular, adjust k and D0 slightly
            k += rand() * delta_k / 100.0
            D0 += rand() * delta_D0 / 100.0
            continue
        end
       
        # Solve the linear system
        Δ = -J \ F
        Δa, Δb, ΔD0 = Δ

        # Update k and D0
        k = (a + Δa) + im * (b + Δb)
        D0 += ΔD0

        # Check for convergence
        if (b + Δb) < tol && abs(sort(eigvals(LMatGain(k, D0)), by=abs)[1]) < tol #abs(detLMatGain(k, D0)) < tol
            break
        end
    
    end

    return k, D0, iter_indx, abs(detLMatGain(k, D0))
end

# Function to calculate the mean of the last two elements of an array of arrays
function mean_last_two_elements(arrays)
    means = []
    for array in arrays
        if length(array) >= 2
            last_two = array[end-1:end]
            mean_value = sum(last_two) / length(last_two)
            push!(means, mean_value)
        elseif length(array) == 1
            push!(array[1], mean_value)
        else
            error("Problem with NR for TLM: each array must have at least one element")
        end
    end
    return means
end

# Function to compute the minimum difference between the last two elements of an array of arrays
function min_last_two_diff(arrays)
    min_diff = Inf
    for array in arrays
        if length(array) >= 2
            last_diff = abs(array[end] - array[end-1])
            if last_diff < min_diff
                min_diff = last_diff
            end
        else
            error("Problem with NR for TLM: each array must have at least 2 elements")
        end
    end
    return min_diff
end


# Function path_threshold, by increasing D_0
function path_threshold(detLMatGain, kroots, delta_D0, sm_eigvec, nij, fibersLength, NEI, NEO, f_values; max_iter_pt=10^23, max_D0 = false, max_D0_value=100)
    all_paths = [[] for i in eachindex(kroots)]
    all_D0s_paths = [[] for i in eachindex(kroots)]
    all_D0s = [0.0 for i in eachindex(kroots)]
    
    @floop for l in eachindex(kroots)
        z0 = kroots[l]
        path = [z0]
        path_D0 = [0.0]
        zn = z0
        D0_init = 0
        sm_eigvec_in = sm_eigvec[l]
        f_mu = f_values[l]
        for _ in 1:max_iter_pt
            zn, _ = newton_raphson_gain(detLMatGain, zn, D0_init, delta_D0, sm_eigvec_in, nij, fibersLength, NEI, NEO, f_mu; max_iter_nr=25, tolerance=1e-5, show_plot=false, meth=1)#(detLMatGain, zn, D0; max_iter_nr=100, tolerance=1e-3, show_plot=false)
            push!(path, zn)
            if imag(zn) > 0 || (max_D0 && D0 > max_D0_value)
                push!(path_D0, D0_init)
                break
            end
            D0_init = D0_init + delta_D0  # Add delta_D0 at each iteration
            push!(path_D0, D0_init)
        end
        all_paths[l] = path
        all_D0s_paths[l] = path_D0
        all_D0s[l] = path_D0[end]
    end
    
    # Calculate the mean of the last two elements for each array in all_paths_th_0
    means_all_paths_th_0 = mean_last_two_elements(all_paths)
    # Calculate the mean of the last two elements for each array in all_D0s_paths
    means_all_D0s_paths = mean_last_two_elements(all_D0s_paths)

    # Calculate mesh_k as the minimum difference of the last two elements for each array in all_paths_th_0
    mesh_k = min_last_two_diff(all_paths)
    # Calculate mesh_D0 as the minimum difference of the last two elements for each array in all_D0s_paths
    mesh_D0 = min_last_two_diff(all_D0s_paths)
    
    success_array = [imag(all_D0s_paths[l][end]) >= 0 for l in eachindex(all_D0s_paths)]
    valid_indices = collect(filter(l -> success_array[l], eachindex(all_paths)))

    @floop for l in valid_indices
        all_paths[l][end], all_D0s_paths[l][end], _, _ = NR_last_point_gain(
            all_paths[l][end], 
            all_D0s_paths[l][end], 
            detLMatGain, 
            mesh_k, 
            mesh_D0; 
            tol=1e-9, max_iter=1000
        )
        all_D0s[l] = all_D0s_paths[l][end]
    end

    return all_paths, all_D0s, all_D0s_paths, success_array
end



# Function path_threshold, by increasing D_0
function path_threshold_bar(detLMatGain, kroots, delta_D0, sm_eigvec, nij, fibersLength, NEI, NEO, f_values; max_D0 = false, max_D0_value=100)
    all_paths = [[] for i in eachindex(kroots)]
    all_D0s_paths = [[] for i in eachindex(kroots)]
    all_D0s = [0.0 for i in eachindex(kroots)]
    
    # # Create a progress bar
    # progress = Progress(length(kroots), 1)  # Length equal to the number of kroots, increment of 1
    
    @floop for l in eachindex(kroots)
        z0 = kroots[l]
        path = [z0]
        path_D0 = [0.0]
        zn = z0
        D0_init = 0.0
        sm_eigvec_in = sm_eigvec[l]
        f_mu = f_values[l]
        while imag(zn) < 0 && (!max_D0 || D0_init < max_D0_value)
            zn, _ = newton_raphson_gain(detLMatGain, zn, D0_init, delta_D0, sm_eigvec_in, nij, fibersLength, NEI, NEO, f_mu; max_iter_nr=25, tolerance=1e-5, show_plot=false, meth=1)
            push!(path, zn)
            D0_init = D0_init + delta_D0  # Add delta_D0 at each iteration
            push!(path_D0, D0_init)
        end
        all_paths[l] = path
        all_D0s_paths[l] = path_D0
        all_D0s[l] = path_D0[end]

        # Update the progress bar at each iteration
        # next!(progress)
    end
    

    # Calculate the mean of the last two elements for each array in all_paths_th_0
    means_all_paths_th_0 = mean_last_two_elements(all_paths)
    # Calculate the mean of the last two elements for each array in all_D0s_paths
    means_all_D0s_paths = mean_last_two_elements(all_D0s_paths)

    # Calculate mesh_k as the minimum difference of the last two elements for each array in all_paths_th_0
    mesh_k = min_last_two_diff(all_paths)
    # Calculate mesh_D0 as the minimum difference of the last two elements for each array in all_D0s_paths
    mesh_D0 = min_last_two_diff(all_D0s_paths)


    success_array = [imag(all_paths[l][end]) >= 0 for l in eachindex(all_D0s_paths)]
    valid_indices = collect(filter(l -> success_array[l], eachindex(all_paths)))

    
    @floop for l in valid_indices
        all_paths[l][end], all_D0s_paths[l][end], _, _ = NR_last_point_gain(
            all_paths[l][end], 
            all_D0s_paths[l][end], 
            detLMatGain, 
            mesh_k, 
            mesh_D0; 
            tol=1e-9, max_iter=1000
        )
        all_D0s[l] = all_D0s_paths[l][end]
    end
    

    return all_paths, all_D0s, all_D0s_paths, success_array
end




# Plot the results
function plot_TLMs(all_paths, success_array, kroots)
    plt = plot(title="TLMs", xlabel="Real Part", ylabel="Imaginary Part", legend=false)
    
    for (path, success) in zip(all_paths, success_array)
        plot!(plt, real.(path), imag.(path), linestyle=:solid)
        if success
            scatter!(plt, [real(path[end])], [imag(path[end])], marker=:star, markersize=7, markercolor=:red)
        end
    end

    rr = real.(kroots)
    ir = imag.(kroots)
    scatter!(plt, rr, ir, legend=false, color=:blue)

    display(plt)
end

# candidate for the smallest D0
small_D0_kroot_index = argmin([abs(z - complex(k_a, 0.0)) for z in kroots]);
# evaluation of D0 for small_D0_kroot_index
begin  ##10 sec
    all_paths_th_0_small, all_D0s_th_0_small, all_D0s_paths_small, success_array_small = path_threshold(detLMatGain, [kroots[small_D0_kroot_index]], delta_D0, [sm_eigvec[small_D0_kroot_index]], nij, fibersLength, NEI, NEO, f_values; max_iter_pt=10^23, max_D0 = false, max_D0_value=100);  #(detLMatGain, [kroots[small_D0_kroot_index]], delta_D0; max_iter_pt=iter_D0,  max_D0 = false, max_D0_value=100)   
end

open(joinpath(dirname(@__FILE__), "out_data.txt"), "w") do file
    write(file, "length(all_D0s_paths_small[1]): $(length(all_D0s_paths_small[1])) must be large\n");
end

max_D0_est = 200*all_D0s_th_0_small[1];
plot_TLMs(all_paths_th_0_small, success_array_small, [kroots[small_D0_kroot_index]]);
delta_D0=all_D0s_th_0_small[1]/2;

# CHECK path_threshold
rr = real.(kroots);
ir = imag.(kroots);
colors = fill(:blue, length(kroots));
# colors[83] = :red
scatter(rr, ir, color=colors, legend=false);
# savefig(joinpath(dirname(@__FILE__), "scatter_plot_julia.pdf"))
# @timed begin # ~45 sec
# all_paths_th_0_try73, all_D0s_th_0_try73, all_D0s_paths_try73, success_array_small_try73 = path_threshold(detLMatGain, [kroots[73]], delta_D0, [sm_eigvec[73]], nij, fibersLength, NEI, NEO, f_values; max_iter_pt=10^23, max_D0 = false, max_D0_value=100)#(detLMatGain, [kroots[200]], delta_D0; max_iter_pt=iter_D0, max_D0 = true, max_D0_value=max_D0_est)   
# end
# all_D0s_th_0_try73[1]
# all_paths_th_0_try73[1][end]
# minimum(abs.(eigvals(LMatGain(all_paths_th_0_try73[1][end],all_D0s_th_0_try73[1]))))
# plot_TLMs(all_paths_th_0_try73, success_array_small_try73, [kroots[73]])
# all_paths_th_0_try73[1]
# all_D0s_paths_try73[1]
# val = [ minimum(abs.(eigvals(LMatGain(all_paths_th_0_try73[1][i],all_D0s_paths_try73[1][i])))) for i in eachindex(all_D0s_paths_try73[1]) ]
# maximum(val)
# minimum(abs.(eigvals(LMatGain(all_paths_th_0_try73[1][end],all_D0s_paths_try73[1][end]))))



##### BEGIN SAVED DATA #####
begin # ~1:40 h, 
    all_paths_th_0_ext, all_D0s_th_0_ext, all_D0s_paths_ext, success_array = path_threshold_bar(detLMatGain, kroots, delta_D0, sm_eigvec, nij, fibersLength, NEI, NEO, f_values;  max_D0 = true, max_D0_value=max_D0_est);  #(detLMatGain, kroots, delta_D0; max_iter_pt=iter_D0, max_D0 = true, max_D0_value=max_D0_est)
end
@save joinpath(dirname(@__FILE__), "out_all_paths_th_0_ext.jld2") all_paths_th_0_ext;
@save joinpath(dirname(@__FILE__), "out_all_D0s_th_0_ext.jld2") all_D0s_th_0_ext;
@save joinpath(dirname(@__FILE__), "out_all_D0s_paths_ext.jld2") all_D0s_paths_ext;
@save joinpath(dirname(@__FILE__), "out_success_array.jld2") success_array;
# @load joinpath(dirname(@__FILE__), "out_all_paths_th_0_ext.jld2") all_paths_th_0_ext
# @load joinpath(dirname(@__FILE__), "out_all_D0s_th_0_ext.jld2") all_D0s_th_0_ext
# @load joinpath(dirname(@__FILE__), "out_all_D0s_paths_ext.jld2") all_D0s_paths_ext
# @load joinpath(dirname(@__FILE__), "out_success_array.jld2") success_array
#### END SAVED DATA ###


restr_indices = collect(filter(l -> success_array[l], eachindex(kroots)));
kroots_ext = deepcopy(kroots);
kroots = [kroots[restr_indices[i]] for i in eachindex(restr_indices)];
all_paths_th_0 = [all_paths_th_0_ext[restr_indices[i]] for i in eachindex(restr_indices)];
all_D0s_th_0 = [all_D0s_th_0_ext[restr_indices[i]] for i in eachindex(restr_indices)];
all_D0s_paths = [all_D0s_paths_ext[restr_indices[i]] for i in eachindex(restr_indices)];


#maximum([abs(detLMatGain(all_paths_th_0[i][end],all_D0s_paths[i][end])) for i in eachindex(all_D0s_paths)])
mm_out_1 = maximum([minimum(abs.(eigvals(LMatGain(all_paths_th_0[i][end],all_D0s_paths[i][end])))) for i in eachindex(all_paths_th_0)]);
mm_out_2 = minimum([imag(all_paths_th_0[i][end]) for i in eachindex(all_paths_th_0)]);

open(joinpath(dirname(@__FILE__), "out_data.txt"), "w") do file
    write(file, "eignv and im_part: $mm_out_1, $mm_out_2 must be small\n");
end

maximum(all_D0s_th_0);
minimum(all_D0s_th_0);
minimum([length(all_D0s_paths[i]) for i in eachindex(all_D0s_paths)]);

plot_TLMs(all_paths_th_0_ext, success_array, kroots_ext);
savefig(joinpath(dirname(@__FILE__), "TLM_scatter_plot_julia.pdf"));

k_th_0 = [all_paths_th_0[i][end] for i in eachindex(all_paths_th_0)];
maximum(imag.(k_th_0));
if minimum(real.(k_th_0)) < 10^5*maximum(imag.(k_th_0))
    @warn "too large imag part of k_th_0"
end
k_th_0 = real.(k_th_0);
k_th_0_and_D0 = [[k_th_0[i], all_D0s_th_0[i]] for i in eachindex(k_th_0)];


######### EIGENVALUES WITH GAIN  ######

sort(eigvals(Array(LMatGain(all_paths_th_0[75][end],all_D0s_th_0[75]))), by = abs)[1:3];

function LMatGain_arr(v)
    LMatGain(v[1],v[2])
end

#= return the two smallest eigenvalues and the eigenvector corresponding to the smallest
 eigenvalue in presence of gain; meth=1 all eigenvalues calc; meth=2 Arnoldi =#
function smallest_eigen_gain(k_th_0_and_D0; meth=1)

    if meth == 1
        LMat_array = Array.(LMatGain_arr.(k_th_0_and_D0))
        E = Vector{LinearAlgebra.Eigen}(undef, length(k_th_0_and_D0))
        eigenvalues = Vector{Vector{ComplexF64}}(undef, length(k_th_0_and_D0))
        sorted_eigenvalues = Vector{Vector{ComplexF64}}(undef, length(k_th_0_and_D0))
        eigenvectors = Vector{Matrix{ComplexF64}}(undef, length(k_th_0_and_D0))
        pairs_eigenvalues = [Vector{ComplexF64}(undef, 2) for _ in eachindex(k_th_0_and_D0)]
        first_eigenvectors = Vector{Vector{ComplexF64}}(undef, length(k_th_0_and_D0))
        
        @threads for i in eachindex(k_th_0_and_D0)
            E[i] = eigen(LMat_array[i])
            eigenvalues[i] = E[i].values
            eigenvectors[i] = E[i].vectors
            sorted_eigenvalues[i] = sort(eigenvalues[i], by=abs)
            pairs_eigenvalues[i] = sorted_eigenvalues[i][1:2]
            first_eigenvectors[i] = eigenvectors[i][:, argmin(abs.(eigenvalues[i]))]
        end
    elseif meth == 2
        LMat_array = LMatGain_arr.(k_th_0_and_D0)
        results = Vector{Any}(undef, length(k_th_0_and_D0))
        
        @threads for i in eachindex(k_th_0_and_D0)
            results[i] = eigs(LMat_array[i], nev=2, which=:SM, ncv=40)
        end
        
        pairs_eigenvalues = [result[1] for result in results]
        first_eigenvectors = [result[2][:, 1] for result in results]
    end
    
    return pairs_eigenvalues, first_eigenvectors
end

#@timed begin 
sm_eigvls_gain, sm_eigvec_gain = smallest_eigen_gain(k_th_0_and_D0, meth=1);
#end

# It gives the norm of the vectors L(k)*eigevector(k)
values_check = [0.0 for i in eachindex(kroots)];
for i in eachindex(kroots)
    values_check[i] = norm(LMatGain_arr(k_th_0_and_D0[i])*sm_eigvec_gain[i]);
end
histogram(values_check, bins=100);   
# Compute the maximum and minimum norms
max_norm = maximum([norm(sm_eigvec_gain[i]) for i in eachindex(kroots)]);
min_norm = minimum([norm(sm_eigvec_gain[i]) for i in eachindex(kroots)]);
# Define the output file path
out_file_ = joinpath(@__DIR__, "out_data.txt");
# Append the results to the file
open(out_file_, "a") do io
    println(io, "Maximum norm: $max_norm");
    println(io, "Minimum norm: $min_norm");
end

###### Intensities calculation ##########

# Define the Psi_mu function
function Psi_mu(x, k_th_0_and_D0, sm_eigvec_gain, fibersLength, N_int, N_ext)
  
    k = k_th_0_and_D0[1]
    D0 = k_th_0_and_D0[2]

    # Calculate nij_g(k, D0)
    nij_values = nij_g(k, D0)
    
    # Initialize the result
    result = zeros(Complex{Float64}, N_int + N_ext)
    
    # Calculate the result for i = 1,...,N_int
    for i in 1:N_int
        lambda_plus = sm_eigvec_gain[i]
        lambda_minus = sm_eigvec_gain[N_int + N_ext + i]
        nij = nij_values[i]
        l_i = fibersLength[i]
        result[i] = lambda_plus * exp(im * k * nij * x) + lambda_minus * exp(im * k * nij * (l_i - x))
    end
    
    # Calculate the result for i = N_int+1,...,N_int+N_ext
    for i in N_int+1:N_int+N_ext
        lambda_plus = sm_eigvec_gain[i]
        nij = nij_values[i]
        result[i] = lambda_plus * exp(im * k * nij * x)
    end
    
    return result
end

### normalization

##### SAVED DATA ####
tol_int_lines = minimum([k_th_0_and_D0[l][1] for l in eachindex(k_th_0_and_D0)])*
    mean(fibersLength)*10^(-7);
n_sm_eigvec_gain = [[0.0+0.0im for i in 1:NE] for mu in eachindex(k_th_0_and_D0)];
begin # 1 min
for mu in eachindex(k_th_0_and_D0)
    norm_vec_ij = Complex{Float64}[Complex(i, 0) for i in 1:NE];
    int_lines_err_ij = Float64[i for i in 1.0:1.0:NE];
    @floop for i in 1:NE
        norm_vec_ij[i], int_lines_err_ij[i] = quadgk(x -> delta_pump[i] *
            Psi_mu(x, k_th_0_and_D0[mu], sm_eigvec_gain[mu], fibersLength, NEI, NEO)[i]^2, 
            0, fibersLength[i], atol=tol_int_lines);
    end
    if maximum(int_lines_err_ij) > tol_int_lines*10
        @warn "Error in the Psi_mu normalization!"
    end 
    norm_val = sum(norm_vec_ij);
    n_sm_eigvec_gain[mu] = sm_eigvec_gain[mu]/sqrt(norm_val);
end
end
# # # # # # # # #### Checking the normalization
check_norm_val = [0.0+0.0im for mu in eachindex(k_th_0_and_D0)];
for mu in eachindex(k_th_0_and_D0)
    norm_vec_ij = Complex{Float64}[Complex(i, 0) for i in 1:NE];
    int_lines_err_ij = Float64[i for i in 1.0:1.0:NE];
    @floop for i in 1:NE
        norm_vec_ij[i], int_lines_err_ij[i] = quadgk(x -> delta_pump[i] * 
            Psi_mu(x, k_th_0_and_D0[mu], n_sm_eigvec_gain[mu], fibersLength, NEI, NEO)[i]^2, 
            0, fibersLength[i], atol=tol_int_lines);
    end
    if maximum(int_lines_err_ij) > tol_int_lines * 10
        @warn "Error in the Psi_mu normalization!";
        
        # Define the output file path
        out_file = joinpath(@__DIR__, "out_data.txt");
        
        # Append the warning message to the file
        open(out_file, "a") do io
            println(io, "Warning: Error in the Psi_mu normalization!");
        end
    end
    check_norm_val[mu] = sum(norm_vec_ij);
end
@save joinpath(dirname(@__FILE__), "out_n_sm_eigvec_gain.jld2") n_sm_eigvec_gain;
@save joinpath(dirname(@__FILE__), "out_check_norm_val.jld2") check_norm_val;
# @load joinpath(dirname(@__FILE__), "out_n_sm_eigvec_gain.jld2") n_sm_eigvec_gain
# @load joinpath(dirname(@__FILE__), "out_check_norm_val.jld2") check_norm_val
#### END SAVED DATA #####


# both min and max ~ 1
# Compute the maximum and minimum of the absolute values
max_check_norm = maximum(abs.(check_norm_val));
min_check_norm = minimum(abs.(check_norm_val));
# Define the output file path
out_file_ = joinpath(@__DIR__, "out_data.txt");
# Append the results to the file
open(out_file_, "a") do io
    println(io, "Maximum absolute value of check_norm_val: $max_check_norm");
    println(io, "Minimum absolute value of check_norm_val: $min_check_norm");
end






######### T_matrix Analytical calculatin ########


function T_mat_calc_3(k_th_0_and_D0, n_sm_eigvec_gain,  k_a, gamma_perp, fibersLength, delta_pump, NEI, NEO, NE)

    N_roots = length(k_th_0_and_D0)

    vLeft = [[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im ] 
        for (nu, i) in product(1:N_roots, 1:NE)]

    @floop for (nu, i) in product(1:N_roots, 1:NEI) 
        vLeft[nu, i][1] = abs(n_sm_eigvec_gain[nu][i])^2
        vLeft[nu, i][2] = n_sm_eigvec_gain[nu][i] * conj(n_sm_eigvec_gain[nu][NEI + NEO + i])
        vLeft[nu, i][3] = conj(n_sm_eigvec_gain[nu][i]) * n_sm_eigvec_gain[nu][NEI + NEO + i]
        vLeft[nu, i][4] = abs(n_sm_eigvec_gain[nu][NEI + NEO + i])^2
    end
    @floop for (nu, i) in product(1:N_roots, (NEI+1):(NEI+NEO)) 
        vLeft[nu, i][1] = abs(n_sm_eigvec_gain[nu][i])^2
        vLeft[nu, i][2] = 0.0 + 0.0im
        vLeft[nu, i][3] = 0.0 + 0.0im
        vLeft[nu, i][4] = 0.0 + 0.0im
    end
    
    vRight = [[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im ] 
        for (nu, i) in product(1:N_roots, 1:NE)]

    @floop for (nu, i) in product(1:N_roots, 1:NEI) 
        vRight[nu, i][1] = n_sm_eigvec_gain[nu][i]^2
        vRight[nu, i][2] = n_sm_eigvec_gain[nu][i] * n_sm_eigvec_gain[nu][NEI + NEO + i]
        vRight[nu, i][3] = n_sm_eigvec_gain[nu][i] * n_sm_eigvec_gain[nu][NEI + NEO + i]
        vRight[nu, i][4] = n_sm_eigvec_gain[nu][NEI + NEO + i]^2
    end
    @floop for (nu, i) in product(1:N_roots, (NEI+1):(NEI+NEO)) 
        vRight[nu, i][1] = n_sm_eigvec_gain[nu][i]^2
        vRight[nu, i][2] = 0.0 + 0.0im
        vRight[nu, i][3] = 0.0 + 0.0im
        vRight[nu, i][4] = 0.0 + 0.0im
    end
    
    k_modif = [0.0+0.0im for (mu, i) in product(1:N_roots, 1:NE)]
    @floop for (mu, i) in product(1:N_roots, 1:NE)
        k_modif[mu, i] = k_th_0_and_D0[mu][1] * nij_g(k_th_0_and_D0[mu][1], k_th_0_and_D0[mu][2])[i]
    end
    
    AAux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    A1Aux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    BAux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    B1Aux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    CAux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    C1Aux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    DAux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    D1Aux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    EAux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    E1Aux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    FAux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    F1Aux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = 
            exp(im * (k_modif[nu, i] - conj(k_modif[nu, i]) + 2 * k_modif[mu, i]) * fibersLength[i]) - 1
        denominator = im * (k_modif[nu, i] - conj(k_modif[nu, i]) + 2 * k_modif[mu, i])
        if denominator == 0
            AAux[mu, nu, i] = fibersLength[i]
        else
            AAux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE) 
        numerator = exp(2im * k_modif[mu, i] * fibersLength[i]) * 
            (exp(im * (k_modif[nu, i]-conj(k_modif[nu, i])-2*k_modif[mu, i]) * fibersLength[i]) - 1)
        denominator = im * (k_modif[nu, i]-conj(k_modif[nu, i])-2*k_modif[mu, i])
        if denominator == 0
            BAux[mu, nu, i] = exp(2im * k_modif[mu, i] * fibersLength[i]) * fibersLength[i]
        else
            BAux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE) 
        numerator = exp(-im * conj(k_modif[nu, i]) * fibersLength[i]) *
            (exp(im*(conj(k_modif[nu,i]) + k_modif[nu,i] +2* k_modif[mu,i])*fibersLength[i]) - 1)
        denominator = im*(conj(k_modif[nu,i]) + k_modif[nu,i] +2* k_modif[mu,i])
        if denominator == 0
            CAux[mu, nu, i] = exp(-im * conj(k_modif[nu, i]) * fibersLength[i]) * fibersLength[i]
        else
            CAux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im * (2*k_modif[mu, i]-conj(k_modif[nu, i])) * fibersLength[i]) *
            (exp(im * (conj(k_modif[nu, i])+k_modif[nu, i] -2*k_modif[mu, i]) * fibersLength[i]) - 1)
        denominator = im * (conj(k_modif[nu, i])+k_modif[nu, i] -2*k_modif[mu, i])
        if denominator == 0
            DAux[mu, nu, i] =  exp(im * (2*k_modif[mu, i]-conj(k_modif[nu, i])) * fibersLength[i]) * fibersLength[i]
        else
            DAux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im  * k_modif[mu, i] * fibersLength[i]) * 
            (exp(im * (k_modif[nu, i]-conj(k_modif[nu, i])) * fibersLength[i]) - 1)
        denominator = im * (k_modif[nu, i]-conj(k_modif[nu, i])) 
        if denominator == 0
            EAux[mu, nu, i] = exp(im  * k_modif[mu, i] * fibersLength[i]) * fibersLength[i]
        else
            EAux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im * (k_modif[mu, i] - conj(k_modif[nu, i])) * fibersLength[i]) *
            (exp(im * (k_modif[nu, i]+conj(k_modif[nu, i])) * fibersLength[i]) - 1)
        denominator = im * (k_modif[nu, i]+conj(k_modif[nu, i]))
        if denominator == 0
            FAux[mu, nu, i] = exp(im * (k_modif[mu, i] - conj(k_modif[nu, i])) * fibersLength[i]) * fibersLength[i]
        else 
            FAux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im * k_modif[nu, i] * fibersLength[i]) *
            (exp(im * (2*k_modif[mu, i]-k_modif[nu, i]-conj(k_modif[nu, i])) * fibersLength[i]) - 1)
        denominator = im * (2*k_modif[mu, i]-k_modif[nu, i]-conj(k_modif[nu, i]))
        if denominator == 0
            D1Aux[mu, nu, i] = exp(im * k_modif[nu, i] * fibersLength[i]) * fibersLength[i]
        else 
            D1Aux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im * (k_modif[mu, i] + k_modif[nu, i]) * fibersLength[i]) *
            (exp(-im * (k_modif[nu, i]+conj(k_modif[nu, i])) * fibersLength[i]) - 1)
        denominator = -im * (k_modif[nu, i]+conj(k_modif[nu, i]))
        if denominator == 0
            F1Aux[mu, nu, i] = exp(im * (k_modif[mu, i] + k_modif[nu, i]) * fibersLength[i])  * fibersLength[i]
        else 
            F1Aux[mu, nu, i] = numerator/denominator
        end
    end


    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im * (2*k_modif[mu, i]+k_modif[nu, i]) * fibersLength[i]) *
            (exp(-im * (2*k_modif[mu, i]+k_modif[nu, i]+conj(k_modif[nu, i])) * fibersLength[i]) - 1)
        denominator = -im * (2*k_modif[mu, i]+k_modif[nu, i]+conj(k_modif[nu, i]))
        if denominator == 0
            C1Aux[mu, nu, i] = exp(im * (2*k_modif[mu, i]+k_modif[nu, i]) * fibersLength[i])  * fibersLength[i]
        else 
            C1Aux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im * (k_modif[nu, i]-conj(k_modif[nu, i])) * fibersLength[i]) *
            (exp(im * (2*k_modif[mu, i]+conj(k_modif[nu, i])-k_modif[nu, i]) * fibersLength[i]) - 1)
        denominator = im * (2*k_modif[mu, i]+conj(k_modif[nu, i])-k_modif[nu, i])
        if denominator == 0
            B1Aux[mu, nu, i] = exp(im * (k_modif[nu, i]-conj(k_modif[nu, i])) * fibersLength[i])  * fibersLength[i]
        else 
            B1Aux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im * (k_modif[mu, i]+k_modif[nu, i]-conj(k_modif[nu, i])) * fibersLength[i]) *
            (exp(im * (conj(k_modif[nu, i])-k_modif[nu, i]) * fibersLength[i]) - 1)
        denominator = im * (conj(k_modif[nu, i])-k_modif[nu, i])
        if denominator == 0
            E1Aux[mu, nu, i] = exp(im * (k_modif[mu, i]+k_modif[nu, i]-conj(k_modif[nu, i])) * fibersLength[i]) * fibersLength[i]
        else 
            E1Aux[mu, nu, i] = numerator/denominator
        end
    end

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        numerator = exp(im * (2*k_modif[mu, i]+k_modif[nu, i]-conj(k_modif[nu, i])) * fibersLength[i]) *
            (exp(im * (conj(k_modif[nu, i])-k_modif[nu, i]-2*k_modif[mu, i]) * fibersLength[i]) - 1)
        denominator = im * (conj(k_modif[nu, i])-k_modif[nu, i]-2*k_modif[mu, i])
        if denominator == 0
            A1Aux[mu, nu, i] = exp(im * (2*k_modif[mu, i]+k_modif[nu, i]-conj(k_modif[nu, i])) * fibersLength[i]) * fibersLength[i]
        else 
            A1Aux[mu, nu, i] = numerator/denominator
        end
    end
    
    IntMatr = [[
        0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;
        0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;
        0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;
        0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im] 
        for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]

    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE) 
        IntMatr[mu, nu, i][1, 1] = AAux[mu, nu, i]
        IntMatr[mu, nu, i][1, 2] = EAux[mu, nu, i]
        IntMatr[mu, nu, i][1, 3] = EAux[mu, nu, i]
        IntMatr[mu, nu, i][1, 4] = BAux[mu, nu, i]
        IntMatr[mu, nu, i][2, 1] = CAux[mu, nu, i]
        IntMatr[mu, nu, i][2, 2] = FAux[mu, nu, i]
        IntMatr[mu, nu, i][2, 3] = FAux[mu, nu, i]
        IntMatr[mu, nu, i][2, 4] = DAux[mu, nu, i]
        IntMatr[mu, nu, i][3, 1] = D1Aux[mu, nu, i]
        IntMatr[mu, nu, i][3, 2] = F1Aux[mu, nu, i]
        IntMatr[mu, nu, i][3, 3] = F1Aux[mu, nu, i]
        IntMatr[mu, nu, i][3, 4] = C1Aux[mu, nu, i]
        IntMatr[mu, nu, i][4, 1] = B1Aux[mu, nu, i]
        IntMatr[mu, nu, i][4, 2] = E1Aux[mu, nu, i]
        IntMatr[mu, nu, i][4, 3] = E1Aux[mu, nu, i]
        IntMatr[mu, nu, i][4, 4] = A1Aux[mu, nu, i]
    end

    TAux = [0.0 + 0.0im for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)]
    @floop for (mu, nu, i) in product(1:N_roots, 1:N_roots, 1:NE)
        TAux[mu, nu, i] = sum(vLeft[nu,i][j]*(IntMatr[mu, nu, i]*vRight[mu, i])[j] for j in 1:4) #dot(vLeft[nu, i], IntMatr[mu, nu, i], vRight[mu, i])
    end

    T_complex = Matrix{Complex{Float64}}(undef, N_roots, N_roots)
    @floop for (mu, nu) in product(1:N_roots, 1:N_roots)
        T_complex[mu,nu] = gamma_perp^2/((k_th_0_and_D0[nu][1] - k_a)^2+gamma_perp^2) * 
            sum(delta_pump[i] * TAux[mu, nu, i] for i in 1:NE)
    end


    return T_complex
end




T_complex = T_mat_calc_3(k_th_0_and_D0, n_sm_eigvec_gain, k_a, gamma_perp, fibersLength, delta_pump, NEI, NEO, NE);




############################



##### real matrx approximation
T_real_complete = real.(T_complex);


T_real = T_real_complete;



minimum(all_D0s_th_0);
# we have to consider all the modes with non-interact D_0 much higher than the first
if maximum(all_D0s_th_0) < 80 * minimum(all_D0s_th_0)
    ratio = maximum(all_D0s_th_0) / minimum(all_D0s_th_0);
    @warn "Low maximum(all_D0s_th_0): max/min = $ratio. More TLMs should be included!";
    
    # Define the output file path
    out_file_ = joinpath(@__DIR__, "out_data.txt");
    
    # Append the warning message to the file
    open(out_file_, "a") do io
        println(io, "Warning: Low maximum(all_D0s_th_0). max/min = $ratio. More TLMs should be included!");
    end
end




function modes_sequence(T_real, T_complex, all_D0s_th_0)

    # Initializations
    # number of modes considered
    M = length(all_D0s_th_0)
    # initialization of partially ordered T-matrix
    T_partially_ordered = deepcopy(T_real)
    # initialization partially ordered non-int threshold 
    all_D0s_th_0_partially_ord = deepcopy(all_D0s_th_0)
    # for the first mode, modes interaction vanishes
    ordered_modes = [argmin(all_D0s_th_0)]
    all_D0s_int = [ all_D0s_th_0[ordered_modes[1]] ]
    N = length(ordered_modes)
    # initialization partially ordered modes: [1,2,...,M]
    partially_ordered_modes = collect(1:M)
    # restricted ordered T-matrix initialization
    T_ordered_restricted = Matrix{Float64}(undef, N, N)
    # ordered complex T-matrix initialization
    T_ordered_restricted_complex = Matrix{Complex{Float64}}(undef, N, N)


    while N < M
        
        N = length(ordered_modes)

        # permutation of N-th mode
        partially_ordered_modes[N] = ordered_modes[N]
        partially_ordered_modes[ordered_modes[N]] = N

        # restricted ordered T-matrix initialization
        T_ordered_restricted = Matrix{Float64}(undef, N, N)
        # ordered T-matrix definition
        for (i,j) in product(eachindex(ordered_modes), eachindex(ordered_modes))
            T_ordered_restricted[i,j] = T_real[ordered_modes[i],ordered_modes[j]]
        end
        # ordered complex T-matrix initialization
        T_ordered_restricted_complex = Matrix{Complex{Float64}}(undef, N, N)

        # checking if T-matrix has large imaginary part or negative real part
        # ordered complex T-matrix
        for (i,j) in product(eachindex(ordered_modes), eachindex(ordered_modes))
            T_ordered_restricted_complex[i,j] = T_complex[ordered_modes[i],ordered_modes[j]]
        end
        check_T_ratio_restr = [0.0 for (_,_) in product(eachindex(ordered_modes), eachindex(ordered_modes))]
        for (mu, nu) in product(eachindex(ordered_modes), eachindex(ordered_modes))
            check_T_ratio_restr[mu, nu] = abs.(imag.(T_ordered_restricted_complex))[mu, nu]/abs.(real.(T_ordered_restricted_complex))[mu, nu]
        end
        if maximum([check_T_ratio_restr[mu, N] for mu in 1:N]) > 0.5 || maximum([check_T_ratio_restr[N, nu] for nu in 1:N]) > 0.5
            warning_message = "T has large imaginary part for the $(N)-th mode: $([T_ordered_restricted_complex[mu, N] for mu in 1:N])"
            @warn warning_message
        
            # Define the output file path
            out_file = joinpath(@__DIR__, "out_data.txt")
        
            # Append the warning message to the file
            open(out_file, "a") do io
                println(io, "Warning: $warning_message")
            end
        end
        
        if minimum(real.([T_ordered_restricted_complex[mu, N] for mu in 1:N])) < 0 || 
            minimum(real.([T_ordered_restricted_complex[N, nu] for nu in 1:N])) < 0
         
             warning_message = "T has negative elements for the $(N)-th mode: $(real.([T_ordered_restricted_complex[mu, N] for mu in 1:N]))"
             @warn warning_message
         
             # Define the output file path
             out_file = joinpath(@__DIR__, "out_data.txt")
         
             # Append the warning message to the file
             open(out_file, "a") do io
                 println(io, "Warning: $warning_message")
             end
         end
         

        # inverse of the ordered T-matrix
        # checking clamping-gain
        if det(T_ordered_restricted) == 0
            pop!(ordered_modes)
            pop!(all_D0s_int)
            # Define the message
            message = "Clamping-gain: N_max = $(N)."
            # Print the message to the console
            println(message)
            # Define the output file path
            out_file = joinpath(@__DIR__, "out_data.txt")
            # Append the message to the file
            open(out_file, "a") do io
                println(io, message)
            end

            break
        end
        T_ordered_restricted_inv = inv(T_ordered_restricted)

        # permutated non-int threshold
        for i in eachindex(partially_ordered_modes)
            all_D0s_th_0_partially_ord[i] = all_D0s_th_0[partially_ordered_modes[i]]
        end
        # permutated T-matrix
        for (i,j) in product(eachindex(partially_ordered_modes), eachindex(partially_ordered_modes))
            T_partially_ordered[i,j] = T_real[partially_ordered_modes[i],partially_ordered_modes[j]]
        end

        # initialization threshold attempt
        D_int_attempt = Vector{Float64}(undef, M-N)
        D_int_attempt_num = Vector{Float64}(undef, M-N)
        D_int_attempt_den = Vector{Float64}(undef, M-N)
        # threshold attempt equation
        for mu in 1:M-N
            D_int_attempt_num[mu] = (1-sum(T_partially_ordered[N+mu,i]*T_ordered_restricted_inv[i,j] for i in 1:N, j in 1:N))
            D_int_attempt_den[mu] = (1-all_D0s_th_0_partially_ord[N+mu]*sum(T_partially_ordered[N+mu,i]*T_ordered_restricted_inv[i,j]/
                all_D0s_th_0_partially_ord[j] for i in 1:N, j in 1:N))
            D_int_attempt[mu] = all_D0s_th_0_partially_ord[N+mu] * D_int_attempt_num[mu] / D_int_attempt_den[mu]     
        end

        # selecting next all_D0s_int (interacting D0s)
        # Get the last value of all_D0s_int
        threshold_value = all_D0s_int[end]
        # Filter the elements of D_int_attempt that are greater than threshold_value
        filtered_elements = filter(x -> x > threshold_value, D_int_attempt)
        # Find the index of the smallest element among the filtered ones
        if !isempty(filtered_elements)
            min_value = minimum(filtered_elements)
            min_index = findfirst(x -> x == min_value, D_int_attempt)
            if D_int_attempt_den[min_index] < 0
                # Define the message
                message = "Clamping-gain: N_max = $(N)."
                # Print the message to the console
                println(message)
                # Define the output file path
                out_file = joinpath(@__DIR__, "out_data.txt")
                # Append the message to the file
                open(out_file, "a") do io
                    println(io, message)
                end
                break
            end
        else
            # Define the message
            message = "Clamping-gain: N_max = $(N)."
            # Print the message to the console
            println(message)
            # Define the output file path
            out_file = joinpath(@__DIR__, "out_data.txt")
            # Append the message to the file
            open(out_file, "a") do io
                println(io, message)
            end

            break
        end
        push!(ordered_modes, N + min_index)
        push!(all_D0s_int, min_value)

    end


    return ordered_modes, all_D0s_int, T_ordered_restricted, T_ordered_restricted_complex

end
 



ordered_modes, all_D0s_int, T_ordered_restricted, T_ordered_restricted_complex =  modes_sequence(T_real, T_complex, all_D0s_th_0);
T_ordered_restricted_complex;
ordered_modes;
all_D0s_int;
T_ordered_restricted;

@save joinpath(@__DIR__, "out_all_D0s_th_0.jld2") all_D0s_th_0;
@save joinpath(@__DIR__, "out_kroots.jld2") kroots;
@save joinpath(@__DIR__, "out_ordered_modes.jld2") ordered_modes;

# Normalize the color values based on min and max of all_D0s_th_0
color_min, color_max = minimum(all_D0s_th_0), maximum(all_D0s_th_0);

# Separate indices for dots and stars
dots_indices = setdiff(1:length(kroots), ordered_modes);  # Indices for dots
stars_indices = ordered_modes;  # Indices for stars

# Create scatter plot for dots
p = scatter(
    real(kroots[dots_indices]), imag(kroots[dots_indices]),  # Real and imaginary parts of dots
    zcolor=all_D0s_th_0[dots_indices],  # Color values based on normalized zcolor
    label="Dots",  # Legend label for dots
    ms=8,  # Marker size for dots
    c=:viridis,  # Color map
    colorbar=true  # Add color bar
);

# Overlay scatter plot for stars
scatter!(
    p,
    real(kroots[stars_indices]), imag(kroots[stars_indices]),  # Real and imaginary parts of stars
    zcolor=all_D0s_th_0[stars_indices],  # Color values based on normalized zcolor
    marker=:star5,  # Marker shape for stars
    label="Stars",  # Legend label for stars
    ms=12,  # Marker size for stars
    c=:viridis  # Color map
);

# Add labels and title
xlabel!("Real Part");  # Label for the x-axis
ylabel!("Imaginary Part");  # Label for the y-axis
title!("Complex Roots with Color-coded all_D0s_th_0");  # Plot title
# Save the plot in the same directory as the script
savefig(joinpath(@__DIR__, "complex_roots_plot.pdf"));  # Save the plot as a PNG file





check_T_ratio_restr = [0.0 for (_,_) in product(eachindex(ordered_modes), eachindex(ordered_modes))];
for (mu, nu) in product(eachindex(ordered_modes), eachindex(ordered_modes))
    check_T_ratio_restr[mu, nu] = abs.(imag.(T_ordered_restricted_complex))[mu, nu]/abs.(real.(T_ordered_restricted_complex))[mu, nu];
end

maximum(check_T_ratio_restr);
heatmap(check_T_ratio_restr, title="im/real Heatmap", xlabel="Column", ylabel="Row", color=:viridis);
heatmap(T_ordered_restricted, title="T real matrx", xlabel="Column", ylabel="Row", color=:viridis);



# example of calculation

all_D0s_int;
# Calculate D0
D0 = 1.05 * maximum(all_D0s_int);
@save joinpath(@__DIR__, "out_D0.jld2") D0;
# Define the message
message = "D0 = 1.05 * max = $D0";
# Print the message to the console
println(message);
# Define the output file path
out_file_ = joinpath(@__DIR__, "out_data.txt");
# Append the message to the file
open(out_file_, "a") do io
    println(io, message);
end




function spectrum_int(D0,all_D0s_int,k_th_0,T_real,ordered_modes)
    # number of modes involved
    N = findlast(x -> x < D0, all_D0s_int)

    # restricted ordered T-matrix initialization
    T_ordered_restricted = Matrix{Float64}(undef, N, N)
    # ordered T-matrix definition
    for (i, j) in product(1:N, 1:N)
        T_ordered_restricted[i,j] = T_real[ordered_modes[i], ordered_modes[j]]
    end
    T_inverse = inv(T_ordered_restricted)
    ordered_D0s_th_0 = [all_D0s_th_0[ordered_modes[i]] for i in 1:N]
    ordered_k_th_0 = [k_th_0[ordered_modes[i]] for i in 1:N]

    # vector of D0 threshold
    D0_vector = Vector{Float64}(undef, N)
    for mu in 1:N
        D0_vector[mu] = D0 / ordered_D0s_th_0[mu] - 1
    end
    I_try0 = T_inverse * D0_vector


    ################# Negative remov
    # Find indices of negative elements in I_try0
    negative_indices = findall(<(0), I_try0)

    # If there are negative indices, process the matrices and I_try0
    if !isempty(negative_indices)
        # Remove rows and columns corresponding to negative indices
        T_inverse = T_inverse[Not(negative_indices), Not(negative_indices)]
        D0_vector = D0_vector[Not(negative_indices)]

        # Recalculate I_try0
        I_try0 = T_inverse * D0_vector

        # Append zeros for previously removed indices
        full_I_try0 = zeros(length(I_try0) + length(negative_indices))
        valid_indices = setdiff(1:length(full_I_try0), negative_indices)
        full_I_try0[valid_indices] = I_try0
        I_try0 = full_I_try0
    end
    ##################################################

    ######### saving data
    @save joinpath(@__DIR__, "out_ordered_k_th_0.jld2") ordered_k_th_0
    @save joinpath(@__DIR__, "out_I_try0.jld2") I_try0

    # Initialize the plot
    plt = plot()
    # Add vertical lines
    for i in 1:length(ordered_k_th_0)
        plot!([ordered_k_th_0[i], ordered_k_th_0[i]], [0, I_try0[i]], label="", color=:blue)
    end
    # Add a vertical dashed red line at x = k_a
    plot!([k_a, k_a], [0, maximum(I_try0)], label="", color=:red, linestyle=:dash)
    # Add label for the series in the legend
    plot!(ordered_k_th_0, I_try0, seriestype=:scatter, label="Int., D0=$D0", color=:blue)
    # Set the x-axis limits to ensure Re_k_min and Re_k_max are displayed
    plot!(xlims=(Re_k_min, Re_k_max))

    # Save the plot in the script directory
    output_path = joinpath(@__DIR__, "spectrum_plot.pdf")
    savefig(plt, output_path)

    # Display the plot explicitly
    display(plt)

    return I_try0
end



I_0 = spectrum_int(D0,all_D0s_int,k_th_0,T_real,ordered_modes);

function plot_fibers_intens(D0,all_D0s_int,k_th_0,T_real,ordered_modes)
    # number of modes involved
    N = findlast(x -> x < D0, all_D0s_int)

    # restricted ordered T-matrix initialization
    T_ordered_restricted = Matrix{Float64}(undef, N, N)
    # ordered T-matrix definition
    for (i, j) in product(1:N, 1:N)
        T_ordered_restricted[i,j] = T_real[ordered_modes[i], ordered_modes[j]]
    end
    T_inverse = inv(T_ordered_restricted)
    ordered_D0s_th_0 = [all_D0s_th_0[ordered_modes[i]] for i in 1:N]
    ordered_k_th_0 = [k_th_0[ordered_modes[i]] for i in 1:N]

    # vector of D0 threshold
    D0_vector = Vector{Float64}(undef, N)
    for mu in 1:N
        D0_ordered_k_th_0vector[mu] = D0 / ordered_D0s_th_0[mu] - 1
    end
    a_mu = ordered_k_th_0sqrt.(T_inverse * D0_vector)

end



ordered_modes;
# number of modes involved
N = findlast(x -> x < D0, all_D0s_int);

# restricted ordered T-matrix initialization
T_ordered_restricted = Matrix{Float64}(undef, N, N);
# ordered T-matrix definition
for (i, j) in product(1:N, 1:N)
    T_ordered_restricted[i,j] = T_real[ordered_modes[i], ordered_modes[j]];
end
T_inverse = inv(T_ordered_restricted);
ordered_D0s_th_0 = [all_D0s_th_0[ordered_modes[i]] for i in 1:N];
ordered_k_th_0 = [k_th_0[ordered_modes[i]] for i in 1:N];

# vector of D0 threshold
D0_vector = Vector{Float64}(undef, N);
for mu in 1:N
    D0_vector[mu] = D0 / ordered_D0s_th_0[mu] - 1;
end
I_0 = T_inverse * D0_vector;
# Combine the data into a matrix
data = hcat(ordered_k_th_0, I_0);

# Save the data to a CSV file in the same directory as the script
csv_path = joinpath(@__DIR__, "out_spectrum.csv");
open(csv_path, "w") do io
    for row in eachrow(data)
        println(io, join(row, ","));
    end
end

# Path to save the D0 data
d0_csv_path = joinpath(@__DIR__, "out_D0.csv");

# Save D0 to a CSV file
open(d0_csv_path, "w") do io
    for value in D0
        println(io, value);
    end
end




function fibers_intens(x,I_0,k_th_0_and_D0,n_sm_eigvec_gain,ordered_modes,fibersLength, NEI, NEO; mod=1)
    i = mod
    
    res = [abs(sqrt(I_0[i])*Psi_mu(x, k_th_0_and_D0[ordered_modes[i]], n_sm_eigvec_gain[ordered_modes[i]], 
    fibersLength, NEI, NEO)[j])^2 for j in eachindex(fibers)]

    return res
end







########################

function plot_fibers_intens(D0, all_D0s_int, k_th_0, T_real, ordered_modes, I_0, k_th_0_and_D0, n_sm_eigvec_gain, fibersLength, NEI, NEO)
    # Create the plot
    p = plot()

    # Find the maximum length among the fibers
    max_fiber_length = maximum(fibersLength)

    # Set the base number of points for the longest fiber
    base_num_points = 150

    # Calculate the maximum and minimum intensity to normalize the color values
    max_intensity = -Inf
    min_intensity = Inf

    for (i, _) in enumerate(fibersPoints)
        xs = range(0, stop=fibersLength[i], length=100)
        for x in xs
            fiber_intensity = fibers_intens(x, I_0, k_th_0_and_D0, n_sm_eigvec_gain, ordered_modes, fibersLength, NEI, NEO)[i]
            max_intensity = max(max_intensity, maximum(fiber_intensity))
            min_intensity = min(min_intensity, minimum(fiber_intensity))
        end
    end

    # Now plot with normalized intensity values using scatter
    for (i, points) in enumerate(fibersPoints)
        # Calculate the number of points proportionally to the fiber length
        num_samples = maximum([2, Int(round(base_num_points * (fibersLength[i] / max_fiber_length)))])

        # Define a set of sample points along the fiber
        xs = range(0, stop=fibersLength[i], length=num_samples)

        # Get the coordinates of the two points defining the fiber
        x_start, y_start = points[1]
        x_end, y_end = points[2]

        # Generate the scatter points along the fiber
        x_vals = [(1 - t) * x_start + t * x_end for t in xs / fibersLength[i]]
        y_vals = [(1 - t) * y_start + t * y_end for t in xs / fibersLength[i]]

        # Calculate intensity for each point along the fiber
        intensities = [fibers_intens(x, I_0, k_th_0_and_D0, n_sm_eigvec_gain, ordered_modes, fibersLength, NEI, NEO)[i][1] for x in xs]

        # Normalize the intensities (for color mapping) but keep real values for the colorbar
        norm_intensities = [(intensity - min_intensity) / (max_intensity - min_intensity) for intensity in intensities]

        # Plot the scatter points with color based on intensity
        scatter!(x_vals, y_vals, zcolor=intensities, c=:inferno, clim=(min_intensity, max_intensity), label="", markerstrokewidth=0, markersize=1)
    end

    # Finalize the plot with equal aspect ratio and automatic colorbar
    plot!(legend=false, aspect_ratio=:equal, colorbar_title="Intensity")

    # Save the plot to the same directory as the script
    output_path = joinpath(@__DIR__, "fibers_intensity_plot.pdf")
    savefig(p, output_path)

    # Optionally display the plot (for interactive use)
    display(p)
end


# Call the function
plot_fibers_intens(D0, all_D0s_int, k_th_0, T_real, ordered_modes, I_0, k_th_0_and_D0, n_sm_eigvec_gain, fibersLength, NEI, NEO);


