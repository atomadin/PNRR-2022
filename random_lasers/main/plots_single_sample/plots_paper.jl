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


using CSV, DataFrames, JLD2
using LinearAlgebra
using Plots
using FLoops, Base.Iterators
using LaTeXStrings
pyplot()


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

# Path to the SAMPLE_NAME.txt file
sample_file = joinpath(@__DIR__, "SAMPLE_NAME.txt")

# Read the file and ignore comment lines starting with "#"
sample_name = filter(!startswith("#"), split(read(sample_file, String), '\n'))
sample_name = strip(sample_name[1])  # first valid line without spaces

# Build the directory path using the selected sample name
directory_name = joinpath(dirname(dirname(@__FILE__)), sample_name)

include(joinpath(directory_name,"InputData/input_data.txt"));



#########    Graph Definition   #############
# Vertices in the spot
vertIn = [row |> Vector for row in 
    eachrow(CSV.read(joinpath(directory_name, "InputData/vertIn.csv"), DataFrame, header=false))];
# Vertices out of the spot
vertOut = [row |> Vector for row in 
    eachrow(CSV.read(joinpath(directory_name, "InputData/vertOut.csv"), DataFrame, header=false))];
# Internal Edges
fibersIn = [row |> Vector for row in 
    eachrow(CSV.read(joinpath(directory_name, "InputData/fibersIn.csv"), DataFrame, header=false))];
# External Edges
fibersOut = [row |> Vector for row in 
    eachrow(CSV.read(joinpath(directory_name, "InputData/fibersOut.csv"), DataFrame, header=false))];
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
NVIn,  NVOut, vert, NV, NEI, NEO, fibers, NE, fibersLength, fibersPoints = 
    derived_params(vertIn,vertOut,fibersIn,fibersOut);

# Write to the file
file_path = joinpath(@__DIR__, "out_numb.txt")
open(file_path, "w") do file
    write(file, "$NVIn\n")
    write(file, "$NE\n")
end


# Number of bins using Sturges' rule
bins_sturges = ceil(Int, log2(length(fibersLength)) + 1)
# Plot histogram
histogram(
    fibersLength,
    bins=bins_sturges,
    title="Fiber Lengths, #bins = $bins_sturges, NE = $NE, NV = $NV",
    xlabel="Length",
    ylabel="Frequency",
    legend=false
)
# Save the plot
output_path = joinpath(dirname(@__FILE__),"fiber_lengths_histogram.pdf")  # Change this to your preferred file name and path
savefig(output_path)

###### Refractive Indices ##############
nij = fill(1.5, NE) ; # refractive indices
delta_pump = fill(1.0, NE)




function plotFiberGraph(fibersPoints, vertIn, vertOut)
    plt = plot()

    # Add lines representing the fiber segments
    for points in fibersPoints
        plot!([points[1][1], points[2][1]], [points[1][2], points[2][2]], 
              color=:black, linewidth=1, label="")
    end

    # Add a red semi-transparent disk of radius 40 centered at (0,0)
    θ = range(0, 2π; length=200)
    r = 40.0
    x_disk = r * cos.(θ)
    y_disk = r * sin.(θ)
    plot!(x_disk, y_disk, fill=(1, :red, 0.3), linecolor=:transparent, label="")

    # Axis labels and limits
    xlabel = L"x\ \mathrm{(\mu m)}"
    ylabel = L"y\ \mathrm{(\mu m)}"
    plot!(xlabel=xlabel, ylabel=ylabel)

    # Axis limits (adjust as needed for your case)
    plot!(xlims=(-43, 43), ylims=(-43, 43))

    # Final plot adjustments
    plot!(legend=false, aspect_ratio=:equal, framestyle=:box)

    # Save and display
    savefig(plt, joinpath(@__DIR__, "fibers.pdf"))
    display(plt)
end


plotFiberGraph(fibersPoints, vertIn, vertOut)
# Path to the PDF file created by the function
pdf_path = joinpath(@__DIR__, "fibers.pdf")
# Apply cropping using pdfcrop (overwrites the original file)
run(`pdfcrop $pdf_path $pdf_path`)


###################################################
## Quasi-Bound
##################################################

@load joinpath(directory_name, "out_kroots.jld2") kroots
@load joinpath(directory_name, "out_ordered_modes.jld2") ordered_modes
@load joinpath(directory_name, "out_all_D0s_th_0.jld2") all_D0s_th_0



# # Separate indices for dots and stars
dots_indices = setdiff(1:length(kroots), ordered_modes)  # Indices for dots
stars_indices = ordered_modes  # Indices for stars




# Scatter plot of all kroots as dark blue dots
p = scatter(
    real(kroots), imag(kroots),
    label = "",
    ms = 8,
    color = :darkblue
)

# Add vertical dashed line at k_a = 10.68
vline!(p, [10.68], linestyle = :dash, color = :darkred, label="")

# Get current xticks (positions only)
xtick_vals = Plots.xticks(p)[1][1]

# Add 10.68 if not already present
xtick_positions = sort(union(xtick_vals, [10.68]))

# Create corresponding labels, replacing 10.68 with "kₐ"
xtick_labels = [abs(x - 10.68) < 1e-3 ? L"k_a" : string(x) for x in xtick_positions]

# Apply updated ticks (now: Vector{Float64}, Vector{String})
xticks!(p, xtick_positions, xtick_labels)

# Axis labels and title in LaTeX
xlabel!(L"\Re\,k^\mu\, (\mu\text{m}^{-1})")
ylabel!(L"\Im\,k^\mu\, (\mu\text{m}^{-1})")
title!(L"\text{(a) Quasi-bound } k^\mu \text{ modes}")

# Save plot as PDF
savefig(p, joinpath(@__DIR__, "quasibound_modes.pdf"))


#######################################
##### TLM
######################################

@load joinpath(directory_name, "out_all_paths_th_0_ext.jld2") all_paths_th_0_ext
@load joinpath(directory_name, "out_success_array.jld2") success_array
restr_indices = collect(filter(l -> success_array[l], eachindex(kroots)))
all_paths_th_0 = [all_paths_th_0_ext[restr_indices[i]] for i in eachindex(restr_indices)]
kroots_ext = deepcopy(kroots)


function plot_TLMs(all_paths, success_array, kroots, all_D0s_th_0, dots_indices, stars_indices)

    # draw TLM paths (thin dark blue lines)
    p = plot()
    for (path, _) in zip(all_paths, success_array)
        plot!(
            p,
            real.(path), imag.(path),
            linestyle = :dash,
            color = :darkblue,
            linewidth = 0.01,  # Thin lines
            label = nothing,   # No legend entry
            legend = :bottomleft
        )
    end

    # Overlay scatter plot for non-activated modes (dots)
    scatter!(
        p,
        real(kroots[dots_indices]), imag(kroots[dots_indices]),
        zcolor = all_D0s_th_0[dots_indices],
        ms = 8,
        marker = :circle,
        label = "Quasi-bound modes",
        c = :viridis,
        colorbar = true,
        colorbar_title = L"D_0^{\text{th}}"
    )

    # Overlay scatter plot for activated modes (stars)
    scatter!(
        p,
        real(kroots[stars_indices]), imag(kroots[stars_indices]),
        zcolor = all_D0s_th_0[stars_indices],
        marker = :star5,
        label = "Activated modes",
        ms = 12,
        c = :viridis
    )

    # Axis labels and plot title (in LaTeX)
    xlabel!(L"\Re\,k^\mu\, (\mu\text{m}^{-1})")
    ylabel!(L"\Im\,k^\mu\, (\mu\text{m}^{-1})")
    title!(L"k^\mu \text{ Modes}")

    first = true
    for (path, success) in zip(all_paths, success_array)
        if success
            scatter!(p,
                     [real(path[end])], [imag(path[end])],
                     marker = :plus,
                     markersize = 5,
                     markercolor = :darkred,
                     label = first ? "TLMs" : "")
            first = false
        end
    end

    # Add vertical dashed line at k_a = 10.68
    vline!(p, [10.68], linestyle = :dash, color = :darkred, label="")


    # Get current xticks (positions only)
    xtick_vals = Plots.xticks(p)[1][1]

    # Add 10.68 if not already present
    xtick_positions = sort(union(xtick_vals, [10.68]))

    # Create corresponding labels, replacing 10.68 with "kₐ"
    xtick_labels = [abs(x - 10.68) < 1e-3 ? L"k_a" : string(x) for x in xtick_positions]

    # Apply updated ticks (now: Vector{Float64}, Vector{String})
    xticks!(p, xtick_positions, xtick_labels)


    display(p)
end




plot_TLMs(all_paths_th_0, success_array, kroots, all_D0s_th_0, dots_indices, stars_indices)


savefig(joinpath(dirname(@__FILE__), "Modes.pdf"))

input_pdf = joinpath(dirname(@__FILE__), "Modes.pdf")
output_pdf = joinpath(dirname(@__FILE__), "Modes.pdf")
run(`pdfcrop $input_pdf $output_pdf`)



#################################àà
### spectrum
###################################


@load joinpath(directory_name, "out_ordered_k_th_0.jld2") ordered_k_th_0
@load joinpath(directory_name, "out_I_try0.jld2") I_try0
function read_D0_from_file(filename::String)
    for line in eachline(filename)
        if occursin("D0 =", line)
            m = match(r"D0\s*=\s*.*=\s*([\d\.eE\+\-]+)", line)
            if m !== nothing
                return parse(Float64, m.captures[1])
            end
        end
    end
    error("D0 not found in file.")
end

# D0 reading
filepath = joinpath(directory_name, "out_data.txt")
D0 = read_D0_from_file(filepath)    # 1.05 * gain clamping D0


function spectrum_plot(k_a, D0, ordered_k_th_0, I_try0, Re_k_min, Re_k_max)

    max_I = maximum(I_try0)

    # Crea il plot base con font ingranditi
    plt = plot(
        xlims=(k_a-1, k_a+1),
        xlabel=L"k \, (\mu \text{m}^{-1})",
        ylabel=L"\text{Intensity (arb.units)}",
        yticks=nothing,
        title = "",# L"(a)~D_0~\text{at 105% of the gain clamping value}",
        legend=false,
        size=(800, 500),
        guidefontsize=18,
        tickfontsize=16,
        titlefontsize=20,
    )

    # Estrai correttamente i tick dell’asse x
    xtick_vals, xtick_labels = Plots.xticks(plt)[1]

    # Aggiungi k_a se non presente
    if !any(x -> isapprox(x, k_a; atol=1e-8), xtick_vals)
        push!(xtick_vals, k_a)
        push!(xtick_labels, L"k_a")
    end

    # Riordina
    idx = sortperm(xtick_vals)
    xtick_vals = xtick_vals[idx]
    xtick_labels = xtick_labels[idx]

    # Applica xticks modificati
    plot!(xticks=(xtick_vals, xtick_labels))

    # Linee verticali blu scuro per lo spettro
    for i in 1:length(ordered_k_th_0)
        plot!([ordered_k_th_0[i], ordered_k_th_0[i]], [0, I_try0[i]], color=:gray, lw=2)
    end

    # Linea tratteggiata su k_a
    plot!([k_a, k_a], [0, max_I], color=:darkred, linestyle=:dash)

    # Scatter
    plot!(ordered_k_th_0, I_try0, seriestype=:scatter, color=:gray)

    # Salva e mostra
    savefig(plt, joinpath(@__DIR__, "spectrum_plot.png"))
    display(plt)
end


spectrum_plot(k_a,D0,ordered_k_th_0,I_try0, k_a-2.0,k_a+2.0)

#############################
## fibre int
##############################

k_th_0 = [all_paths_th_0[i][end] for i in eachindex(all_paths_th_0)]
k_th_0_and_D0 = [[k_th_0[i], all_D0s_th_0[i]] for i in eachindex(k_th_0)]
@load joinpath(directory_name, "out_n_sm_eigvec_gain.jld2") n_sm_eigvec_gain


function nij_g(k,D0)
    n = length(nij)
    nij_g = [0.0 + 0.0im for i in 1:n]
    for i in eachindex(nij)
        nij_g[i] = sqrt(nij[i]^2 + delta_pump[i] * D0*gamma_perp/(k - k_a + im*gamma_perp))
    end
    return nij_g
end

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


T_complex = T_mat_calc_3(k_th_0_and_D0, n_sm_eigvec_gain, k_a, gamma_perp, fibersLength, delta_pump, NEI, NEO, NE)
T_real_complete = real.(T_complex)
T_real = T_real_complete


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
 

ordered_modes, all_D0s_int, T_ordered_restricted, T_ordered_restricted_complex =  modes_sequence(T_real, T_complex, all_D0s_th_0)



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



    return I_try0
end

I_0 = spectrum_int(D0,all_D0s_int,k_th_0,T_real,ordered_modes)

function fibers_intens(x,I_0,k_th_0_and_D0,n_sm_eigvec_gain,ordered_modes,fibersLength, NEI, NEO; mod=1)
    i = mod
    
    res = [abs(sqrt(I_0[i])*Psi_mu(x, k_th_0_and_D0[ordered_modes[i]], n_sm_eigvec_gain[ordered_modes[i]], 
    fibersLength, NEI, NEO)[j])^2 for j in eachindex(fibers)]

    return res
end



pyplot()



function plot_fibers_intens(D0, k_th_0, T_real, ordered_modes, I_0, k_th_0_and_D0,
                             n_sm_eigvec_gain, fibersLength, NEI, NEO)

    p = plot(
        guidefontsize = 18,
        tickfontsize = 16,
        titlefontsize = 20,
        colorbar_titlefontsize = 18,
        size = (800, 600),
        xlims=(-43,43),
        ylims=(-43,43)
    )

    i_mod = argmax(I_0)

    max_fiber_length = maximum(fibersLength)
    base_num_points = 150

    max_intensity = -Inf
    min_intensity = Inf

    # Intensity interval
    for (i, _) in enumerate(fibersPoints)
        xs = range(0, stop=fibersLength[i], length=100)
        for x in xs
            fiber_intensity = fibers_intens(x, I_0, k_th_0_and_D0,
                                            n_sm_eigvec_gain, ordered_modes,
                                            fibersLength, NEI, NEO, mod=i_mod)[i]
            max_intensity = max(max_intensity, maximum(fiber_intensity))
            min_intensity = min(min_intensity, minimum(fiber_intensity))
        end
    end

    # Plot intensità su ogni fibra
    for (i, points) in enumerate(fibersPoints)
        num_samples = max(2, Int(round(base_num_points * (fibersLength[i] / max_fiber_length))))
        xs = range(0, stop=fibersLength[i], length=num_samples)

        x_start, y_start = points[1]
        x_end, y_end = points[2]

        x_vals = [(1 - t) * x_start + t * x_end for t in xs / fibersLength[i]]
        y_vals = [(1 - t) * y_start + t * y_end for t in xs / fibersLength[i]]

        intensities = [
            fibers_intens(x, I_0, k_th_0_and_D0,
                          n_sm_eigvec_gain, ordered_modes,
                          fibersLength, NEI, NEO, mod=i_mod)[i][1]
            for x in xs
        ]


        # (smooth line)
        plot!(x_vals, y_vals,
              line_z = intensities,
              c = :turbo,
              clim = (min_intensity, max_intensity),
              label = "",
              linewidth = 2)
    end

    xlabel!(p, L"x \; (\mu\mathrm{m})")
    ylabel!(p, L"y \; (\mu\mathrm{m})")

    plot!(p,
        aspect_ratio = :equal,
        colorbar_title = L"\text{Intensity (}E_c^2\text{ units)}"
    )

    output_path = joinpath(@__DIR__, "fibers_intensity_plot.png")
    savefig(p, output_path)

    display(p)
end




# Call the function
plot_fibers_intens(D0, k_th_0, T_real, ordered_modes, I_0, k_th_0_and_D0, n_sm_eigvec_gain, fibersLength, NEI, NEO)



