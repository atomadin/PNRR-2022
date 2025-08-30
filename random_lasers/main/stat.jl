using CSV
using DataFrames
using JLD2
using LinearAlgebra
using StatsBase, Distributions, HypothesisTests, Random, Statistics
using Plots
using Dates
using DelimitedFiles
using FilePaths, FileIO, Printf
using LaTeXStrings
using FLoops
using Distributions
using Statistics
using PyCall

pyplot()

# Check order function
function check_order(fibersIn)
    for k in 1:length(fibersIn)-1
        if fibersIn[k][1] > fibersIn[k+1][1] ||
           (fibersIn[k][1] == fibersIn[k+1][1] && fibersIn[k][2] >= fibersIn[k+1][2])
            error("The array is not correctly ordered: element $k ($(fibersIn[k])) precedes element $(k+1) ($(fibersIn[k+1]))")
        end
    end
end

# Base directory containing sample folders
function restr_sample_ind(base_dir::String)
    # Get all directories within the base directory
    sample_dirs = filter(isdir, readdir(base_dir, join=true))
    
    # Extract just the names of the directories
    dir_names = basename.(sample_dirs)
    
    # Filter to keep only directories named "sample_xxx"
    valid_dirs = dir_names[match.(r"^sample_\d+$", dir_names) .!== nothing]
    n_dir = length(valid_dirs)
    
    # Find indices of directories containing the file "out_I_try0.jld2"
    #valid_indices = findall(dir -> isfile(joinpath(base_dir, dir, "out_I_try0.jld2")), valid_dirs)
    valid_indices = findall(
        dir -> (isfile(joinpath(base_dir, dir, "out_I_try0.jld2")) &&
            isfile(joinpath(base_dir, dir, "out_D0.csv"))),
        valid_dirs)

    return valid_indices, n_dir
end

function get_full_path(n, base_dir)
    if n < 1 || n > 999
        error("n must be between 1 and 999")
    end
    dir_name = string("sample_", lpad(n, 3, '0'))
    full_path = joinpath(base_dir, dir_name)
    return full_path
end

base_dir = @__DIR__;
samples_with_data, n_dir = restr_sample_ind(base_dir);

# Function to read and process data for a given sample directory
function process_sample(sample_dir)
    # Construct paths for each file
    vertIn_file = joinpath(sample_dir, "InputData/vertIn.csv")
    vertOut_file = joinpath(sample_dir, "InputData/vertOut.csv")
    fibersIn_file = joinpath(sample_dir, "InputData/fibersIn.csv")
    fibersOut_file = joinpath(sample_dir, "InputData/fibersOut.csv")
    I_try0_file = joinpath(sample_dir, "out_I_try0.jld2")
    ordered_k_th_0 = joinpath(sample_dir, "out_ordered_k_th_0.jld2")
    #spectrum_file = joinpath(sample_dir, "out_spectrum.csv")
    D0_file = joinpath(sample_dir, "out_D0.csv")
    
    # Read data from CSV files
    vertIn = [row |> Vector for row in eachrow(CSV.read(vertIn_file, DataFrame, header=false))]
    vertOut = [row |> Vector for row in eachrow(CSV.read(vertOut_file, DataFrame, header=false))]
    fibersIn = [row |> Vector for row in eachrow(CSV.read(fibersIn_file, DataFrame, header=false))]
    fibersOut = [row |> Vector for row in eachrow(CSV.read(fibersOut_file, DataFrame, header=false))]
    @load I_try0_file I_try0
    @load ordered_k_th_0 ordered_k_th_0
    ####
    # spectrum = [row |> Vector for row in eachrow(CSV.read(spectrum_file, DataFrame, header=false))]
    # k_values = [spectrum[i][1] for i in eachindex(spectrum)]
    # I_try0 = [spectrum[i][2] for i in eachindex(spectrum)]
    ####
    D0 = CSV.read(D0_file, DataFrame, header=false)[1, 1]

    # Check order of fibers
    check_order(fibersIn)
    check_order(fibersOut)

    # Derive parameters
    NVIn = length(vertIn)
    NVOut = length(vertOut)
    vert = vcat(vertIn, vertOut)
    NV = NVIn + NVOut
    NEI = length(fibersIn)
    NEO = length(fibersOut)
    fibers = vcat(fibersIn, fibersOut)
    NE = length(fibers)
    fibersLength = [norm(vert[fibers[i][1]] - vert[fibers[i][2]]) for i in 1:NE]

    return NVIn, NVOut, vert, NV, NEI, NEO, fibers, NE, fibersLength, I_try0, #=k_values=# ordered_k_th_0, D0, vertIn, fibersIn
end


function process_all_samples(base_dir, samples_with_data, n_dir)
    # Initialize arrays
    NVIn_array = Vector{Union{Nothing, Int}}(nothing, n_dir)
    NVOut_array = Vector{Union{Nothing, Int}}(nothing, n_dir)
    vert_array = Vector{Union{Nothing, Vector{Vector{Float64}}}}(nothing, n_dir)
    NV_array = Vector{Union{Nothing, Int}}(nothing, n_dir)
    NEI_array = Vector{Union{Nothing, Int}}(nothing, n_dir)
    NEO_array = Vector{Union{Nothing, Int}}(nothing, n_dir)
    fibers_array = Vector{Union{Nothing, Vector{Vector{Int}}}}(nothing, n_dir)
    NE_array = Vector{Union{Nothing, Int}}(nothing, n_dir)
    fibersLength_array = Vector{Union{Nothing, Vector{Float64}}}(nothing, n_dir)
    I_try0_array = Vector{Union{Nothing, Any}}(nothing, n_dir)
    k_values_array = Vector{Union{Nothing, Any}}(nothing, n_dir)
    D0_array = Vector{Union{Nothing, Float64}}(nothing, n_dir)
    vertIn_array = Vector{Union{Nothing, Vector{Vector{Float64}}}}(nothing, n_dir)
    fibersIn_array = Vector{Union{Nothing, Vector{Vector{Int}}}}(nothing, n_dir)

    for index in samples_with_data
        sample_dir = get_full_path(index, base_dir)
        NVIn, NVOut, vert, NV, NEI, NEO, fibers, NE, fibersLength, I_try0, k_values, D0, vertIn, fibersIn = process_sample(sample_dir)

        # Check for valid data
        if !isnan(mean(I_try0)) && !any(x -> x < 0, I_try0)
            NVIn_array[index] = NVIn
            NVOut_array[index] = NVOut
            vert_array[index] = vert
            NV_array[index] = NV
            NEI_array[index] = NEI
            NEO_array[index] = NEO
            fibers_array[index] = fibers
            NE_array[index] = NE
            fibersLength_array[index] = fibersLength
            I_try0_array[index] = I_try0
            k_values_array[index] = k_values
            D0_array[index] = D0
            vertIn_array[index] = vertIn
            fibersIn_array[index] = fibersIn
        end
    end

    return NVIn_array, NVOut_array, vert_array, NV_array, NEI_array, NEO_array, fibers_array, NE_array, fibersLength_array, I_try0_array, k_values_array, D0_array, vertIn_array, fibersIn_array
end

NVIn_array, NVOut_array, vert_array, NV_array, NEI_array, NEO_array, fibers_array, NE_array, fibersLength_array, I_try0_array, k_values_array, D0_array, vertIn_array, fibersIn_array = process_all_samples(base_dir, samples_with_data, n_dir)

samples_with_data = []
for i in 1:n_dir
    if D0_array[i] != nothing
        push!(samples_with_data,i)
    end
end

@save joinpath(base_dir, "out_samples_with_data.jld2") samples_with_data



#####################
#### WD
####################


function compare_distributions(
    fibersLength::Vector{Float64};
    plt::Bool = true,
    save::Bool = false,
    plot_title::String = "(a)",
    remove_y_label::Bool = false
)
    mean_length = mean(fibersLength)
    normalized_lengths = sort(fibersLength ./ mean_length)

    λ = 1.0
    poisson_cdf = x -> cdf(Exponential(λ), x)
    wigner_dyson_cdf = x -> 1 - exp(-π * x^2 / 4)

    n = length(normalized_lengths)
    empirical_cdf = [i / n for i in 1:n]

    poisson_differences = [abs(empirical_cdf[i] - poisson_cdf(normalized_lengths[i])) for i in 1:n]
    wigner_dyson_differences = [abs(empirical_cdf[i] - wigner_dyson_cdf(normalized_lengths[i])) for i in 1:n]

    sup_poisson = maximum(poisson_differences)
    sup_wigner_dyson = maximum(wigner_dyson_differences)

    closer_distribution = if sup_poisson < sup_wigner_dyson
        "Poisson"
    else
        "Wigner-Dyson"
    end

    if plt
        x_vals = normalized_lengths
        poisson_cdf_vals = [poisson_cdf(x) for x in x_vals]
        wigner_dyson_cdf_vals = [wigner_dyson_cdf(x) for x in x_vals]

        p = plot(x_vals, empirical_cdf, label="Empirical CDF", lw=2, color=:blue)
        plot!(p, x_vals, poisson_cdf_vals, label="Poisson-like CDF", lw=2, linestyle=:dash, color=:red)
        plot!(p, x_vals, wigner_dyson_cdf_vals, label="Wigner-Dyson CDF", lw=2, linestyle=:dot, color=:green)

        plot!(p;
            xlabel="Normalized Length",
            ylabel=remove_y_label ? "" : "CDF",
            yformatter=remove_y_label ? :none : :auto,
            title=save ? plot_title : "Comparison of Empirical and Theoretical CDFs",
            legend=remove_y_label ? false : :bottomright
        )

        if save
            # Estrai contenuto tra parentesi, es. "a" o "2"
            m = match(r"\((\w)\)", plot_title)
            suffix = "x"  # default se parsing fallisce
            if m !== nothing
                s = m.captures[1]
                if occursin(r"^\d+$", s)
                    suffix = s
                elseif occursin(r"^[a-zA-Z]$", s)
                    suffix = string(Int(lowercase(s[1])) - Int('a') + 1)
                end
            end
            filename = joinpath(@__DIR__, "cdf_comparison_plot_$(suffix).pdf")
            savefig(p, filename)
        else
            display(p)
        end
    end

    return sup_poisson, sup_wigner_dyson, closer_distribution
end



# Evaluate distribution type for each sample
P_WD = Vector{Union{Nothing, Tuple{Float64, Float64, String}}}(nothing, n_dir)
for index in samples_with_data
    P_WD[index] = compare_distributions(fibersLength_array[index]; plt = false)
end

closer_distributions = Vector{Union{Nothing, String}}(nothing, n_dir)
for i in eachindex(P_WD)
    if P_WD[i] != nothing
        closer_distributions[i] = P_WD[i][3]
    end
end

poisson_indices = findall(==("Poisson"), closer_distributions)
wigner_dyson_indices = findall(==("Wigner-Dyson"), closer_distributions)
poisson_count = length(poisson_indices)
wigner_dyson_count = length(wigner_dyson_indices)


poisson_indices
wigner_dyson_indices
fibersLength_array[poisson_indices[1]]
compare_distributions(fibersLength_array[poisson_indices[1]])
compare_distributions(fibersLength_array[wigner_dyson_indices[1]])

wigner_dyson_indices[2]
compare_distributions(fibersLength_array[445], save=true)
compare_distributions(fibersLength_array[wigner_dyson_indices[2]], save=true, remove_y_label=true, plot_title="(b)")


############################################
# IPR
############################################

p_array = Vector{Union{Nothing, Vector{Float64}}}(undef, length(I_try0_array))
for i in eachindex(p_array)
    if I_try0_array[i] == nothing
        p_array[i] = nothing
    else
        p_array[i] = I_try0_array[i] ./ sum(I_try0_array[i])
    end
end

IPR_array = Vector{Union{Nothing, Float64}}(undef, length(I_try0_array))
for i in eachindex(p_array)
    if p_array[i] == nothing
        IPR_array[i] = nothing
    else
        # IPR = 1 / sum( p_j^2 )
        IPR_array[i] = 1/sum(p_array[i][j]^2 for j in eachindex(p_array[i]))
    end
end

IPR_indx = []
for i in eachindex(IPR_array)
    if IPR_array[i] != nothing && !isnan(IPR_array[i])
        push!(IPR_indx, i)
    end
end

IPR_indx = sort(union(IPR_indx, poisson_indices, wigner_dyson_indices))

filtered_IPR_array = IPR_array[IPR_indx]
indices_range = 1:length(filtered_IPR_array)

# Find indices in IPR_array that are NaN, ignoring Nothing values
nan_positions = findall(x -> x !== nothing && isnan(x), IPR_array)
maximum(I_try0_array[106])
maximum(I_try0_array[492])


# Split IPR by distribution
IPR_Poisson = [IPR_array[i] for i in poisson_indices if IPR_array[i] !== nothing]
IPR_WD      = [IPR_array[i] for i in wigner_dyson_indices if IPR_array[i] !== nothing]


### save IPR in samples
function save_ipr_to_file(IPR_array, base_dir)
    for i in 1:length(IPR_array)
        if IPR_array[i] !== nothing
            dir_name = string("sample_", lpad(i, 3, '0'))
            file_path = joinpath(base_dir, dir_name, "IPR_value.txt")
            open(file_path, "w") do file
                write(file, string(IPR_array[i], "\n"))
            end
        end
    end
end

save_ipr_to_file(IPR_array, base_dir)


###################################
# Tests & p-values
###################################

function test_poisson_distribution(fibersLength::Vector{Float64})
    mean_length = mean(fibersLength)
    # Small random shift to avoid ties in KS
    normalized_lengths = sort(fibersLength ./ mean_length .+ 1e-6 * randn(length(fibersLength)))
    d = Exponential(1.0)  
    test = ExactOneSampleKSTest(normalized_lengths, d)
    return pvalue(test)
end

wigner_dyson_cdf = x -> 1 - exp(-π * x^2 / 4)
struct WignerDysonDist <: ContinuousUnivariateDistribution end
Distributions.cdf(::WignerDysonDist, x::Real) = wigner_dyson_cdf(x)

function test_wigner_dyson_distribution(fibersLength::Vector{Float64})
    mean_length = mean(fibersLength)
    # Small random shift to avoid ties in KS
    normalized_lengths = sort(fibersLength ./ mean_length .+ 1e-6 * randn(length(fibersLength)))
    test = ExactOneSampleKSTest(normalized_lengths, WignerDysonDist())
    return pvalue(test)
end

p_values_array = Vector{Union{Nothing, Float64}}(nothing, n_dir)
for i in 1:n_dir
    if closer_distributions[i] == "Poisson"
        p_values_array[i] = test_poisson_distribution(fibersLength_array[i])
    elseif closer_distributions[i] == "Wigner-Dyson"
        p_values_array[i] = test_wigner_dyson_distribution(fibersLength_array[i])
    else
        p_values_array[i] = nothing
    end
end

# For reference, p-values specifically tested for the assigned distribution
p_value_poisson_for_poisson = map(i -> test_poisson_distribution(fibersLength_array[i]), poisson_indices)
p_value_WD_for_wd = map(i -> test_wigner_dyson_distribution(fibersLength_array[i]), wigner_dyson_indices)

###########################
# Write out "stat_P_WD.txt"
###########################

function print_sta_file(closer_distributions, p_values_array, base_dir)
    for i in 1:length(closer_distributions)
        dist = closer_distributions[i]
        if dist === nothing
            continue
        end
        dir_name = string("sample_", lpad(i, 3, '0'))
        file_path = joinpath(base_dir, dir_name, "stat_P_WD.txt")

        pval_str = p_values_array[i] === nothing ? "N/A" : string(p_values_array[i])

        open(file_path, "w") do file
            if dist == "Poisson"
                write(file, "P\n$pval_str\n")
            elseif dist == "Wigner-Dyson"
                write(file, "WD\n$pval_str\n")
            end
        end
    end
end

print_sta_file(closer_distributions, p_values_array, base_dir)


# ---------------------------------------------------
# Mean value + Detuning of the max intensity
# ---------------------------------------------------
function save_weighted_spectrum(I_try0_array, k_values_array, base_dir; k_a=10.68)
    """
    For each sample i:
      - Compute the weighted mean of the spectrum and its difference from k_a
      - Find the index k_max where the intensity is maximal
        and compute (k_max - k_a).
      - Save them to a file spectrum_mean.txt 
        (4 lines: weighted mean, difference, k_max, difference(k_max - k_a))
    """
    
    for i in 1:length(I_try0_array)
        if I_try0_array[i] !== nothing
            I = I_try0_array[i]
            # Define k_array = 1,2,3,...,N
            k_array = k_values_array[i]
            
            # WeightedMean = sum( k_j * I[j] ) / sum(I[j])
            total_intensity = sum(I)
            weighted_mean = sum(k_array[j]*I[j] for j in eachindex(I)) / total_intensity
            difference_wm = weighted_mean - k_a

            # Find the index of maximum intensity
            # (in Julia, argmax returns the index of the first maximum).
            k_max_ind = argmax(I)
            k_max = k_array[k_max_ind]
            difference_max = k_max - k_a

            # Write all to file
            dir_name = string("sample_", lpad(i, 3, '0'))
            file_path = joinpath(base_dir, dir_name, "spectrum_mean.txt")
            open(file_path, "w") do f
                # 1) Weighted mean
                # 2) Weighted mean difference (detuning)
                # 3) k_max
                # 4) difference of k_max from k_a
                write(f, "$weighted_mean\n$difference_wm\n$k_max\n$difference_max\n")
            end
        end
    end
end


save_weighted_spectrum(I_try0_array, k_values_array, base_dir; k_a=10.68)


##################################################
# Read weighted mean, difference from k_a, 
#   k_max and difference for k_max
##################################################
function read_spectrum_means(base_dir, n_dir)
    """
    Reads 'spectrum_mean.txt' from each sample directory.
    Returns four arrays of length n_dir:
      weighted_mean_array[i], difference_array[i],
      kmax_array[i], difference_max_array[i]
    or 'nothing' if not available.
    """
    weighted_mean_array     = Vector{Union{Nothing,Float64}}(undef, n_dir)
    difference_array        = Vector{Union{Nothing,Float64}}(undef, n_dir)
    kmax_array              = Vector{Union{Nothing,Float64}}(undef, n_dir)
    difference_max_array    = Vector{Union{Nothing,Float64}}(undef, n_dir)

    for i in 1:n_dir
        dir_name = string("sample_", lpad(i, 3, '0'))
        file_path = joinpath(base_dir, dir_name, "spectrum_mean.txt")
        if isfile(file_path)
            vals = readlines(file_path)
            if length(vals) >= 4
                # 1) weighted mean
                # 2) difference from k_a
                # 3) k_max
                # 4) difference from k_a for k_max
                wm = parse(Float64, vals[1])
                diff_wm = parse(Float64, vals[2])
                km = parse(Float64, vals[3])
                diff_km = parse(Float64, vals[4])
                weighted_mean_array[i]  = wm
                difference_array[i]     = diff_wm
                kmax_array[i]           = km
                difference_max_array[i] = diff_km
            else
                weighted_mean_array[i]  = nothing
                difference_array[i]     = nothing
                kmax_array[i]           = nothing
                difference_max_array[i] = nothing
            end
        else
            weighted_mean_array[i]  = nothing
            difference_array[i]     = nothing
            kmax_array[i]           = nothing
            difference_max_array[i] = nothing
        end
    end

    return weighted_mean_array, difference_array, kmax_array, difference_max_array
end

weighted_mean_array, difference_array, kmax_array, difference_max_array = read_spectrum_means(base_dir, n_dir)

weighted_mean_Poisson = [weighted_mean_array[i] for i in poisson_indices if weighted_mean_array[i] !== nothing]
weighted_mean_WD      = [weighted_mean_array[i] for i in wigner_dyson_indices if weighted_mean_array[i] !== nothing]

difference_Poisson      = [difference_array[i] for i in poisson_indices if difference_array[i] !== nothing]
difference_WD           = [difference_array[i] for i in wigner_dyson_indices if difference_array[i] !== nothing]

kmax_Poisson           = [kmax_array[i] for i in poisson_indices if kmax_array[i] !== nothing]
kmax_WD                = [kmax_array[i] for i in wigner_dyson_indices if kmax_array[i] !== nothing]

difference_max_Poisson = [difference_max_array[i] for i in poisson_indices if difference_max_array[i] !== nothing]
difference_max_WD      = [difference_max_array[i] for i in wigner_dyson_indices if difference_max_array[i] !== nothing]



##################################################
# PLOTS: 
#   1) IPR vs p-value
#   2) Weighted mean vs p-value
#   3) (weighted mean - k_a) vs p-value
#   4) (k_max - k_a) vs p-value [the new detuning of max]
##################################################

# We'll build arrays for each distribution type to keep them separate/colors.
function get_scatter_data(IPR_array, p_values_array, 
                          weighted_mean_array, difference_array,
                          kmax_array, difference_max_array,
                          dist_label)
    """
    Returns vectors of 
        x_ipr, y_p, x_wm, x_diff, x_kmax, x_diffmax
    for samples labeled with dist_label 
    (either "Poisson" or "Wigner-Dyson").
    """
    idx = findall(i -> closer_distributions[i] == dist_label &&
                        IPR_array[i] != nothing &&
                        p_values_array[i] != nothing &&
                        weighted_mean_array[i] != nothing &&
                        difference_array[i] != nothing &&
                        kmax_array[i] != nothing &&
                        difference_max_array[i] != nothing,
                  1:length(IPR_array))

    x_ipr       = [IPR_array[i] for i in idx]
    y_p         = [p_values_array[i] for i in idx]
    x_wm        = [weighted_mean_array[i] for i in idx]
    x_diff      = [difference_array[i] for i in idx]
    x_kmax      = [kmax_array[i] for i in idx]
    x_diffmax   = [difference_max_array[i] for i in idx]

    return x_ipr, y_p, x_wm, x_diff, x_kmax, x_diffmax
end

x_ipr_p, y_p_p, x_wm_p, x_diff_p, x_kmax_p, x_diffmax_p = get_scatter_data(
    IPR_array, p_values_array, weighted_mean_array, difference_array, 
    kmax_array, difference_max_array, "Poisson"
)
x_ipr_wd, y_p_wd, x_wm_wd, x_diff_wd, x_kmax_wd, x_diffmax_wd = get_scatter_data(
    IPR_array, p_values_array, weighted_mean_array, difference_array, 
    kmax_array, difference_max_array, "Wigner-Dyson"
)



################ PLOTs
#############################

# 1) IPR


# Function to compute the number of bins using Sturges' rule
function sturges_rule(n::Int)
    #Sturges' Rule number of bins
    return ceil(Int, log2(n) + 1)
end




function plot_histogram(data, name, label, n_bins; prefix="", save=true, xlabel_on=true, ylabel_on=true)
    if isempty(data)
        println("Skipping histogram for $name ($label) due to empty data.")
        return
    end

    title_str = "$prefix Histogram of $label for $name"

    xlabel_str = xlabel_on ? "$label" : ""
    ylabel_str = ylabel_on ? L"\mathrm{Frequency}" : ""

    histogram(data,
              bins = n_bins,
              xlabel = xlabel_str,
              ylabel = ylabel_str,
              fillcolor = :grey,
              #title = title_str,
              legend = false,
              guidefontsize = 32,
              tickfontsize = 30,
              titlefontsize = 40,
              size = (800, 600)
    )

    if save
        filename = joinpath(@__DIR__, "plot_paper", "$(label)_$(name)_histogram.png")
        savefig(filename)
        println("Saved histogram to $filename")
    end

    display(current())  # mostra il grafico
end


n_bins = sturges_rule(length(IPR_Poisson))
n_bins = sturges_rule(length(IPR_WD))

plot_histogram(IPR_Poisson, "Poisson", "IPR", 10, prefix="(a)")
plot_histogram(IPR_WD, "Wigner-Dyson", "IPR", 10, prefix="(b)")

#mean(IPR_Poisson)
mean(filter(!isnan, IPR_Poisson))
std(filter(!isnan, IPR_Poisson))
mean(filter(!isnan, IPR_WD))
std(filter(!isnan, IPR_WD))


x_p_dist = [P_WD[IPR_indx[i]][1] for i in eachindex(IPR_indx)]
y_IPR = [IPR_array[IPR_indx[i]] for i in eachindex(IPR_indx)]
scatter(x_p_dist, y_IPR, markersize=4, markercolor=:blue, label="IPR data", xlabel="distance from Poisson", ylabel="IPR", title="IPR Scatter Plot")

x_wd_dist = [P_WD[IPR_indx[i]][2] for i in eachindex(IPR_indx)]
scatter(x_wd_dist, y_IPR, markersize=4, markercolor=:blue, label="IPR data", xlabel="distance from WD", ylabel="IPR", title="IPR Scatter Plot")



function plot_histogram_with_exponential(data, name, label, n_bins; prefix="", save=true)
    data = filter(!isnan, data)  # remove NaNs
    
    if isempty(data)
        println("Skipping histogram for $name ($label) due to empty or invalid data.")
        return
    end
    
    #lamb_est = mean(data) # 3.787
    lamb_est = 3.5

    title_str = "$prefix Histogram of $label for $name"

    # Plot histogram normalized as probability density
    hist = histogram(data,
                     bins = n_bins,
                     normalize = true,
                     xlabel = "$label $name",
                     ylabel = L"\mathrm{Relative\ Frequency}",
                     fillalpha = 0.5,
                     fillcolor = :black,
                     linecolor = :black,
                     #title = title_str,
                     legend = false,
                     #label = "Data",
                     guidefontsize = 18,
                     tickfontsize = 16,
                     titlefontsize = 20,
                     legendfontsize = 16,
                     #size = (800, 600)
                     size = (800, 400) 
    )

    # Fit the Exponential distribution
    exp_dist = Exponential(lamb_est)

    min_val, max_val = minimum(data), maximum(data)
    if !isfinite(min_val) || !isfinite(max_val)
        println("Cannot plot Exponential overlay: data range contains non-finite values.")
        return
    end

    x_vals = range(min_val, stop = max_val, length = 500)
    y_vals = pdf.(exp_dist, x_vals)

    # Plot PDF
    plot!(x_vals, y_vals,
          linewidth = 3,
          linecolor = :black,
          linestyle = :dash)
          #label = "Exponential(λ=0.269)")

    # Extend y-axis limits to avoid clipping
    max_y = maximum(vcat(y_vals, hist.series_list[1][:y]))
    ylims!(0, 0.21)

    if save
        filename = joinpath(@__DIR__, "plot_paper", "$(label)_$(name)_histogram_exponential.pdf")
        savefig(filename)
        println("Saved histogram with Exponential overlay to $filename")
    end

    display(current())

    return lamb_est
end



plot_histogram_with_exponential(IPR_Poisson, "", "", 10, prefix="")




function plot_histogram_with_WD(data, name, label, n_bins; prefix="", save=true)
    # Rimuovi eventuali NaN
    data = filter(!isnan, data)

    if isempty(data)
        println("Skipping histogram for $name ($label) due to empty or invalid data.")
        return
    end

    # Calcola stima di lambda per Wigner-Dyson
    function estimate_lambda_wigner_dyson(data)
        n = length(data)
        s2 = sum(x^2 for x in data)
        return sqrt(4n / (π * s2))
    end

    # Definisci PDF Wigner-Dyson
    function wigner_dyson_pdf(x::Real, lambda::Real)
        if x < 0
            return 0.0
        end
        return (π * lambda^2 * x / 2) * exp(-π * lambda^2 * x^2 / 4)
    end

    # Stima lambda
    # stima_lambda = estimate_lambda_wigner_dyson(data) # 0.25528
    stima_lambda = 0.31


    title_str = "$prefix Histogram of $label for $name"

    # Istogramma normalizzato
    hist = histogram(data,
                     bins = n_bins,
                     normalize = true,
                     xlabel = "$label $name",
                     ylabel = L"\mathrm{Relative\ Frequency}",
                     fillalpha = 0.5,
                     fillcolor = :black,
                     linecolor = :black,
                     title = "",#title_str,
                     legend = false,
                     label = "Data",
                     guidefontsize = 18,
                     tickfontsize = 16,
                     titlefontsize = 20,
                     legendfontsize = 16,
                     #size = (800, 600)
                     size = (800, 400) 
    )

    min_val, max_val = minimum(data), maximum(data)
    if !isfinite(min_val) || !isfinite(max_val)
        println("Cannot plot Wigner-Dyson overlay: data range contains non-finite values.")
        return
    end

    # Costruisci asse x per il fit
    x_start = max(min_val, eps())
    x_vals = range(x_start, stop = max_val, length = 500)
    y_vals = wigner_dyson_pdf.(x_vals, stima_lambda)

    # Sovrapponi la curva Wigner-Dyson
    plot!(x_vals, y_vals,
          linewidth = 3,
          linecolor = :black,
          linestyle = :dashdot,
          label = @sprintf("Wigner-Dyson(λ = %.3f)", stima_lambda))

    # Estendi asse y per evitare tagli
    max_y = maximum(vcat(y_vals, hist.series_list[1][:y]))
    ylims!(0, 0.25)

    if save
        filename = joinpath(@__DIR__, "plot_paper", "$(label)_$(name)_histogram_WD.pdf")
        savefig(filename)
        println("Saved histogram with Wigner-Dyson overlay to $filename")
    end

    display(current())

    return stima_lambda

end




plot_histogram_with_WD(IPR_WD, "", "Modal Intensity Distribution IPR", 10, prefix="")

# 2) Weighted Mean


plot_histogram(weighted_mean_Poisson, "Poisson", "Weighted Mean Spectrum", 10, save=false)
plot_histogram(weighted_mean_WD, "Wigner-Dyson", "Weighted Mean Spectrum", 10, save=false)

mean(filter(!isnan, weighted_mean_Poisson))
std(filter(!isnan, weighted_mean_Poisson))
mean(filter(!isnan, weighted_mean_WD))
std(filter(!isnan, weighted_mean_WD))


y_weighted_mean = [weighted_mean_array[IPR_indx[i]] for i in eachindex(IPR_indx)]
scatter(x_p_dist, y_weighted_mean, markersize=4, markercolor=:blue, label="Mean spectrum data", xlabel="distance from Poisson", ylabel="mean spectrum", title="Mean Scatter Plot")
scatter(x_wd_dist, y_weighted_mean, markersize=4, markercolor=:blue, label="Mean spectrum data", xlabel="distance from WD", ylabel="mean spectrum", title="Mean Scatter Plot")


# 3) (weighted mean - k_a) (detuning)


n_bins = sturges_rule(length(difference_Poisson))
n_bins = sturges_rule(length(difference_WD))


plot_histogram(difference_Poisson, "Poisson", "Spectrum Mean Detuning",10, prefix="(a)")
plot_histogram(difference_WD, "Wigner-Dyson", "Spectrum Mean Detuning",10, prefix="(b)", ylabel_on=false)


m_mean_det_Poiss = mean(filter(!isnan, difference_Poisson))
std(filter(!isnan, difference_Poisson))
m_mean_det_WD = mean(filter(!isnan, difference_WD))
std(filter(!isnan, difference_WD))

(m_mean_det_WD - m_mean_det_Poiss)/ m_mean_det_Poiss * 100




y_detuning_mean = [difference_array[IPR_indx[i]] for i in eachindex(IPR_indx)]
scatter(x_p_dist, y_detuning_mean, markersize=4, markercolor=:blue, label="Mean Detuning data", xlabel="distance from Poisson", ylabel="IPR", title="Mean Detuning Scatter Plot")
scatter(x_wd_dist, y_detuning_mean, markersize=4, markercolor=:blue, label="IPR data", xlabel="distance from WD", ylabel="IPR", title="Mean Detuning Scatter Plot")


# 4) (k_max - k_a) max detuning


n_bins = sturges_rule(length(difference_max_Poisson))
n_bins = sturges_rule(length(difference_max_WD))


plot_histogram(difference_max_Poisson, "Poisson", "Spectrum Max Detuning", 10, prefix="(c)")
plot_histogram(difference_max_WD, "Wigner-Dyson", "Spectrum Max Detuning", 10, prefix="(d)", ylabel_on = false)

m_max_det_Poiss = mean(filter(!isnan, difference_max_Poisson))
std(filter(!isnan, difference_max_Poisson))
m_max_det_WD = mean(filter(!isnan, difference_max_WD))
std(filter(!isnan, difference_max_WD))

(m_max_det_WD - m_max_det_Poiss)/ m_max_det_Poiss * 100



y_detuning_max = [difference_max_array[IPR_indx[i]] for i in eachindex(IPR_indx)]
scatter(x_p_dist, y_detuning_max, markersize=4, markercolor=:blue, label="Mean Detuning data", xlabel="distance from Poisson", ylabel="IPR", title="Mean Detuning Scatter Plot")
scatter(x_wd_dist, y_detuning_max, markersize=4, markercolor=:blue, label="IPR data", xlabel="distance from WD", ylabel="IPR", title="Mean Detuning Scatter Plot")


#####################################################
#####################################################


# TRAILS


function trail_add_level(input, lati)
    # STEP 1: Seleziona da edges solo le coppie non ordinate **diverse** da quelle presenti in input consecutive
    # Costruiamo le coppie consecutive da input
    input_pairs = [[input[i], input[i+1]] for i in 1:length(input)-1]
    # Consideriamo anche la versione non ordinata
    input_pairs_unordered = Set([sort(p) for p in input_pairs])

    # Costruiamo edges_restr rimuovendo le coppie che corrispondono (non ordinate)
    edges_restr = [e for e in lati if sort(e) ∉ input_pairs_unordered]

    #println("edges_restr = ", edges_restr)

    # STEP 2: Controlla se gli estremi di input (inizio o fine) coincidono con uno degli estremi delle edges_restr
    first_node = input[1]
    last_node = input[end]
    new_paths = []

    for e in edges_restr
        if first_node in e
            # Prepend all’altro nodo
            new_node = e[1] == first_node ? e[2] : e[1]
            push!(new_paths, [new_node; input])
        end
        if last_node in e
            # Append all’altro nodo
            new_node = e[1] == last_node ? e[2] : e[1]
            push!(new_paths, vcat(input, new_node))
        end
    end

    return new_paths

end





function remove_redund_unordered(all_trails)
    function trail_to_edge_set(trail)
        edges = Set([join(sort([trail[i], trail[i+1]]), "-") for i in 1:length(trail)-1])
        return edges
    end

    seen = Set{Set{String}}()
    unique_trails = Vector{Vector{Int}}()

    for trail in all_trails
        edge_set = trail_to_edge_set(trail)
        if !(edge_set in seen)
            push!(unique_trails, trail)
            push!(seen, edge_set)
        end
    end

    return unique_trails
end




function trails_calc(edges, n_max)

    n_cycles = n_max - 1
    # trails contiene tutti gli edges
    grouped_trails = [deepcopy(edges)]

    for i in 1:n_cycles
        # considero i trails con 3 edges
        all_trails = vcat(trail_add_level.(grouped_trails[i],Ref(edges))...)
        # rimuovo ridondanti
        trails_uniq = [remove_redund_unordered(all_trails)]
        # aggiungo a trails
        append!(grouped_trails,trails_uniq)
    end

    trails = vcat(grouped_trails...)
    #trails = grouped_trails

    return trails 

end

function group_trails_by_unique_nodes(trails::Vector{Vector{Int}}, max_vert::Int; skip_one::Bool = true)
    groups = Dict{Int, Vector{Vector{Int}}}()

    for trail in trails
        n_unique = length(unique(trail))
        push!(get!(groups, n_unique, Vector{Vector{Int}}()), trail)
    end

    max_len = maximum(keys(groups))
    result_p = [get(groups, i, []) for i in 1:max_len]

    if max_len > max_vert
        result_pp = result_p[1:max_vert]
    else
        result_pp = result_p
    end

    if skip_one
        result =  result_pp[2:end]
    else
        result =  result_p[1:end]
    end

    
    return result
end


function remove_loops(fibers::Vector{Vector{Int}})
    return [edge for edge in fibers if edge[1] != edge[2]]
end


fibre_filtrate_array = Vector{Union{Vector{Vector{Int64}}, Nothing}}(undef, 999)
fill!(fibre_filtrate_array, nothing)
@floop for i in samples_with_data
    fibre_filtrate_array[i] = remove_loops(fibers_array[i])
end


## Saved TRIALS
tutti_trails_array = Vector{Union{Vector{Vector{Int64}}, Nothing}}(undef, 999)
fill!(tutti_trails_array, nothing)
@timed begin
    #Threads.@threads 
    for i in samples_with_data
        tutti_trails_array[i] = trails_calc(fibre_filtrate_array[i], 9)
    end
end

@save joinpath(@__DIR__,  "tutti_trails_array.jld2") tutti_trails_array
# @load joinpath(@__DIR__,  "tutti_trails_array.jld2") tutti_trails_array

trails_by_nodes_array = Vector{Union{Vector{Vector{Vector{Int}}}, Nothing}}(undef, 999)
fill!(trails_by_nodes_array, nothing)
@timed begin
    #Threads.@threads 
    for i in samples_with_data
        trails_by_nodes_array[i] = group_trails_by_unique_nodes(tutti_trails_array[i],7)
    end
end

@save joinpath(@__DIR__,  "trails_by_nodes_array.jld2") trails_by_nodes_array
# @load joinpath(@__DIR__,  "trails_by_nodes_array.jld2") trails_by_nodes_array

function trail_length(trail, positions)
    total = 0.0
    for i in 1:length(trail)-1
        p1 = positions[trail[i]]
        p2 = positions[trail[i+1]]
        total += norm([p2[1] - p1[1], p2[2] - p1[2]])
    end
    return total
end


lunghezze_per_nodi_array = Vector{Union{Vector{Vector{Float64}}, Nothing}}(fill(nothing, 999))
for j in samples_with_data
    lunghezze_per_nodi_array[j] = [map(trail -> trail_length(trail, vert_array[j]), trails_by_nodes_array[j][i]) for i in eachindex(trails_by_nodes_array[j])]
end



histogram(lunghezze_per_nodi_array[1][4]; bins=20, xlabel="Lunghezza", ylabel="Frequenza", title="Distribuzione lunghezze trail")


# mean value and standard deviation
media_lunghezze_per_nodi_array = Vector{Union{Vector{Float64}, Nothing}}(undef, 999)
deviazione_std_per_nodo_array = Vector{Union{Vector{Float64}, Nothing}}(undef, 999)

#@floop 
for i in samples_with_data
    media_lunghezze_per_nodi_array[i] = mean.(lunghezze_per_nodi_array[i])
end
#@floop 
for i in samples_with_data
    deviazione_std_per_nodo_array[i] = std.(lunghezze_per_nodi_array[i])
end


#############
# Isolation Forest: Anomaly Detection Using Decision Trees

# === 1. Isolation Forest ===
@pyimport sklearn.ensemble as sk_ensemble


maximum(poisson_indices)

features = Vector{Vector{Float64}}()
valid_indices = Int[]

for i in samples_with_data
    if i < 868
        media = media_lunghezze_per_nodi_array[i]
        deviazione = deviazione_std_per_nodo_array[i]
        if media !== nothing && deviazione !== nothing
            push!(features, vcat(media, deviazione))
            push!(valid_indices, i)
        end
    end
end

X = hcat(features...)'  # (n_samples × 12)

clf = sk_ensemble.IsolationForest(n_estimators=1000, contamination="auto")#, random_state=42)
clf.fit(X)
scores = clf.decision_function(X)
threshold = quantile(scores, 0.10)  # 10% considered anomalous

#######
samples_with_data_restr = [i for i in samples_with_data if i < 868]

# check
valid_indices == samples_with_data_restr

poisson_indices_restr =  [i for i in poisson_indices if i < 868]
wigner_dyson_indices_restr = [i for i in wigner_dyson_indices if i < 868]


scores_complete = Vector{Union{Nothing, Float64}}(fill(nothing, 999))
scores_complete[valid_indices] .= scores

mean(scores_complete[poisson_indices_restr])
mean(scores_complete[wigner_dyson_indices_restr])

(mean(scores_complete[poisson_indices_restr])- mean(scores_complete[wigner_dyson_indices_restr]))/mean(scores_complete[poisson_indices_restr])*100


######

data_spectrum = IPR_array[valid_indices]

function plot_anomaly_vs_spectrum_data(data_spectrum::Vector{Union{Nothing, Float64}}, name::String;
    scores=scores::Vector{Float64}, threshold=threshold::Float64, valid_indices=valid_indices::Vector{Int64},
    save_fig=true::Bool)

    data_spectrum_partial = data_spectrum[valid_indices]
    valid_idx = [i for i in eachindex(data_spectrum_partial) if !isnan(data_spectrum_partial[i])]

    data_spectrum_clean = data_spectrum_partial[valid_idx]
    scores_clean = scores[valid_idx]

    # Determina chi è anomalo: is_outlier vale true se è oltre la soglia
    is_outlier = scores_clean .< threshold
    data_spectrum_normali = data_spectrum_clean[.!is_outlier]
    data_spectrum_anomali = data_spectrum_clean[is_outlier]

    mu_norm = mean(data_spectrum_normali)
    mu_anom = mean(data_spectrum_anomali)

    scatter(data_spectrum_normali, scores_clean[.!is_outlier], 
            label = "Normal", color = :gray,
            xlabel = "$(name)", ylabel = "Anomaly Score",
            legend = :topright,
            dpi=300)

    scatter!(data_spectrum_anomali, scores_clean[is_outlier], 
            label = "Anomalous", color = :red,markershape = :cross )

    display(current())
    
    if save_fig
        filename = joinpath(@__DIR__, "plot_paper", "$(name)_anomalies_histogram.png")
        savefig(filename)
    end

    return mu_norm, mu_anom

end


mu_norm_IPR, mu_anom_IPR = plot_anomaly_vs_spectrum_data(IPR_array,"IPR")

(mu_norm_IPR- mu_anom_IPR)/mu_norm_IPR*100

########

mu_norm_IWD, mu_anom_IWD = plot_anomaly_vs_spectrum_data(difference_array,"IWD")

(mu_norm_IWD - mu_anom_IWD)/mu_norm_IWD*100
mu_norm_IWD
mu_anom_IWD


#############

mu_norm_MID, mu_anom_MID = plot_anomaly_vs_spectrum_data(difference_max_array,"MID")

(mu_norm_MID - mu_anom_MID)/mu_norm_MID*100
mu_norm_IWD
mu_anom_IWD
