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


using CSV, DataFrames
using Printf
using Optim
using FLoops
using Measures
using Plots
using BlackBoxOptim
using Distributed

#######################################################################################################################
#######################################################################################################################


function getFermiDiracHoles(Energy::Float64, Mu::Float64, KT::Float64)
    
    KT_CUTOFF = 5.0    # Thermal energy cutoff

    if (Mu - Energy) >= (KT_CUTOFF * KT)
        return 0.0
    else
        return 1.0 / (exp((Mu - Energy) / KT) + 1.0)
    end
end



#######################################################################################################################
#######################################################################################################################

function N2D_calc(Ef::Float64,T::Float64)
    
    NumPopStates= 4.0
    NumEmpStates=24.0
    d3K=0.0000001
    Aperp=117.5

    EDOUBLET = 5.0                   # Doublet energy
    KBOLTZMANN = 0.086173423         # Boltzmann constant (meV/K)
    EnArraySize = 500                # Energy array size
    ETrMax = 1000.0                  # Maximum energy transfer (meV)
    DIM = 1e14 / pi / 137.0 / 4.0    
    # 10^14/PI*costante di struttura fine / indice di rifrazione (germanio) con dipolo da codice e energie in meV  e d3k in A^-3 =#


    KT = KBOLTZMANN.*T
    d3K=d3K.*1.0 # here I consider the degeneracy of 1 due to the grid3dfull full positive and negative values for a1 a2 and z 
    #  d3K=d3K*8.d0 ! here I consider the degeneracy of 8 due to the grid3d not full with positive values only for a1 a2 and z 


    # Define the file path
    file_path_input = joinpath(@__DIR__, "InputData", "QWIP_VSx80_APTprofile_highGe_CMP3D.dat")
    # Attempt to read the file without a header
    data_df = CSV.read(file_path_input, DataFrame; delim=' ', ignorerepeated=true, header=false)
    # Assign default names if necessary
    col_names = [:k_cnt, :kx_cur, :ky_cur, :fill_cnt, :other_cnt, :EFill, :EOther, :PXX2, :PYY2, :PZZ2]
    rename!(data_df, col_names[1:ncol(data_df)])
    # Access arrays directly
    # k_cnt_arr = data_df[:, :k_cnt]
    # kx_cur_arr = data_df[:, :kx_cur]
    # ky_cur_arr = data_df[:, :ky_cur]
    # fill_cnt_arr = data_df[:, :fill_cnt]
    other_cnt_arr = data_df[:, :other_cnt]
    EFill_arr = data_df[:, :EFill] .* 1000.0
    EOther_arr = data_df[:, :EOther] .* 1000.0
    # PXX2_arr = data_df[:, :PXX2]
    # PYY2_arr = data_df[:, :PYY2]
    # PZZ2_arr = data_df[:, :PZZ2]


    absEtr_arr = abs.(EFill_arr .- EOther_arr)
    # Find indices of elements greater than the threshold
    indices = findall(x -> x > ETrMax*1000.0, absEtr_arr)
    if !isempty(indices)
        @error("Error abs(EFill - EOther) > ETrMax.")
    end

    ### N2D calc
    # Calculate the target value
    target_value = NumEmpStates + NumPopStates
    # Find indices where the condition is satisfied
    indices = findall(x -> x == target_value, other_cnt_arr)
    N2D = d3K / (8.0 * π * π * π) * 1.0e24 * 1.0e-8 * Aperp * sum(getFermiDiracHoles(EFill_arr[indices[i]], Ef, KT) for i in eachindex(indices))

    return N2D 

end

# check
#N2D_calc(132.9,19.0, 11.0)

############################################################################
##############################################################################

function N2D_HH_LH(Ef::Float64, T::Float64)#N2D_HH_LH    # separa il contributo di carica iniziale tra HH LH (in genere è tutto HH anche a T=300K)
    
    NumPopStates= 4.0
    NumEmpStates=24.0
    d3K=0.0000001
    Aperp=117.5

    EDOUBLET = 5.0                   # Doublet energy
    KBOLTZMANN = 0.086173423         # Boltzmann constant (meV/K)
    EnArraySize = 500                # Energy array size
    ETrMax = 1000.0                  # Maximum energy transfer (meV)
    DIM = 1e14 / pi / 137.0 / 4.0    
    # 10^14/PI*costante di struttura fine / indice di rifrazione (germanio) con dipolo da codice e energie in meV  e d3k in A^-3 =#


    KT = KBOLTZMANN.*T
    d3K=d3K.*1.0 # here I consider the degeneracy of 1 due to the grid3dfull full positive and negative values for a1 a2 and z 
    #  d3K=d3K*8.d0 ! here I consider the degeneracy of 8 due to the grid3d not full with positive values only for a1 a2 and z 


    # Define the file path
    file_path_input = joinpath(@__DIR__, "InputData", "QWIP_VSx80_APTprofile_highGe_CMP3D.dat")
    # Attempt to read the file without a header
    data_df = CSV.read(file_path_input, DataFrame; delim=' ', ignorerepeated=true, header=false)
    # Assign default names if necessary
    col_names = [:k_cnt, :kx_cur, :ky_cur, :fill_cnt, :other_cnt, :EFill, :EOther, :PXX2, :PYY2, :PZZ2]
    rename!(data_df, col_names[1:ncol(data_df)])
    # Access arrays directly
    # k_cnt_arr = data_df[:, :k_cnt]
    # kx_cur_arr = data_df[:, :kx_cur]
    # ky_cur_arr = data_df[:, :ky_cur]
    fill_cnt_arr = data_df[:, :fill_cnt]
    other_cnt_arr = data_df[:, :other_cnt]
    EFill_arr = data_df[:, :EFill] .* 1000.0
    EOther_arr = data_df[:, :EOther] .* 1000.0
    # PXX2_arr = data_df[:, :PXX2]
    # PYY2_arr = data_df[:, :PYY2]
    # PZZ2_arr = data_df[:, :PZZ2]


    absEtr_arr = abs.(EFill_arr .- EOther_arr)
    # Find indices of elements greater than the threshold
    indices = findall(x -> x > ETrMax*1000.0, absEtr_arr)
    if !isempty(indices)
        @error("Error abs(EFill - EOther) > ETrMax.")
    end

    ### N2D calc
    # Calculate the target value
    target_value = NumEmpStates + NumPopStates
    # Find indices where the condition is satisfied
    indices = findall(x -> x == target_value, other_cnt_arr)
    N2D = d3K / (8.0 * π * π * π) * 1.0e24 * 1.0e-8 * Aperp * sum(getFermiDiracHoles(EFill_arr[indices[i]], Ef, KT) for i in eachindex(indices))

    ### N2D HH0 calc
    # Calculate the target value
    target_value = NumEmpStates + NumPopStates
    # Find indices where the condition is satisfied
    indices_0 = findall(x -> x == target_value, other_cnt_arr)
    indices_HH0 = findall(x -> ( (x == 28 || x==27) ), fill_cnt_arr)
    indices = intersect(indices_0, indices_HH0)
    N2D_HH0 = d3K / (8.0 * π * π * π) * 1.0e24 * 1.0e-8 * Aperp * sum(getFermiDiracHoles(EFill_arr[indices[i]], Ef, KT) for i in eachindex(indices))

    ### N2D LH calc
    # Calculate the target value
    target_value = NumEmpStates + NumPopStates
    # Find indices where the condition is satisfied
    indices_0 = findall(x -> x == target_value, other_cnt_arr)
    indices_LH = findall(x -> ( (x == 26 || x==25) ), fill_cnt_arr)
    indices = intersect(indices_0, indices_LH)
    if isempty(indices)
        N2D_LH = 0.0
    else
        N2D_LH = d3K / (8.0 * π * π * π) * 1.0e24 * 1.0e-8 * Aperp * sum(getFermiDiracHoles(EFill_arr[indices[i]], Ef, KT) for i in eachindex(indices))
    end
    

    # Return the calculated values
    return N2D, N2D_HH0, N2D_LH

end





#########################################################################################################################
######################################################################################################

function TM_max(Ef::Float64,T::Float64, HWHM::Float64; normalization=1.0, ETrMin = 105.0, ETrMax = 300.0)  # trova il picco di assorbimento
    
    NumPopStates= 4.0
    NumEmpStates=24.0
    d3K=0.0000001
    Aperp=117.5

    EDOUBLET = 5.0                   # Doublet energy
    KBOLTZMANN = 0.086173423         # Boltzmann constant (meV/K)
    EnArraySize = Int(round((ETrMax-ETrMin)/2))  #500                # Energy array size
    DIM = 1e14 / pi / 137.0 / 4.0    
    # 10^14/PI*costante di struttura fine / indice di rifrazione (germanio) con dipolo da codice e energie in meV  e d3k in A^-3 =#


    KT = KBOLTZMANN.*T
    d3K=d3K.*1.0 # here I consider the degeneracy of 1 due to the grid3dfull full positive and negative values for a1 a2 and z 
    #  d3K=d3K*8.d0 ! here I consider the degeneracy of 8 due to the grid3d not full with positive values only for a1 a2 and z 


    # Define the file path
    file_path_input = joinpath(@__DIR__, "InputData", "QWIP_VSx80_APTprofile_highGe_CMP3D.dat")
    # Attempt to read the file without a header
    data_df = CSV.read(file_path_input, DataFrame; delim=' ', ignorerepeated=true, header=false)
    # Assign default names if necessary
    col_names = [:k_cnt, :kx_cur, :ky_cur, :fill_cnt, :other_cnt, :EFill, :EOther, :PXX2, :PYY2, :PZZ2]
    rename!(data_df, col_names[1:ncol(data_df)])
    # Access arrays directly
    # k_cnt_arr = data_df[:, :k_cnt]
    # kx_cur_arr = data_df[:, :kx_cur]
    # ky_cur_arr = data_df[:, :ky_cur]
    # fill_cnt_arr = data_df[:, :fill_cnt]
    # other_cnt_arr = data_df[:, :other_cnt]
    EFill_arr = data_df[:, :EFill] .* 1000.0
    EOther_arr = data_df[:, :EOther] .* 1000.0
    # PXX2_arr = data_df[:, :PXX2]
    # PYY2_arr = data_df[:, :PYY2]
    PZZ2_arr = data_df[:, :PZZ2]


    absEtr_arr = abs.(EFill_arr .- EOther_arr)
    # Find indices of elements greater than the threshold
    indices = findall(x -> x > ETrMax*1000.0, absEtr_arr)
    if !isempty(indices)
        @error("Error abs(EFill - EOther) > ETrMax.")
    end


    # AbsorptionTE, AbsorptionTM calc
    indices = findall(x -> x > EDOUBLET, absEtr_arr)
    E_arr = zeros(EnArraySize)
    # calc
    absEtr_arr_restr = [absEtr_arr[indices[i]] for i in eachindex(indices)]
    PZZ2_arr_restr = [PZZ2_arr[indices[i]] for i in eachindex(indices)]
    EFill_arr_restr = [EFill_arr[indices[i]] for i in eachindex(indices)]
    EOther_arr_restr = [EOther_arr[indices[i]] for i in eachindex(indices)]
    FermiEFill_restr = getFermiDiracHoles.(EFill_arr_restr, Ef, KT)
    FermiEOther_restr = getFermiDiracHoles.(EOther_arr_restr, Ef, KT)
    E_arr = [ETrMin + EBinCnt * ETrMax / EnArraySize for EBinCnt in 1:EnArraySize]

    # Precompute invariant terms
    const_factor = DIM * d3K * HWHM / π
    HWHM2 = HWHM^2


    AbsorptionTM = [
        const_factor * sum(
            @. 1 / (HWHM2 + (E_arr[EBinCnt] - absEtr_arr_restr).^2) *
            PZZ2_arr_restr * FermiEFill_restr * (1.0 - FermiEOther_restr) / 
            copysign(E_arr[EBinCnt], EFill_arr_restr - EOther_arr_restr)
        )
        for EBinCnt in eachindex(E_arr)
    ]

    

    return maximum(AbsorptionTM)*normalization
end

# check
#TM_max(132.9,19.0, 11.0)

#############################################################################################
#####################################################################################

function update_or_create_csv(file_path::String, T::Float64, HWHM::Float64, N2D::Float64, Ef::Float64)  # crea o aggiorna il file 
    # Ensure the file is correctly handled
    if isfile(file_path)
        # Read the CSV with proper settings
        df = CSV.read(file_path, DataFrame; header=true, delim=',')
    else
        # Create a new DataFrame if the file does not exist
        df = DataFrame(T=Float64[], HW=Float64[], n=Float64[], mu=Float64[])
    end

    # Check if a row with the given T exists
    row_index = findfirst(row -> row.T == T, eachrow(df))
    
    if row_index !== nothing
        # Update the existing row with the new values
        df[row_index, :] .= (T, HWHM, N2D, Ef)
    else
        # Add a new row if T is not present
        push!(df, (T, HWHM, N2D, Ef))
    end

    # Write the updated DataFrame back to the file with correct formatting
    CSV.write(file_path, df; writeheader=true, delim=',')
end

function Absorption_calc(Ef::Float64,T::Float64, HWHM::Float64; normalization=1.0)
    
    NumPopStates= 4.0
    NumEmpStates=24.0
    d3K=0.0000001
    Aperp=117.5

    EDOUBLET = 5.0                   # Doublet energy
    KBOLTZMANN = 0.086173423         # Boltzmann constant (meV/K)
    EnArraySize = 500                # Energy array size
    ETrMax = 1000.0                  # Maximum energy transfer (meV)
    DIM = 1e14 / pi / 137.0 / 4.0    
    # 10^14/PI*costante di struttura fine / indice di rifrazione (germanio) con dipolo da codice e energie in meV  e d3k in A^-3 =#


    KT = KBOLTZMANN.*T
    d3K=d3K.*1.0 # here I consider the degeneracy of 1 due to the grid3dfull full positive and negative values for a1 a2 and z 
    #  d3K=d3K*8.d0 ! here I consider the degeneracy of 8 due to the grid3d not full with positive values only for a1 a2 and z 


    # Define the file path
    file_path_input = joinpath(@__DIR__, "InputData", "QWIP_VSx80_APTprofile_highGe_CMP3D.dat")
    # Attempt to read the file without a header
    data_df = CSV.read(file_path_input, DataFrame; delim=' ', ignorerepeated=true, header=false)
    # Assign default names if necessary
    col_names = [:k_cnt, :kx_cur, :ky_cur, :fill_cnt, :other_cnt, :EFill, :EOther, :PXX2, :PYY2, :PZZ2]
    rename!(data_df, col_names[1:ncol(data_df)])
    # Access arrays directly
    # k_cnt_arr = data_df[:, :k_cnt]
    # kx_cur_arr = data_df[:, :kx_cur]
    # ky_cur_arr = data_df[:, :ky_cur]
    # fill_cnt_arr = data_df[:, :fill_cnt]
    # other_cnt_arr = data_df[:, :other_cnt]
    EFill_arr = data_df[:, :EFill] .* 1000.0
    EOther_arr = data_df[:, :EOther] .* 1000.0
    PXX2_arr = data_df[:, :PXX2]
    # PYY2_arr = data_df[:, :PYY2]
    PZZ2_arr = data_df[:, :PZZ2]


    absEtr_arr = abs.(EFill_arr .- EOther_arr)
    # Find indices of elements greater than the threshold
    indices = findall(x -> x > ETrMax*1000.0, absEtr_arr)
    if !isempty(indices)
        @error("Error abs(EFill - EOther) > ETrMax.")
    end


    # AbsorptionTE, AbsorptionTM calc
    indices = findall(x -> x > EDOUBLET, absEtr_arr)
    E_arr = zeros(EnArraySize)
    # calc
    absEtr_arr_restr = [absEtr_arr[indices[i]] for i in eachindex(indices)]
    PXX2_arr_restr = [PXX2_arr[indices[i]] for i in eachindex(indices)]
    PZZ2_arr_restr = [PZZ2_arr[indices[i]] for i in eachindex(indices)]
    EFill_arr_restr = [EFill_arr[indices[i]] for i in eachindex(indices)]
    EOther_arr_restr = [EOther_arr[indices[i]] for i in eachindex(indices)]
    FermiEFill_restr = getFermiDiracHoles.(EFill_arr_restr, Ef, KT)
    FermiEOther_restr = getFermiDiracHoles.(EOther_arr_restr, Ef, KT)

    # energies
    E_arr = [0.0 + EBinCnt * ETrMax / EnArraySize for EBinCnt in 1:EnArraySize]

    # Precompute invariant terms
    const_factor = normalization * DIM * d3K * HWHM / π
    HWHM2 = HWHM^2


    AbsorptionTM = [
        const_factor * sum(
            @. 1 / (HWHM2 + (E_arr[EBinCnt] - absEtr_arr_restr).^2) *
            PZZ2_arr_restr * FermiEFill_restr * (1.0 - FermiEOther_restr) / 
            copysign(E_arr[EBinCnt], EFill_arr_restr - EOther_arr_restr)
        )
        for EBinCnt in eachindex(E_arr)
    ]

    AbsorptionTE = [
        const_factor * sum(
            @. 1 / (HWHM2 + (E_arr[EBinCnt] - absEtr_arr_restr).^2) *
            PXX2_arr_restr * FermiEFill_restr * (1.0 - FermiEOther_restr) / 
            copysign(E_arr[EBinCnt], EFill_arr_restr - EOther_arr_restr)
        )
        for EBinCnt in eachindex(E_arr)
    ]
    
    
    # Define the output path with the dynamic filename
    #output_file_dat = joinpath(@__DIR__, "OutputData", string("QWIP_VSx80_APTprofile_highGe_CMP3D.dat", ".Ef", Ef, "meV.", T, "K.alpha3D.dat"))
    output_file_dat = joinpath(@__DIR__, "OutputData", string("QWIP_VSx80_APTprofile_highGe_CMP3D.dat.", "T",T, "K" ,".alpha3D.dat"))

    # Open the file in write mode and write the data for each EBinCnt
    open(output_file_dat, "w") do file
        for EBinCnt in 1:EnArraySize
            # Calculate E based on EBinCnt
            E = 0.0 + EBinCnt * ETrMax / EnArraySize
            # Calculate values based on the given format
            absorption_tm = AbsorptionTM[EBinCnt]
            absorption_te = AbsorptionTE[EBinCnt]
            difference = absorption_tm - absorption_te
            # Write the line in the specified format
            write(file, "$E, $absorption_tm, $absorption_te, $difference\n")
        end
    end

    N2D = N2D_calc(Ef,T)
    file_path = joinpath(@__DIR__, "OutputData", "11229_data.csv")
    update_or_create_csv(file_path, T, HWHM, N2D, Ef)
    
end 

# check
#Absorption_calc(132.9,19.0, 11.0)

#################################################################################
##############################################################################


function plot_check(files_theo::Vector{String}, files_exp::Vector{String}, save_file::String = ""; shift::Bool = true)
    # Create a single-panel layout
    p = plot(layout = (1, 1), size = (1200, 600))  # Single plot layout

    # Initialize variables to store experimental maxima for alignment
    max_exp_val = -Inf
    max_exp_x = nothing

    # Plot experimental data
    for (i, file) in enumerate(files_exp)
        if isfile(file)
            temp_match = match(r"(\d+)K", file)
            temp_label = temp_match !== nothing ? temp_match.captures[1] * " K" : "Unknown T"

            energy = Float64[]
            absorption = Float64[]

            open(file, "r") do f
                readline(f)  # Skip the header line
                for line in eachline(f)
                    fields = split(line)
                    if length(fields) >= 2
                        try
                            energy_val = parse(Float64, fields[1])
                            absorption_val = parse(Float64, fields[2])
                            if 100 <= energy_val <= 400#300
                                push!(energy, energy_val)
                                push!(absorption, absorption_val)
                            end
                        catch
                            println("Skipping invalid line in file $file: $line")
                        end
                    end
                end
            end

            if !isempty(energy) && !isempty(absorption)
                # Update experimental maximum for alignment
                max_idx = argmax(absorption)
                if absorption[max_idx] > max_exp_val
                    max_exp_val = absorption[max_idx]
                    max_exp_x = energy[max_idx]
                end

                plot!(p, energy, absorption, label = temp_label, color = :black, linewidth = 2, linestyle = :dash)
            else
                println("No valid data found in file $file")
            end
        else
            println("File not found: $file")
        end
    end

    # Plot theoretical data
    for (i, file) in enumerate(files_theo)
        if isfile(file)
            x = Float64[]
            y = Float64[]

            open(file, "r") do f
                for line in eachline(f)
                    fields = split(line, ',')
                    if length(fields) >= 2
                        try
                            x_val = parse(Float64, fields[1])  # First column is x
                            y_val = parse(Float64, fields[2])  # Second column is y
                            if 100 <= x_val <= 400#300
                                push!(x, x_val)
                                push!(y, y_val)
                            end
                        catch
                            println("Skipping invalid line in file $file: $line")
                        end
                    end
                end
            end

            if !isempty(x) && !isempty(y)
                temp_label = "Theo $(i)"  # Label for theoretical files

                # Calculate shift along x if needed
                if shift && max_exp_x !== nothing
                    max_theo_idx = argmax(y)
                    max_theo_x = x[max_theo_idx]
                    x_shift = max_exp_x - max_theo_x
                    x .= x .+ x_shift  # Apply the shift
                end

                plot!(p, x, y, label = temp_label, color = :black, linewidth = 2)
            else
                println("No valid data found in file $file")
            end
        else
            println("File not found: $file")
        end
    end

    # Set labels for axes
    plot!(p, xlabel = "Energy (meV)", ylabel = "Absorption", title = "Theoretical and Experimental Data")

    # Save the plot if a save_file is provided
    if !isempty(save_file)
        output_path = joinpath(@__DIR__, save_file)
        savefig(p, output_path)
        println("Plot saved to: $output_path")
    end

    return p
end



#################################################################################
######################################################################################

default_plot_path = joinpath(@__DIR__, "OutputData", "11229_data.png")

function plot_theoretical_experimental(
    files_theo::Vector{String}, 
    files_exp::Vector{String}; 
    save_plot::Bool = false, 
    output_path::String = default_plot_path
)
    # Colors for each temperature (expand as needed)
    colors = [:blue, :cyan, :green, :gold, :red]
    
    # Create a 2-panel layout: top for theoretical, bottom for experimental
    p = plot(layout = (2, 1), size = (1200, 900), left_margin = 10mm)  # Increased left margin
    
    # Plot theoretical data (Top Panel)
    for (i, file) in enumerate(files_theo)
        if isfile(file)
            x = Float64[]
            y = Float64[]
            
            open(file, "r") do f
                for line in eachline(f)
                    # Split by comma
                    fields = split(line, ',')
                    if length(fields) >= 2
                        try
                            x_val = parse(Float64, fields[1])  # First column is x
                            y_val = parse(Float64, fields[2])  # Second column is y
                            if 100 <= x_val <= 300
                                push!(x, x_val)
                                push!(y, y_val)
                            end
                        catch
                            println("Skipping invalid line in file $file: $line")
                        end
                    end
                end
            end
            
            if !isempty(x) && !isempty(y)
                temp_label = ""#"File $(i)"  # Example label
                plot!(p[1], x, y, label = temp_label, color = colors[mod1(i, length(colors))], linewidth = 2)
            else
                println("No valid data found in file $file")
            end
        else
            println("File not found: $file")
        end
    end
    plot!(p[1], xlabel = "Energy (meV)", ylabel = "Absorption", title = "Theoretical")
    
    # Plot experimental data (Bottom Panel)
    for (i, file) in enumerate(files_exp)
        if isfile(file)
            temp_match = match(r"(\d+)K", file)
            temp_label = temp_match !== nothing ? temp_match.captures[1] * " K" : "Unknown T"
            
            energy = Float64[]
            absorption = Float64[]
            
            open(file, "r") do f
                readline(f)  # Skip the header line
                for line in eachline(f)
                    fields = split(line)
                    if length(fields) >= 2
                        try
                            energy_val = parse(Float64, fields[1])
                            absorption_val = parse(Float64, fields[2])
                            if 100 <= energy_val <= 300
                                push!(energy, energy_val)
                                push!(absorption, absorption_val)
                            end
                        catch
                            println("Skipping invalid line in file $file: $line")
                        end
                    end
                end
            end
            
            if !isempty(energy) && !isempty(absorption)
                plot!(p[2], energy, absorption, label = temp_label, color = colors[mod1(i, length(colors))], linewidth = 2)
            else
                println("No valid data found in file $file")
            end
        else
            println("File not found: $file")
        end
    end
    plot!(p[2], xlabel = "Energy (meV)", ylabel = "Absorption", title = "Experimental")
    
    # Save plot if requested
    if save_plot
        if isempty(output_path)
            println("Output path is empty. Please provide a valid path.")
        else
            savefig(p, output_path)
            println("Plot saved to: $output_path")
        end
    end
    
    return p
end



####################################################################################################
#############################################################################################

function TM_area(Ef::Float64,T::Float64, HWHM::Float64; normalization=1.0, ETrMin = 105.0, ETrMax = 285.0)
    
    NumPopStates= 4.0
    NumEmpStates=24.0
    d3K=0.0000001
    Aperp=117.5

    EDOUBLET = 5.0                   # Doublet energy
    KBOLTZMANN = 0.086173423         # Boltzmann constant (meV/K)
    EnArraySize = Int(round((ETrMax-ETrMin)/2))               # Energy array size
    DIM = 1e14 / pi / 137.0 / 4.0    
    # 10^14/PI*costante di struttura fine / indice di rifrazione (germanio) con dipolo da codice e energie in meV  e d3k in A^-3 =#


    KT = KBOLTZMANN.*T
    d3K=d3K.*1.0 # here I consider the degeneracy of 1 due to the grid3dfull full positive and negative values for a1 a2 and z 
    #  d3K=d3K*8.d0 ! here I consider the degeneracy of 8 due to the grid3d not full with positive values only for a1 a2 and z 


    # Define the file path
    file_path_input = joinpath(@__DIR__, "InputData", "QWIP_VSx80_APTprofile_highGe_CMP3D.dat")
    # Attempt to read the file without a header
    data_df = CSV.read(file_path_input, DataFrame; delim=' ', ignorerepeated=true, header=false)
    # Assign default names if necessary
    col_names = [:k_cnt, :kx_cur, :ky_cur, :fill_cnt, :other_cnt, :EFill, :EOther, :PXX2, :PYY2, :PZZ2]
    rename!(data_df, col_names[1:ncol(data_df)])
    # Access arrays directly
    # k_cnt_arr = data_df[:, :k_cnt]
    # kx_cur_arr = data_df[:, :kx_cur]
    # ky_cur_arr = data_df[:, :ky_cur]
    # fill_cnt_arr = data_df[:, :fill_cnt]
    # other_cnt_arr = data_df[:, :other_cnt]
    EFill_arr = data_df[:, :EFill] .* 1000.0
    EOther_arr = data_df[:, :EOther] .* 1000.0
    # PXX2_arr = data_df[:, :PXX2]
    # PYY2_arr = data_df[:, :PYY2]
    PZZ2_arr = data_df[:, :PZZ2]


    absEtr_arr = abs.(EFill_arr .- EOther_arr)
    # Find indices of elements greater than the threshold
    indices = findall(x -> x > ETrMax*1000.0, absEtr_arr)
    if !isempty(indices)
        @error("Error abs(EFill - EOther) > ETrMax.")
    end


    # AbsorptionTE, AbsorptionTM calc
    indices = findall(x -> x > EDOUBLET, absEtr_arr)
    E_arr = zeros(EnArraySize)
    # calc
    absEtr_arr_restr = [absEtr_arr[indices[i]] for i in eachindex(indices)]
    PZZ2_arr_restr = [PZZ2_arr[indices[i]] for i in eachindex(indices)]
    EFill_arr_restr = [EFill_arr[indices[i]] for i in eachindex(indices)]
    EOther_arr_restr = [EOther_arr[indices[i]] for i in eachindex(indices)]
    FermiEFill_restr = getFermiDiracHoles.(EFill_arr_restr, Ef, KT)
    FermiEOther_restr = getFermiDiracHoles.(EOther_arr_restr, Ef, KT)

    # energies
    E_arr = [ETrMin + EBinCnt * ETrMax / EnArraySize for EBinCnt in 1:EnArraySize]

    # Precompute invariant terms
    const_factor = normalization * DIM * d3K * HWHM / π
    HWHM2 = HWHM^2


    AbsorptionTM = [
        const_factor * sum(
            @. 1 / (HWHM2 + (E_arr[EBinCnt] - absEtr_arr_restr).^2) *
            PZZ2_arr_restr * FermiEFill_restr * (1.0 - FermiEOther_restr) / 
            copysign(E_arr[EBinCnt], EFill_arr_restr - EOther_arr_restr)
        )
        for EBinCnt in eachindex(E_arr)
    ]
    
    # trapezoidal method
    area = sum(0.5 * (E_arr[i+1] - E_arr[i]) * (AbsorptionTM[i] + AbsorptionTM[i+1]) for i in 1:length(E_arr)-1)


    return area
end

############################################################################################
####################################################################################

function area_exp(energy_exp, absorption_coefficient_exp; ETrMax=285.0)
    # Filter the indices where energy_exp is less than or equal to ETrMax
    indices = findall(x -> x <= ETrMax, energy_exp)

    # Ensure the last point doesn't exceed ETrMax
    last_index = indices[end]
    if energy_exp[last_index] < ETrMax
        # Add a point at ETrMax using linear interpolation
        ΔE = energy_exp[last_index + 1] - energy_exp[last_index]
        ΔAbsorption = absorption_coefficient_exp[last_index + 1] - absorption_coefficient_exp[last_index]
        interpolated_absorption = absorption_coefficient_exp[last_index] + (ΔAbsorption * (ETrMax - energy_exp[last_index]) / ΔE)
        energy_exp = vcat(energy_exp[1:last_index], ETrMax)
        absorption_coefficient_exp = vcat(absorption_coefficient_exp[1:last_index], interpolated_absorption)
        last_index += 1
    end

    # Compute the area under the curve using the trapezoidal rule
    area = sum(0.5 * (energy_exp[i+1] - energy_exp[i]) * (absorption_coefficient_exp[i] + absorption_coefficient_exp[i+1]) for i in 1:last_index-1)

    return area
end

#############################################################
#############################################################################


# Define a function that finds Ef given T, target_N2D, and tolerance percentage
function find_optimal_Ef(T::Float64, HWHM::Float64, target_N2D::Float64, Ef::Float64, tolerance_percent::Float64)
    # Calculate Ef range based on tolerance percentage
    Ef_min = Ef * (1.0 - tolerance_percent / 100.0)
    Ef_max = Ef * (1.0 + tolerance_percent / 100.0)

    # Define the objective function to minimize the difference between dens(Ef, T) and target N2D
    function objective(Ef)
        N2D = N2D_calc(Ef, T)
        return abs(N2D - target_N2D)  # Calculate the absolute difference
    end

    # Use an optimization method to find the Ef that minimizes the objective function within the range
    result = optimize(objective, Ef_min, Ef_max)

    # Extract the optimal Ef and calculate the corresponding N2D
    optimal_Ef = Optim.minimizer(result)
    optimal_N2D = N2D_calc(optimal_Ef, T)

    return optimal_Ef, optimal_N2D
end

###############################################################
##################################################################

function update_N2D_res(T::Float64, N2D_check::Float64, N2D_HH0::Float64, N2D_LH::Float64; file_path=joinpath(@__DIR__, "11229_N2D_res.csv"))
    # Check if the file exists
    if isfile(file_path)
        # Read the existing CSV file, explicitly specifying the delimiter
        df = CSV.read(file_path, DataFrame; delim=',', ignorerepeated=true)
    else
        # Create a new DataFrame if the file doesn't exist
        df = DataFrame(T = Float64[], N2D_check = Float64[], N2D_HH0 = Float64[], N2D_LH = Float64[])
    end

    # Check if a row with the given T_10K already exists
    existing_row = findfirst(row -> isapprox(row[:T], T, atol=1e-6), eachrow(df))

    if isnothing(existing_row)
        # Add a new row if it doesn't exist
        push!(df, (T, N2D_check, N2D_HH0, N2D_LH))
    else
        # Update the existing row
        df[existing_row, :] = (T, N2D_check, N2D_HH0, N2D_LH)
    end

    # Write the updated DataFrame back to the file
    CSV.write(file_path, df; writeheader=true)
end


###############################################################
###########################################

### 11229
### We use FTIR values at low temperatures: N2D = 5.16e11 cm-3 and HWHM = 11.3

### Optimisation to find Ef


T_10K = 10.0
HWHM_10K = 10.0
N2D_10K = 5.16e11

opt_Ef_10K, opt_N2D_10 = find_optimal_Ef(T_10K, HWHM_10K, N2D_10K, 146.0, 30.0)

### res
N2D_check, N2D_HH0, N2D_LH = N2D_HH_LH(opt_Ef_10K ,T_10K) 
update_N2D_res(T_10K, N2D_check, N2D_HH0, N2D_LH)

## Maximum experimantal data
# Define the file path
file_path = joinpath(@__DIR__, "11229_exp", "11229 absorption 10K.txt")
# Read the file into a DataFrame
data = CSV.read(file_path, DataFrame; header=true, delim='\t')
# Extract columns as arrays
energy_exp_10K = reverse(data[!, 1])  # Column "energy (meV)"
absorption_coefficient_exp_10K = reverse(data[!, 2])  # Column "absorption coefficient"
maximum_TM_exp = maximum(absorption_coefficient_exp_10K)

## Maximum theoretical data
maximum_TM_theo = TM_max(opt_Ef_10K,10.0, HWHM_10K)

norm_factor = maximum_TM_exp/maximum_TM_theo

### Theoretical data computation for 10.0K
norm_factor
Absorption_calc(opt_Ef_10K, 10.0, HWHM_10K; normalization=norm_factor)


### Plot
files_theo = [
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T10.0K.alpha3D.dat"),
]
files_exp = [
    joinpath(@__DIR__, "11229_exp/11229 absorption 10K.txt"),
]

plot_check(files_theo, files_exp, "comp.10K.png")
plot_check(files_theo, files_exp, "comp.10K.png", shift=false)



#########################################################################################################################
######################################################################################################

# Experimental areas  ETrMin = 105.44, ETrMax = 285.0

area_exp_10K = area_exp(energy_exp_10K, absorption_coefficient_exp_10K, ETrMax=285.0)
area_theo_10K = TM_area(opt_Ef_10K,T_10K, HWHM_10K; normalization=norm_factor, ETrMin = 105.0, ETrMax = 285.0)
# target ratio
target_area_ratio_theo_exp = area_theo_10K/area_exp_10K



#######################################################################à
#################################################à


### 50K
T_50K = 50.0

# Define the file path
file_path = joinpath(@__DIR__, "11229_exp", "11229 absorption 50K.txt")
# Read the file into a DataFrame
data = CSV.read(file_path, DataFrame; header=true, delim='\t')
# Extract columns as arrays
energy_exp_50K = reverse(data[!, 1])  # Column "energy (meV)"
absorption_coefficient_exp_50K = reverse(data[!, 2])  # Column "absorption coefficient"

area_exp_50K = area_exp(energy_exp_50K, absorption_coefficient_exp_50K, ETrMax=285.0)
area_target_50K = area_exp_50K * target_area_ratio_theo_exp

####

# Initial guess for Ef and HWHM; Ef estimation from FTIR
N2D_FTIR_50K = 5.10e11
Ef_initial_guess, target_N2D = find_optimal_Ef(T_50K, HWHM_10K, N2D_FTIR_50K, 100.0, 70.0)
opt_Ef_10K

##
maximum_target = maximum(absorption_coefficient_exp_50K)

# This function computes the absolute difference between the calculated area and the target area,
# and ensures TM_max(Ef, T_50K, HWHM) is close to maximum_target
function objective(params)
    Ef, HWHM = params
    TM_area_value = TM_area(Ef, T_50K, HWHM; normalization = norm_factor)  # Compute the area using the given parameters
    TM_max_value = TM_max(Ef, T_50K, HWHM; normalization=norm_factor)  # Compute the maximum value using the given parameters
    
    # Compute the differences
    area_difference = abs(TM_area_value - area_target_50K)  # Difference to the target area
    max_difference = abs(TM_max_value - maximum_target)  # Difference to the target maximum
    
    # Combine the two differences (you can adjust the weight if necessary)
    return area_difference + max_difference
end

# Initial guess for Ef and HWHM
initial_guess = [Ef_initial_guess, 9.0]  # Example: Ef = 150.0, HWHM = 10.0
#initial_guess = [149.58, 11.3]

# Perform optimization using the Nelder-Mead method
result = optimize(objective, initial_guess, NelderMead())

# Extract the optimal values for Ef and HWHM
opt_Ef_50K, opt_HWHM_50K = Optim.minimizer(result)
opt_Ef_10K
#opt_HWHM_10K
# generaation file
Absorption_calc(opt_Ef_50K, T_50K, opt_HWHM_50K, normalization=norm_factor)

opt_HWHM_50K
# check
TM_area(opt_Ef_50K, T_50K, opt_HWHM_50K; normalization=norm_factor) - area_target_50K
### Plot
files_theo = [
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T50.0K.alpha3D.dat"),
]
files_exp = [
    joinpath(@__DIR__, "11229_exp/11229 absorption 50K.txt"),
]

plot_check(files_theo, files_exp, "comp.50K.png")
plot_check(files_theo, files_exp, "comp.50K.png",shift=false)

### res
N2D_check, N2D_HH0, N2D_LH = N2D_HH_LH(opt_Ef_50K ,T_50K) 
update_N2D_res(T_50K, N2D_check, N2D_HH0, N2D_LH)



### 100K
T_100K = 100.0

# Define the file path
file_path = joinpath(@__DIR__, "11229_exp", "11229 absorption 100K.txt")
# Read the file into a DataFrame
data = CSV.read(file_path, DataFrame; header=true, delim='\t')
# Extract columns as arrays
energy_exp_100K = reverse(data[!, 1])  # Column "energy (meV)"
absorption_coefficient_exp_100K = reverse(data[!, 2])  # Column "absorption coefficient"

area_exp_100K = area_exp(energy_exp_100K, absorption_coefficient_exp_100K, ETrMax=285.0)
area_target_100K = area_exp_100K * target_area_ratio_theo_exp

####

# Initial guess for Ef and HWHM; Ef estimation from FTIR
N2D_FTIR_100K = 5.19e11
Ef_initial_guess, target_N2D = find_optimal_Ef(T_100K, opt_HWHM_50K, N2D_FTIR_100K, 100.0, 200.0)
opt_Ef_50K

##
maximum_target = maximum(absorption_coefficient_exp_100K)


# This function computes the absolute difference between the calculated area and the target area,
# and ensures TM_max(Ef, T_50K, HWHM) is close to maximum_target
function objective(params)
    Ef, HWHM = params
    TM_area_value = TM_area(Ef, T_100K, HWHM; normalization = norm_factor)  # Compute the area using the given parameters
    TM_max_value = TM_max(Ef, T_100K, HWHM; normalization=norm_factor)  # Compute the maximum value using the given parameters
    
    # Compute the differences
    area_difference = abs(TM_area_value - area_target_100K)  # Difference to the target area
    max_difference = abs(TM_max_value - maximum_target)  # Difference to the target maximum
    
    # Combine the two differences (you can adjust the weight if necessary)
    return area_difference + max_difference
end

# Initial guess for Ef and HWHM
initial_guess = [Ef_initial_guess+10, 6.0]  # Example: Ef = 150.0, HWHM = 10.0
#initial_guess = [149.58, 11.3]
opt_HWHM_50K
# Perform optimization using the Nelder-Mead method
result = optimize(objective, initial_guess, NelderMead())

# Extract the optimal values for Ef and HWHM
opt_Ef_100K, opt_HWHM_100K = Optim.minimizer(result)
opt_Ef_100K
opt_HWHM_100K
# generaation file
Absorption_calc(opt_Ef_100K, T_100K, opt_HWHM_100K, normalization=norm_factor)


# check
TM_area(opt_Ef_100K, T_100K, opt_HWHM_100K; normalization=norm_factor) - area_target_100K
### Plot
files_theo = [
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T100.0K.alpha3D.dat"),
]
files_exp = [
    joinpath(@__DIR__, "11229_exp/11229 absorption 100K.txt"),
]

plot_check(files_theo, files_exp, "comp.100K.png")
plot_check(files_theo, files_exp, "comp.100K.png", shift=false)


### res
N2D_check, N2D_HH0, N2D_LH = N2D_HH_LH(opt_Ef_100K ,T_100K) 
update_N2D_res(T_100K, N2D_check, N2D_HH0, N2D_LH)


### 200K
T_200K = 200.0

# Define the file path
file_path = joinpath(@__DIR__, "11229_exp", "11229 absorption 200K.txt")
# Read the file into a DataFrame
data = CSV.read(file_path, DataFrame; header=true, delim='\t')
# Extract columns as arrays
energy_exp_200K = reverse(data[!, 1])  # Column "energy (meV)"
absorption_coefficient_exp_200K = reverse(data[!, 2])  # Column "absorption coefficient"

area_exp_200K = area_exp(energy_exp_200K, absorption_coefficient_exp_200K, ETrMax=285.0)
area_target_200K = area_exp_200K * target_area_ratio_theo_exp

####

# Initial guess for Ef and HWHM; Ef estimation from FTIR
N2D_FTIR_200K = 4.06e11
Ef_initial_guess, target_N2D = find_optimal_Ef(T_200K, 5.0, N2D_FTIR_200K, 140.0, 200.0)
opt_Ef_100K

##
maximum_target = maximum(absorption_coefficient_exp_200K)


# This function computes the absolute difference between the calculated area and the target area,
# and ensures TM_max(Ef, T_50K, HWHM) is close to maximum_target
function objective(params)
    Ef, HWHM = params
    TM_area_value = TM_area(Ef, T_200K, HWHM; normalization = norm_factor)  # Compute the area using the given parameters
    TM_max_value = TM_max(Ef, T_200K, HWHM; normalization=norm_factor)  # Compute the maximum value using the given parameters
    
    # Compute the differences
    area_difference = abs(TM_area_value - area_target_200K)  # Difference to the target area
    max_difference = abs(TM_max_value - maximum_target)  # Difference to the target maximum
    
    # Combine the two differences (you can adjust the weight if necessary)
    return area_difference + max_difference
end

# Initial guess for Ef and HWHM
#initial_guess = [Ef_initial_guess, opt_HWHM_100K]  # Example: Ef = 150.0, HWHM = 10.0
initial_guess = [Ef_initial_guess+20, 5.0]

# Perform optimization using the Nelder-Mead method
result = optimize(objective, initial_guess, NelderMead())

# Extract the optimal values for Ef and HWHM
opt_Ef_200K, opt_HWHM_200K = Optim.minimizer(result)
opt_Ef_100K
opt_HWHM_100K
# generaation file
Absorption_calc(opt_Ef_200K, T_200K, opt_HWHM_200K, normalization=norm_factor)


# check
TM_area(opt_Ef_200K, T_200K, opt_HWHM_200K; normalization=norm_factor) - area_target_200K
### Plot
files_theo = [
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T200.0K.alpha3D.dat"),
]
files_exp = [
    joinpath(@__DIR__, "11229_exp/11229 absorption 200K.txt"),
]

plot_check(files_theo, files_exp, "comp.200K.png",shift=false)
### res
N2D_check, N2D_HH0, N2D_LH = N2D_HH_LH(opt_Ef_200K ,T_200K) 
update_N2D_res(T_200K, N2D_check, N2D_HH0, N2D_LH)



### 300K
T_300K = 300.0

# Define the file path
file_path = joinpath(@__DIR__, "11229_exp", "11229 absorption 300K.txt")
# Read the file into a DataFrame
data = CSV.read(file_path, DataFrame; header=true, delim='\t')
# Extract columns as arrays
energy_exp_300K = reverse(data[!, 1])  # Column "energy (meV)"
absorption_coefficient_exp_300K = reverse(data[!, 2])  # Column "absorption coefficient"

area_exp_300K = area_exp(energy_exp_300K, absorption_coefficient_exp_300K, ETrMax=285.0)
area_target_300K = area_exp_300K * target_area_ratio_theo_exp

####

# Initial guess for Ef and HWHM; Ef estimation from FTIR
N2D_FTIR_300K = 4.08e11
Ef_initial_guess, target_N2D = find_optimal_Ef(T_300K, 5.0, N2D_FTIR_300K, 140.0, 200.0)
opt_Ef_200K

##
maximum_target = maximum(absorption_coefficient_exp_300K)



# This function computes the absolute difference between the calculated area and the target area,
# and ensures TM_max(Ef, T_50K, HWHM) is close to maximum_target
function objective(params)
    Ef, HWHM = params
    TM_area_value = TM_area(Ef, T_300K, HWHM; normalization = norm_factor)  # Compute the area using the given parameters
    TM_max_value = TM_max(Ef, T_300K, HWHM; normalization=norm_factor)  # Compute the maximum value using the given parameters
    
    # Compute the differences
    area_difference = abs(TM_area_value - area_target_300K)  # Difference to the target area
    max_difference = abs(TM_max_value - maximum_target)  # Difference to the target maximum
    
    # Combine the two differences (you can adjust the weight if necessary)
    return area_difference + max_difference
end

# Initial guess for Ef and HWHM
#initial_guess = [Ef_initial_guess, opt_HWHM_100K]  # Example: Ef = 150.0, HWHM = 10.0
initial_guess = [Ef_initial_guess+20.0, 5.0]

# Perform optimization using the Nelder-Mead method
result = optimize(objective, initial_guess, NelderMead())

# Extract the optimal values for Ef and HWHM
opt_Ef_300K, opt_HWHM_300K = Optim.minimizer(result)
opt_Ef_200K
opt_HWHM_200K
# generaation file
Absorption_calc(opt_Ef_300K, T_300K, opt_HWHM_300K, normalization=norm_factor)


# check
TM_area(opt_Ef_300K, T_300K, opt_HWHM_300K; normalization=norm_factor) - area_target_300K
### Plot
files_theo = [
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T300.0K.alpha3D.dat"),
]
files_exp = [
    joinpath(@__DIR__, "11229_exp/11229 absorption 300K.txt"),
]

plot_check(files_theo, files_exp, "comp.300K.png",shift=false)

### res
N2D_check, N2D_HH0, N2D_LH = N2D_HH_LH(opt_Ef_300K ,T_300K) 
update_N2D_res(T_300K, N2D_check, N2D_HH0, N2D_LH)


##############################################################
############################################################
# Plot

files_theo = [
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T10.0K.alpha3D.dat"),
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T50.0K.alpha3D.dat"),
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T100.0K.alpha3D.dat"),
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T200.0K.alpha3D.dat"),
    joinpath(@__DIR__, "OutputData/QWIP_VSx80_APTprofile_highGe_CMP3D.dat.T300.0K.alpha3D.dat")
]

files_exp = [
    joinpath(@__DIR__, "11229_exp/11229 absorption 10K.txt"),
    joinpath(@__DIR__, "11229_exp/11229 absorption 50K.txt"),
    joinpath(@__DIR__, "11229_exp/11229 absorption 100K.txt"),
    joinpath(@__DIR__, "11229_exp/11229 absorption 200K.txt"),
    joinpath(@__DIR__, "11229_exp/11229 absorption 300K.txt")
]


plot_theoretical_experimental(files_theo, files_exp; save_plot = true, output_path= default_plot_path)


function plot_theoretical_scaled(
    files_theo::Vector{String},
    scale_factor::Float64;                 # Fattore di scala per le y
    save_plot::Bool = false,
    output_path::String = "theoretical_non_scaled_plot.png"
)
    # Colori ciclici per i plot
    colors = [:blue, :cyan, :green, :gold, :red]

    # Crea un layout singolo
    p = plot(layout = (1, 1), size = (1200, 600), left_margin = 10mm)

    # Loop sui file teorici
    for (i, file) in enumerate(files_theo)
        if isfile(file)
            x = Float64[]
            y = Float64[]

            open(file, "r") do f
                for line in eachline(f)
                    fields = split(line, ',')
                    if length(fields) >= 2
                        try
                            x_val = parse(Float64, fields[1])
                            y_val = parse(Float64, fields[2]) * scale_factor  # Applica il fattore di scala
                            if 100 <= x_val <= 300
                                push!(x, x_val)
                                push!(y, y_val)
                            end
                        catch
                            println("Skipping invalid line in file $file: $line")
                        end
                    end
                end
            end

            if !isempty(x) && !isempty(y)
                temp_label = "Theo $(i)"
                plot!(p[1], x, y, label = temp_label, color = colors[mod1(i, length(colors))], linewidth = 2)
            else
                println("No valid data found in file $file")
            end
        else
            println("File not found: $file")
        end
    end

    plot!(p[1], xlabel = "Energy (meV)", ylabel = "Non-Scaled Absorption", title = "Theoretical")

    # Salva il grafico se richiesto
    if save_plot
        savefig(p, output_path)
        println("Plot saved to: $output_path")
    end

    return p
end

#plot_theoretical_experimental(files_theo, files_exp; save_plot = true, output_path= default_plot_path)
default_plot_path_scaled = joinpath(@__DIR__, "OutputData", "11229_data_scaled.png")

plot_theoretical_scaled(files_theo,1.0/norm_factor,save_plot = true, output_path= default_plot_path_scaled)