"""
This script compares parameters across multiple cases
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel, init_results_2Layers

#### Setup IO
# Input case names - Retrofit
model_dir       = "../ModelSaved"
# Input case names - Ethanol
key_sized_Eth       = "acOptimized_BatOptEth" #180Pass_Opt: "acOptimized_BatOptEth", 3000nmiJet2EthRetroCase_but_optimized: "OptimizedJetToEth3000_"
des_range_sized_Eth = float.(collect(300:100:2900)) #Has to match with existings
# Input case names - Jet fuel
key_sized_Jet       = "acOptimized_BatOptJet" #180Pass_Opt: "CenterFuelTank_BatOptEth", 3000nmiJet2EthRetroCase_but_optimizedWithCenTank: "OptCenTankJetToEth3000_"
des_range_sized_Jet = float.(collect(300:100:3000)) #Has to match with existings
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "Compare_at_design" #sub_folder will be created
# Fields to read out
const fields = (:range_nmi,   :PFEI_JJ, :massTO_Ton, :massFuelTot_Ton, :massPayload_Ton,
                :massEmpty_Ton, :rhoFuel_kgm3,  :volFuel_m3, :LD_cru, :LHV_cru_Jkg, :TSFC_cru_kgsN,
                :vel_cru_ms, :eta_total_cru, :ene_fli_J, :frac_rese, :PFEI_cru_JJ, :range_cru_m,
                :eta_therm_cru, :spe_power_cru_Jkg, :eta_propu_cru, :OPR_cru, :Tt_turbin_cru_K, :FuelVolumeFraction)

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
results_sized_Eth = init_results_2Layers(length(des_range_sized_Eth), fields)
results_sized_Jet = init_results_2Layers(length(des_range_sized_Jet), fields)

#### Extract data for the sized missions
for (i, des_ranges_cur) in enumerate(des_range_sized_Eth)
    #### Read in the case
    ac_dir_sized = joinpath(model_dir,key_sized_Eth,key_sized_Eth*"$(round(Int,des_ranges_cur)).jld2")
    ac_sized = quickload_aircraft(ac_dir_sized)
    println("File, $(ac_dir_sized), read successfully")
    out_dict = extract_acModel(ac_sized,1)
    
    #### Output the data
    for f in fields
        results_sized_Eth[f][i] = out_dict[String(f)]
    end
end

for (i, des_ranges_cur) in enumerate(des_range_sized_Jet)
    #### Read in the case
    ac_dir_sized = joinpath(model_dir,key_sized_Jet,key_sized_Jet*"$(round(Int,des_ranges_cur)).jld2")
    ac_sized = quickload_aircraft(ac_dir_sized)
    println("File, $(ac_dir_sized), read successfully")
    out_dict = extract_acModel(ac_sized,1)
    
    #### Output the data
    for f in fields
        results_sized_Jet[f][i] = out_dict[String(f)]
    end
end

#### Plotting
# Create style
linestyles = repeat([:solid, :dash, :dot, :dashdot, :dashdotdot],1000)
linecolors = repeat([:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :black, :cyan,
                     :magenta, :teal, :navy, :maroon, :olive, :gold, :coral, :turquoise, :lime, :indigo], 1000) 
markers = repeat([:rect, :circle, :diamond, :utriangle, :dtriangle],1000)
# Plot PFEI
p1_1 = plot(xlabel="Design Ranges (nmi)", ylabel="PFEI (J/J)", dpi=800)
plot!(p1_1,ylims=(0.6,1.4))
global il = 1
plot!(p1_1, results_sized_Jet[:range_nmi], results_sized_Jet[:PFEI_JJ], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Optimized Jet Fuel Aircraft")
global il += 1
plot!(p1_1, results_sized_Eth[:range_nmi], results_sized_Eth[:PFEI_JJ], marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Optimized Ethanol Aircraft")
savefig(p1_1, joinpath(save_dir_sub, "PFEI.png"))

# Plot Fuel Fraction
p1_2 = plot(xlabel="Design Ranges (nmi)", ylabel="Fuel Tank Used (%)", dpi=800, legend=:bottomright)
global il = 1
plot!(p1_2, results_sized_Jet[:range_nmi], results_sized_Jet[:FuelVolumeFraction] .* 100.0, marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Optimized Jet Fuel Aircraft")
global il += 1
plot!(p1_2, results_sized_Eth[:range_nmi], results_sized_Eth[:FuelVolumeFraction] .* 100.0, marker=markers[il], mc=linecolors[il], msc=linecolors[il], color=linecolors[il], lw=2, linestyle=linestyles[il], label="Optimized Ethanol Aircraft")
savefig(p1_2, joinpath(save_dir_sub, "FuelTankUsed.png"))