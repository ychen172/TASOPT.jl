"""
This script compares parameters across multiple cases
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel, init_results_2Layers, fill_results!, plot_cases

#### Setup IO
# Input case names - Retrofit
model_dir  = "../ModelSaved"
caseKeys   = ["Opti_Jet_NoACT_",   "Opti_Eth_NoACT_",  "Opti_Eth_NoACT_NoWSpan_", "Opti_Eth_ACT_NoComp_UniEta_NoWSpan_"]
caseNames  = ["Jet Fuel Baseline", "Ethanol Baseline", "Ethanol Extended Wing",   "Ethanol Extended Wing Ideal ACT"]
idxMiss    = [1,1,1,1]
ranges     = collect(300:100:3000) 
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "Wing_span_effects" #sub_folder will be created
# Fields to read out
const fields = (:range_nmi,:PFEI_JJ,:lenFuseCyl_m,:FuelVolumeFractionACT,:lenACT_m,:volFuelTot_m3,:massEmpty_Ton,
                :massFuelTot_Ton,:massPayload_Ton,:ene_fli_J,:volFuelMax_m3,:FuelVolumeFraction,:massTO_Ton,:span_wing_m,
                :LD_cru,:AR_wing,:CL_cruise,:Alt_cruise_ft)

#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
results = [init_results_2Layers(length(ranges), fields) for _ in caseKeys]

#### Extract data for the sized missions
for (i, caseKey) in enumerate(caseKeys)
    for (j, ran) in enumerate(ranges)
        # Read in the case
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(round(Int,ran)).jld2")
        ac = quickload_aircraft(ac_dir)
        println("File, $(ac_dir), read successfully")
        # Extract data
        out = extract_acModel(ac,idxMiss[i])
        # Store the data
        fill_results!(results[i],out,j)
    end
end

#### Plotting
plot_cases("Design Range","PFEI (J/J)",results,:range_nmi,:PFEI_JJ,caseNames,joinpath(save_dir_sub, "PFEI.png"))
plot_cases("Design Range","Fuselage Cylinder Length (m)",results,:range_nmi,:lenFuseCyl_m,caseNames,joinpath(save_dir_sub, "LCylinder.png"))
plot_cases("Design Range","Fractional Fuel in ACT",results,:range_nmi,:FuelVolumeFractionACT,caseNames,joinpath(save_dir_sub, "FractionFuelACT.png"))
plot_cases("Design Range","Length of ACT (m)",results,:range_nmi,:lenACT_m,caseNames,joinpath(save_dir_sub, "ACTLength.png"))
plot_cases("Design Range","Volume of Fuel (m³)",results,:range_nmi,:volFuelTot_m3,caseNames,joinpath(save_dir_sub, "FuelVolume.png"))
plot_cases("Design Range","Mass of Fuel (Ton)",results,:range_nmi,:massFuelTot_Ton,caseNames,joinpath(save_dir_sub, "FuelMass.png"))
plot_cases("Design Range","Empty Mass (Ton)",results,:range_nmi,:massEmpty_Ton,caseNames,joinpath(save_dir_sub, "EmptyMass.png"))
plot_cases("Design Range","Payload Mass (Ton)",results,:range_nmi,:massPayload_Ton,caseNames,joinpath(save_dir_sub, "PayloadMass.png"))
plot_cases("Design Range","Fuel Energy (J)",results,:range_nmi,:ene_fli_J,caseNames,joinpath(save_dir_sub, "FuelEnergy.png"))
plot_cases("Design Range","Wing Inner Volume (m³)",results,:range_nmi,:volFuelMax_m3,caseNames,joinpath(save_dir_sub, "FuelVolumeWing.png"))
plot_cases("Design Range","Fractional Fuel Volume Used",results,:range_nmi,:FuelVolumeFraction,caseNames,joinpath(save_dir_sub, "FracFuelVolumeUsed.png"))
plot_cases("Design Range","Takeoff Mass (Ton)",results,:range_nmi,:massTO_Ton,caseNames,joinpath(save_dir_sub, "TakeeoffMass.png"))
plot_cases("Design Range","Wing Span (m)",results,:range_nmi,:span_wing_m,caseNames,joinpath(save_dir_sub, "WingSpan.png"))
plot_cases("Design Range","Wing Aspect Ratio",results,:range_nmi,:AR_wing,caseNames,joinpath(save_dir_sub, "WingAR.png"))
plot_cases("Design Range","Cruise LD Ratio",results,:range_nmi,:LD_cru,caseNames,joinpath(save_dir_sub, "CruiseLD.png"))
plot_cases("Design Range","Cruise CL",results,:range_nmi,:CL_cruise,caseNames,joinpath(save_dir_sub, "CruiseCL.png"))
plot_cases("Design Range","Cruise Altitude (ft)",results,:range_nmi,:Alt_cruise_ft,caseNames,joinpath(save_dir_sub, "CruiseAlt.png"))