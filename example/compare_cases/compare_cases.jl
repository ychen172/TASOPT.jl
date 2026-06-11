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
caseKeys   = ["Ethanol_ACTPlane_Opt","Ethanol_NoACT_Opt","JetFuel_ACTPlane_Opt","JetFuel_NoACT_Opt"]
caseNames  = ["Ethanol with ACT","Ethanol without ACT","Jet Fuel with ACT","Jet Fuel without ACT"]
idxMiss    = [1,1,1,1]
ranges     = collect(300:100:2900) 
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "Compare_at_design" #sub_folder will be created
# Fields to read out
const fields = (:range_nmi,:PFEI_JJ,:lenFuseCyl_m,:FuelVolumeFractionACT,:lenACT_m)

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