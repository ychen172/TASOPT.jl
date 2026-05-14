"""
This script compare the design point performance between ethanol and jet
Fixed design bounds and Fixed constraints across different cases assumed
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(@__DIR__, "postprocess.jl"))
using .PostProcess: ExtractDes
using .PostProcess.OptimizeRangeFuel: BoundsOpt, ConstraintsOpt

#### Setup IO
model_dir   = "ModelSaved"
model_name  = "acOptimized_BatOpt" #Frontal key name for the models(FuelRange)
save_dir    = "ModelProcessed" #Outer Directory for saving models
save_name   = "EthVsJet" #The total save name for jet and ethanol comparison
mkpath(model_dir)
mkpath(save_dir)
save_dirSub = joinpath(save_dir,save_name) #Sub-directory to save comparison data
mkpath(save_dirSub)

#### Setup parameters
Fuels = ["Eth","Jet"]
Ranges = collect(300:100:3000)
# flg_save = false #Whether to save individual model outputs

#### Extract data
desParamRows = [] #Each row is one case
for (i,curFuel) in enumerate(Fuels)
    for (j,curRange) in enumerate(Ranges)
        case_id = "$(curFuel)$(curRange)" #Current Case
        model_file = "$(model_name)$(case_id).jld2"
        println("Reading case: Fuel: $(curFuel), Range: $(curRange)")
        # Try to load and extract the design parameters
        try
            ac = quickload_aircraft(joinpath(model_dir,model_file))
            designParam, _  = ExtractDes(ac, save_dir, case_id; flg_save=false)
            push!(desParamRows, (; case = case_id, designParam...))
        catch err
            println("$(model_file) not readable")
        end
    end
end

#### Save the comparison design parameters
df_desPara = DataFrame(desParamRows)
CSV.write(joinpath(save_dirSub, "$(save_name)_DesPara.csv"), df_desPara)

#### Also save the fixed constraints and fixed optimization bounds
bounds_opt = BoundsOpt()
constraints_opt = ConstraintsOpt()
# saving
function save_struct(datStr,filename)
    numfields = fieldcount(typeof(datStr))
    lenfields = length(getfield(datStr,1))
    namfields = collect(fieldnames(typeof(datStr)))
    data = Matrix{Float64}(undef, lenfields, numfields)
    for i = 1:numfields
        for j = 1:lenfields
            data[j,i] = getfield(datStr,i)[j]
        end
    end
    df = DataFrame(data, Symbol.(namfields))
    CSV.write(filename, df)
    return df
end
save_struct(bounds_opt, joinpath(save_dirSub, "$(save_name)_OptBounds.csv"))
save_struct(constraints_opt, joinpath(save_dirSub, "$(save_name)_Constraints.csv"))