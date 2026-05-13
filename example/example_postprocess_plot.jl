"""
This script compare the design point performance between ethanol and jet
Adaptive design parameter bounds across different cases
Fixed and default constraints across different cases
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(@__DIR__, "postprocess.jl"))
using .PostProcess: ExtractDes
using .PostProcess.OptimizeRangeFuel: BoundsOpt, ConstraintsOpt
include(joinpath(@__DIR__, "optimize_rangefuel.jl"))
using .OptimizeRangeFuel.extract_opt_para

#### Setup IO
model_dir    = "ModelSavedTest"
model_prefix = "acOptimized_BatOpt" #Frontal key name for the models(FuelRange)
save_dir     = "ModelProcessed" #Outer Directory for saving models
save_prefix  = "Eth_Jet_Comparison" #The total save name for jet and ethanol comparison
mkpath(model_dir)
mkpath(save_dir)
save_dirSub = joinpath(save_dir,save_prefix) #Sub-directory to save comparison data
mkpath(save_dirSub)

#### Missions Parameters to Compare
Fuels = ["Eth","Jet"] #These corresponding to the model file name
Ranges = collect(300:100:3000) #Prefixex+Fuel+Range.jld2
flgSaveIndividual = false #Whether to save an output for individual case

#### Load the default constraints which is believed at this point to be common to call cases
constraints_opt = ConstraintsOpt()

#### Extract data for each case
paraSet = [] # Each element is one fuel
boundSet = [] # Each element is one fuel
optParSet = []
for (i,curFuel) in enumerate(Fuels)
    paraSetSub = [] # Each element is one range
    boundSetSub = [] # Each element is one range
    optParSetSub = []
    for (j,curRange) in enumerate(Ranges)
        caseCur = "$(curFuel)$(round(Int,curRange))"
        println("Attempt to read case: Fuel type: $(curFuel), Range: $(curRange)")
        try
            # Load the aircraft model
            ac = quickload_aircraft(joinpath(model_dir,"$(model_prefix)$(caseCur).jld2"))
            # Load the parameter bounds
            bdcsv = CSV.read(joinpath(model_dir,"$(model_prefix)$(caseCur)_BoundLocal.csv"), DataFrame)
            bounds_opt = BoundsOpt()
            for fName in fieldnames(typeof(bounds_opt))
                fNameSym = Symbol(fName)
                setfield!(bounds_opt, fName, (Float64(bdcsv[1, fNameSym]), Float64(bdcsv[2, fNameSym]), Float64(bdcsv[3, fNameSym])))
            end
            # Extract the design parameters and optionally save individual case
            design_para = ExtractDes(ac, save_dir, "$(model_prefix)$(caseCur)"; flg_save=flgSaveIndividual,
                          bounds_opt=bounds_opt, constraints_opt=constraints_opt)
            # Extract the subset of design parameters that have bounds corresponds to
            optimi_para = extract_opt_para(ac) #vector{Float64} Those parameters having bounds. Same order as the bound
            # Append the data
            push!(paraSetSub, (; range_id = string(round(Int,curRange)), fuel_id = curFuel, design_para...))
            push!(boundSetSub, (bounds_opt))
            push!(optParSetSub, (optimi_para))
        catch err
            println("Reading case failed for Fuel type: $(curFuel), Range: $(curRange)")
        end
    end
    push!(paraSet, (paraSetSub))
    push!(boundSet, (boundSetSub))
    push!(optParSet, (optParSetSub))
end