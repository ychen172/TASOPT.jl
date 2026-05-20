"""
This script take a design point and sweep an off-design payload-range envelope around it
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(@__DIR__, "offdesign.jl"))
using .PRD: off_design_PRD

#### File Name Creation
function modelFileName(model_prefix::String, curFuel::String, curRange::String)
    return "$(model_prefix)$(curFuel)$(curRange)"
    # return "$(curFuel)$(model_prefix)$(curRange)"
end

#### Setup IO
model_dir    = "ModelSaved"
model_prefix = "acOptimized_BatOpt" #Frontal key name for the models(FuelRange)
save_dir_prefix = "ModelSaved" #Outer Directory for saving the offdesign models
save_prefix  = "off_design" #added with fuel and designed_range to create a sub folder name. Each then contain sub-model for off-design sweep
mkpath(model_dir)
Fuels = ["Eth"] #["Eth", "Jet"] #These corresponding to the model file name
Ranges_design = [300] #collect(300:100:3000) #Prefixex+Fuel+Range.jld2
#### Offdesign parameter
Ranges_sweep = Float64.(collect(100:100:18000)) #off design range to sweep [nmi]
fuel_idx = [32] #[32 , 24]
rho_fuel = [789.0] #[789.0 , 817.0] #kg/m3
hvap_fuel = [918187.9] #[918187.9 , 358694.0] #J/kg

#### Extract data for each case
for (i,curFuel) in enumerate(Fuels)
    for (j,curRange) in enumerate(Ranges_design)
        model_name = modelFileName(model_prefix,curFuel,string(round(Int,curRange)))
        println("Compute off design for: Fuel type: $(curFuel), Range: $(curRange)")
        try
            #### Load the aircraft model
            ac = quickload_aircraft(joinpath(model_dir,"$(model_name).jld2"))
            #### run off-design
            # Create a sub-folder as save Directory
            save_dir_sub = joinpath(save_dir_prefix,"$(save_prefix)$(curFuel)$(string(round(Int,curRange)))")
            mkpath(save_dir_sub)
            out_off = off_design_PRD(ac, fuel_idx[i], rho_fuel[i], hvap_fuel[i], Ranges_sweep; 
                      save_dir = save_dir_sub, save_name = "$(save_prefix)$(curFuel)$(string(round(Int,curRange)))_")
            println("Get off-design maximum ranges: $(out_off["range_nmi"])")
        catch err
            println(err)
            println("Reading case failed for Fuel type: $(curFuel), Range: $(curRange)")
        end
    end
end