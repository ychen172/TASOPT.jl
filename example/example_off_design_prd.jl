"""
This script take a design point and sweep an off-design payload-range envelope around it
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(@__DIR__, "offdesign.jl"))
using .PRD: off_design_PRD, off_design_R1R2

#### File Name Creation
function modelFileName(model_prefix::String, curFuel::String, curRange::String)
    return "$(model_prefix)$(curFuel)$(curRange)"
    # return "$(curFuel)$(model_prefix)$(curRange)"
end

#### Setup IO
model_dir       = "ModelSaved"
input_prefix    = "acOptimized_BatOptJet" #Output name will be appended from this input name
save_key        = "Jet_" #addtional keyword for saving (Can be empty)
Ranges_design   = collect(300:100:3000) #Must match
#### Offdesign parameter
Ranges_sweep    = Float64.(collect(300:100:18000)) #off design range to test [nmi]
fuel_idx        = 24       #Eth: 32 , Jet: 24
rho_fuel        = 817.0    #Eth: 789.0 , Jet: 817.0 #kg/m3
hvap_fuel       = 358694.0 #Eth: 918187.9 , Jet: 358694.0 #J/kg

#### Initialize new folder for saved aircraft model
save_dir = joinpath(model_dir, input_prefix*"OffDes_"*save_key)
mkpath(save_dir)

#### Extract data for each case
for (i,curRange) in enumerate(Ranges_design)
    #### File path setup
    des_ran_str = string(round(Int,curRange)) #String form of current rounded range
    design_file_path = joinpath(model_dir, input_prefix, input_prefix*des_ran_str*".jld2")
    
    #### Load the aircraft model at design point
    ac = quickload_aircraft(design_file_path)
    
    #### run off-design
    # Create save folder
    save_name_cur = input_prefix*"_OffDes_"*save_key*des_ran_str*"_"
    save_dir_cur = joinpath(save_dir, save_name_cur)
    mkdir(save_dir_cur)
    
    # Compute the Envelop
    out_off = off_design_PRD(ac, fuel_idx, rho_fuel, hvap_fuel, Ranges_sweep; 
                             save_dir = save_dir_cur, save_name = save_name_cur,  flg_save_ac = true)
    println("For design range: $(des_ran_str) find feasible offdesign range from $(out_off["range_nmi"][1]) to $(out_off["range_nmi"][end])")

    #### Compute the R1 R2
    # determine the bounds for R1 and R2
    RLB = minimum(out_off["range_nmi"])
    RUB = maximum(out_off["range_nmi"]) #Lazy initialization 
    rangeBounds = [RLB,RUB,RLB,RUB]
    epsRange = 1e-6 #Need stronger convergence critria for lazy initialization
    out_R1, out_R2 = off_design_R1R2(ac, fuel_idx, rho_fuel, hvap_fuel, rangeBounds; 
                               epsRange = epsRange, save_dir = save_dir_cur, save_name = save_name_cur,  flg_save_ac = true)
    println("For design range: $(des_ran_str) find R1 and R2 ranges of $(out_R1["range_nmi"][1]) to $(out_R1["range_nmi"][1])")
end