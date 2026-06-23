"""
This script take a series of design model and identify their R1, R2, R3 missions requirement.
Off-design allows the use of different fuel for retrofitting case.
"""

using TASOPT
include(__TASOPTindices__)
using DataFrames, CSV
include(joinpath(@__DIR__, "offdesign.jl"))
using .PRD: findR1R2R3

#### Setup IO
model_dir     = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_/Opti_Jet_NoACT_") #Need design range extension
save_dir_pre  = joinpath(__TASOPTroot__,"../example/ModelSaved") #A subfolder with save_key name will be created
save_key      = "R1R2R3_Jet_NoACT_to_Eth"
# Design ranges to test
ranges_design = collect(300:100:3000) #match with thr ranges in model_dir
# Off-design fuel type and range bound
range_LB      = 300.0 #[nmi] Lower bound range. Assured model convergence
range_UB      = 6000.0 #[nmi] Upper bound range. may fail the mission
fuel_idx      = 32       #Eth: 32 , Jet: 24
rho_fuel      = 789.0    #Eth: 789.0 , Jet: 817.0 #kg/m3
hvap_fuel     = 918187.9 #Eth: 918187.9 , Jet: 358694.0 #J/kg

#### Initialize new folder for saved aircraft model
save_dir = joinpath(save_dir_pre,save_key*"_") #Consistent naming between content and folder
mkpath(save_dir)

#### Initialize container to save ranges output into CSV
out_R1_collect = Dict() #Collect range outputs
out_R2_collect = Dict()
out_R3_collect = Dict()

#### Find R1 R2 R3 for each design mission
for (i,range_cur) in enumerate(ranges_design)
    println("Find R1, R2, R3 for $(range_cur) nmi design mission")
    # Read in the aircraft model
    ac = quickload_aircraft(model_dir*"$(round(Int,range_cur)).jld2")    
    
    # Find R1, R2, R3
    out_R1 = findR1R2R3(:R1, range_LB, range_UB, ac, fuel_idx, rho_fuel, hvap_fuel;
                          flg_save_ac = true, save_name = save_key*"_$(round(Int,range_cur))_R1_", save_dir = save_dir)
    out_R2 = findR1R2R3(:R2, range_LB, range_UB, ac, fuel_idx, rho_fuel, hvap_fuel;
                          flg_save_ac = true, save_name = save_key*"_$(round(Int,range_cur))_R2_", save_dir = save_dir)
    out_R3 = findR1R2R3(:R3, range_LB, range_UB, ac, fuel_idx, rho_fuel, hvap_fuel;
                          flg_save_ac = true, save_name = save_key*"_$(round(Int,range_cur))_R3_", save_dir = save_dir)
    (length(out_R1["payload_weight_N"])>0 && length(out_R2["payload_weight_N"])>0 && length(out_R3["payload_weight_N"])>0) ||
    throw(ErrorException("One of the len(R1)$(length(out_R1["payload_weight_N"])), len(R2)$(length(out_R2["payload_weight_N"])), len(R3)$(length(out_R3["payload_weight_N"])), finding cases has zero element (no solution found)"))

    # Collect the R1 R2 R3 Information
    if i == 1
        global out_R1_collect = out_R1
        global out_R2_collect = out_R2
        global out_R3_collect = out_R3
    else
        for (key, value) in out_R1
            push!(out_R1_collect[key], value[1])
        end
        for (key, value) in out_R2
            push!(out_R2_collect[key], value[1])
        end
        for (key, value) in out_R3
            push!(out_R3_collect[key], value[1])
        end
    end
end

#### Preprending the design range Information
out_R1_collect["design_range_nmi"] = ranges_design
out_R2_collect["design_range_nmi"] = ranges_design
out_R3_collect["design_range_nmi"] = ranges_design

#### Print out the collect R1 R2 R3 data into csv a csv file
R1_df = DataFrame(; (Symbol(k) => v for (k,v) in out_R1_collect)...)
CSV.write(joinpath(save_dir,save_key*"_R1.csv"), R1_df)
R2_df = DataFrame(; (Symbol(k) => v for (k,v) in out_R2_collect)...)
CSV.write(joinpath(save_dir,save_key*"_R2.csv"), R2_df)
R3_df = DataFrame(; (Symbol(k) => v for (k,v) in out_R3_collect)...)
CSV.write(joinpath(save_dir,save_key*"_R3.csv"), R3_df)