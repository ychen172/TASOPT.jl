"""
This script take a design point and sweep an off-design payload-range envelope around it 
and also compute the R1 R2 R3 ranges requirements
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_for_offdesign/offdesign.jl"))
using .PRD: off_design_PRD, findR1R2R3

#### Load directory
model_dir       = joinpath(__TASOPTroot__,"../example/ModelSaved")
input_prefix    = "Opti_Jet_NoACT_V3_" #Output name will be appended from this input name
ranges_des_nmi  = collect(300.0:100.0:3000.0) #Must match with available model data

#### Offdesign parameter
ranges_off_nmi  = collect(300.0:100.0:8000.0) #off design range to test [nmi]
fuel_idx        = 24       #Eth: 32 , Jet: 24
rho_fuel        = 817.0    #Eth: 789.0 , Jet: 817.0 #kg/m3
hvap_fuel       = 358694.0 #Eth: 918187.9 , Jet: 358694.0 #J/kg

#### Save parameter
name_save       = "Opti_Jet_NoACT_V3_PRD_"

#### Create directory
load_dir        = joinpath(model_dir,input_prefix,input_prefix)
save_dir        = joinpath(model_dir,name_save)

#### Setup job-array for ORCD cluster run
@assert length(ARGS) >= 2 "Usage: julia xxx.jl task_id num_tasks"
task_id = parse(Int,ARGS[1]) #Current task id
num_tasks = parse(Int,ARGS[2]) #Total number of tasks
@assert 1 <= task_id <= num_tasks

#### Extract data for each case
for i in task_id:num_tasks:length(ranges_des_nmi)
    #### Load the baseline aircraft model at design point
    ac = quickload_aircraft(load_dir*"$(round(Int , ranges_des_nmi[i])).jld2")
    
    #### Setup sub-directory for saving
    name_save_cur = name_save*"$(round(Int , ranges_des_nmi[i]))_"
    save_dir_cur = joinpath(save_dir,name_save_cur)
    mkpath(save_dir_cur)
    
    #### Compute the PRD Envelop
    out_off = off_design_PRD(ac, fuel_idx, rho_fuel, hvap_fuel, ranges_off_nmi; 
                             save_dir = save_dir_cur, save_name = name_save_cur,  flg_save_ac = true)
    println("For design range: $(ranges_des_nmi[i]) find feasible offdesign range from $(out_off["range_nmi"][1]) to $(out_off["range_nmi"][end])")

    #### Compute the R1 R2 R3 characteristic mission
    RLB = minimum(out_off["range_nmi"])
    RUB = ranges_off_nmi[findfirst(x -> x>maximum(out_off["range_nmi"]), ranges_off_nmi)] #Shortest infeasible range [nmi]
    out_R1 = findR1R2R3(:R1, RLB, RUB, ac, fuel_idx, rho_fuel, hvap_fuel;
                        flg_save_ac = true, save_name = name_save_cur*"R1_", save_dir = save_dir_cur,
                        epsRange = 1e-4, epsWpay = 1e-8)
    out_R2 = findR1R2R3(:R2, RLB, RUB, ac, fuel_idx, rho_fuel, hvap_fuel;
                        flg_save_ac = true, save_name = name_save_cur*"R2_", save_dir = save_dir_cur,
                        epsRange = 1e-4, epsWpay = 1e-8)
    out_R3 = findR1R2R3(:R3, RLB, RUB, ac, fuel_idx, rho_fuel, hvap_fuel;
                        flg_save_ac = true, save_name = name_save_cur*"R3_", save_dir = save_dir_cur,
                        epsRange = 1e-4, epsWpay = 1e-8)
    println("For design range: $(ranges_des_nmi[i]) find R1, R2, R3 ranges of $(out_R1["range_nmi"][1]), $(out_R2["range_nmi"][1]), $(out_R3["range_nmi"][1])")
end