"""
This script takes the R1 R2 R3 missions of another model, and offdesign the current model on the same mission.
"""

using Glob
using TASOPT
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_for_offdesign/offdesign.jl"))
using .PRD: off_design_specified!

#### Load model directory
read_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
model_key = "Opti_Eth_NoACT_V2_2_MatchR1R2R3_" #Model to read
miss_key  = "Opti_Jet_NoACT_V2_1_ULByEng_Eth_PRD_" #R1R2R3 key
save_key  = "Opti_Eth_NoACT_V2_2_OffdesR1R2R3_" #Save key
ran_des_nmi = collect(300.0:100.0:3000.0) #Must match with available model data

#### Offdesign parameter
fuel_idx  = 32       #Eth: 32 , Jet: 24
rho_fuel  = 789.0    #Eth: 789.0 , Jet: 817.0 #kg/m3
hvap_fuel = 918187.9 #Eth: 918187.9 , Jet: 358694.0 #J/kg

#### Setup job-array for ORCD cluster run
@assert length(ARGS) >= 2 "Usage: julia xxx.jl task_id num_tasks"
task_id = parse(Int,ARGS[1]) #Current task id
num_tasks = parse(Int,ARGS[2]) #Total number of tasks
@assert 1 <= task_id <= num_tasks

#### Extract data for each case
for i in task_id:num_tasks:length(ran_des_nmi)
    #### Load aircraft model
    model_dir = joinpath(read_dir,model_key,"$(model_key)$(round(Int,ran_des_nmi[i])).jld2")
    ac = quickload_aircraft(model_dir)
        
    #### Load mission requirements
    miss_dir = joinpath(read_dir,miss_key,"$(miss_key)$(round(Int,ran_des_nmi[i]))_")
    R1_key   = "$(miss_key)$(round(Int,ran_des_nmi[i]))_R1_*.jld2"
    R1_key   = only(glob(R1_key, miss_dir)) #Find the only matched file
    miss     = quickload_aircraft(R1_key)
    Range_R1 = miss.parm[imRange, 2]/1852.0 #[nmi]
    Wpay_R1  = miss.parm[imWpay, 2] #[N]
    #
    R2_key   = "$(miss_key)$(round(Int,ran_des_nmi[i]))_R2_*.jld2"
    R2_key   = only(glob(R2_key, miss_dir)) #Find the only matched file
    miss     = quickload_aircraft(R2_key)
    Range_R2 = miss.parm[imRange, 2]/1852.0 #[nmi]
    Wpay_R2  = miss.parm[imWpay, 2] #[N]
    #
    R3_key   = "$(miss_key)$(round(Int,ran_des_nmi[i]))_R3_*.jld2"
    R3_key   = only(glob(R3_key, miss_dir)) #Find the only matched file
    miss     = quickload_aircraft(R3_key)
    Range_R3 = miss.parm[imRange, 2]/1852.0 #[nmi]
    Wpay_R3  = miss.parm[imWpay, 2] #[N]

    #### Setup sub-directory for saving
    save_dir_cur = joinpath(read_dir,save_key,"$(save_key)$(round(Int,ran_des_nmi[i]))_")
    mkpath(save_dir_cur)

    #### Perform off-design on the three missions
    off_design_specified!(ac,fuel_idx,rho_fuel,hvap_fuel,[Range_R1],[Wpay_R1];
                          mod_ac_inplace=false,save_model=true,save_dir=save_dir_cur,save_name="$(save_key)$(round(Int,ran_des_nmi[i]))_"*"R1_")
    off_design_specified!(ac,fuel_idx,rho_fuel,hvap_fuel,[Range_R2],[Wpay_R2];
                          mod_ac_inplace=false,save_model=true,save_dir=save_dir_cur,save_name="$(save_key)$(round(Int,ran_des_nmi[i]))_"*"R2_")
    off_design_specified!(ac,fuel_idx,rho_fuel,hvap_fuel,[Range_R3],[Wpay_R3];
                          mod_ac_inplace=false,save_model=true,save_dir=save_dir_cur,save_name="$(save_key)$(round(Int,ran_des_nmi[i]))_"*"R3_")
end