"""
This script take a design point and sweep an off-design payload-range envelope around it 
    - Compute payload range envelope
    - Compute R1 R2 R3 missions
    - Works for both OAG and R1S data
"""

using TASOPT
using DataFrames, CSV
using Glob
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/UtilityFunctions/offdesign.jl"))
using .PRD: off_design_PRD, findR1R2R3

#### Load directory
model_dir       = joinpath(__TASOPTroot__,"../example/ModelSaved")
model_folder    = "Opti_Eth_NoACT_OAG_MI_Tail_V3_"
name_front      = "Opti_Eth_NoACT_OAG_Ml_*" #Use "*" to represent outside of pattern region
name_back       = "" #led by range number, and followed by ".jld2"
ranges_or_seats = collect(50:10:240) # matching section must be "_(range number)" or else 50.jld2 cannot be distinguished from 150.jld2
                                     # this has to be from small to large for bisection bracketing 
#### Offdesign parameter
ranges_off_nmi  = collect(300.0:100.0:8000.0) #Regions of off-design ranges to test [nmi]
fuel_idx        = 32       #Eth: 32 , Jet: 24
rho_fuel        = 789.0    #Eth: 789.0 , Jet: 817.0 #kg/m3
hvap_fuel       = 918187.9 #Eth: 918187.9 , Jet: 358694.0 #J/kg

#### Save parameter (Save inside model_dir)
save_folder     = "Opti_Eth_NoACT_OAG_MI_Tail_V3_PRD"

#### Create directory
load_dir        = joinpath(model_dir,model_folder)
save_dir        = joinpath(model_dir,save_folder)

#### Setup job-array for ORCD cluster run
@assert length(ARGS) >= 2 "Usage: julia xxx.jl task_id num_tasks"
task_id = parse(Int,ARGS[1]) #Current task id
num_tasks = parse(Int,ARGS[2]) #Total number of tasks
@assert 1 <= task_id <= num_tasks

#### Extract data for each case
for i in task_id:num_tasks:length(ranges_or_seats)
    #### Load the baseline aircraft model at design point
    pattern = name_front * "_$(ranges_or_seats[i])" * name_back * ".jld2"
    matches = glob(pattern, load_dir)
    @assert length(matches)==1 "Expected exactly 1 match for $(pattern) in $(load_dir), found $(length(matches)): $(matches)"
    ac = quickload_aircraft(matches[1])
    matched_name = splitext(basename(matches[1]))[1]

    #### Setup sub-directory for saving
    name_save_cur = matched_name*"_PRD"
    save_dir_cur = joinpath(save_dir,name_save_cur)
    mkpath(save_dir_cur)
    
    #### Compute the PRD Envelop
    out_off = off_design_PRD(ac, fuel_idx, rho_fuel, hvap_fuel, ranges_off_nmi; 
                             save_dir = save_dir_cur, save_name = name_save_cur,  flg_save_ac = true)
    println("For design range: $(ranges_or_seats[i]) find feasible offdesign range from $(out_off["range_nmi"][1]) to $(out_off["range_nmi"][end])")

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
    println("For design range: $(ranges_or_seats[i]) find R1, R2, R3 ranges of $(out_R1["range_nmi"][1]), $(out_R2["range_nmi"][1]), $(out_R3["range_nmi"][1])")
end