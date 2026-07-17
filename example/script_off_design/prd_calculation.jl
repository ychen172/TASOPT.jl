"""
This script take a design point and sweep an off-design payload-range envelope around it 
and also compute the R1 R2 R3 ranges requirements
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_for_offdesign/offdesign.jl"))
using .PRD: off_design_PRD, findR1R2R3

#### Setup IO
model_dir       = "ModelSaved"
input_prefix    = "Opti_Jet_NoACT_" #Output name will be appended from this input name
name_added      = "to_Eth_OffDes_"
Ranges_design   = collect(300:100:3000) #Must match

#### Offdesign parameter
Ranges_sweep    = Float64.(collect(300:100:8000)) #off design range to test [nmi]
fuel_idx        = 32       #Eth: 32 , Jet: 24
rho_fuel        = 789.0    #Eth: 789.0 , Jet: 817.0 #kg/m3
hvap_fuel       = 918187.9 #Eth: 918187.9 , Jet: 358694.0 #J/kg

#### Initialize new folder for saved aircraft model
save_dir = joinpath(model_dir, input_prefix*name_added)
mkpath(save_dir)

#### Extract data for each case
for (i,curRange) in enumerate(Ranges_design)
    #### File path setup
    des_ran_str = string(round(Int,curRange)) #String form of current rounded range
    design_file_path = joinpath(model_dir, input_prefix, input_prefix*des_ran_str*".jld2")
    
    #### Load the aircraft model at design point
    ac = quickload_aircraft(design_file_path)
    
    #### run off-design
    save_name_cur = input_prefix*name_added*des_ran_str*"_"
    save_dir_cur = joinpath(save_dir, save_name_cur)
    mkdir(save_dir_cur) # sub-sub directory for each design case
    
    #### Compute the PRD Envelop
    out_off = off_design_PRD(ac, fuel_idx, rho_fuel, hvap_fuel, Ranges_sweep; 
                             save_dir = save_dir_cur, save_name = save_name_cur,  flg_save_ac = true)
    println("For design range: $(des_ran_str) find feasible offdesign range from $(out_off["range_nmi"][1]) to $(out_off["range_nmi"][end])")

    #### Compute the R1 R2 R3 characteristic mission
    RLB = minimum(out_off["range_nmi"])
    R_max_fea = maximum(out_off["range_nmi"]) #Maximum feasible range found from rough sweeping
    idx_min_infea = findfirst(x->x>R_max_fea, Ranges_sweep)
    RUB = Ranges_sweep[idx_min_infea] #Smallest infeasible range [nmi]
    out_R1 = findR1R2R3(:R1, RLB, RUB, ac, fuel_idx, rho_fuel, hvap_fuel;
                        flg_save_ac = true, save_name = save_name_cur*"R1_", save_dir = save_dir_cur,
                        epsRange = 1e-4, epsWpay = 1e-8)
    out_R2 = findR1R2R3(:R2, RLB, RUB, ac, fuel_idx, rho_fuel, hvap_fuel;
                        flg_save_ac = true, save_name = save_name_cur*"R2_", save_dir = save_dir_cur,
                        epsRange = 1e-4, epsWpay = 1e-8)
    out_R3 = findR1R2R3(:R3, RLB, RUB, ac, fuel_idx, rho_fuel, hvap_fuel;
                        flg_save_ac = true, save_name = save_name_cur*"R3_", save_dir = save_dir_cur,
                        epsRange = 1e-4, epsWpay = 1e-8)
    println("For design range: $(des_ran_str) find R1, R2, R3 ranges of $(out_R1["range_nmi"][1]), $(out_R2["range_nmi"][1]), $(out_R3["range_nmi"][1])")
end