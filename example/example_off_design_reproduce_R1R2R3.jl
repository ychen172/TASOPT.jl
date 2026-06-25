"""
This script takes an optimally sized ethanol aircraft and run the R1 R2 R3 missions fo a retroffited aircraft
"""

using CSV
using DataFrames
using TASOPT
include(__TASOPTindices__)
include(joinpath(@__DIR__,"utilities_for_offdesign/offdesign.jl"))
using .PRD:off_design_specified!

#### Optimization parameters
# CSV for R1 R2 R3 Missions requirements
R1_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R1.csv")
R2_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R2.csv")
R3_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R3.csv")
# Prefix for baseline aircraft parameters as a starting point
par_path_base_prefix = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_/Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_") #Also the prefixed for aircraft model
# Path to save the models from optimization
save_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/")
save_key = "Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_OffDes"
# 1.1) Optimization configuration parameters
max_iter_sizing = 150 # Maximum iterations for TASOPT sizing

#### Save the directory to save the optimized solution
save_dir_actual = joinpath(save_dir,save_key*"_") #Actual directory to save data
mkpath(save_dir_actual)

#### Read in the R1 R2 R3 requirements
idx_start = 1
df_R1 = CSV.read(R1_csv_path, DataFrame)
df_R2 = CSV.read(R2_csv_path, DataFrame)
df_R3 = CSV.read(R3_csv_path, DataFrame)
ran_des_nmi  = Vector{Float64}(df_R1.design_range_nmi[idx_start:end])
ran_des_nmi  = round.(Int, ran_des_nmi) #Round the original jet aircraft design range just for name in saving
wei_pay_N_R1 = Vector{Float64}(df_R1.payload_weight_N[idx_start:end])
ran_nmi_R1   = Vector{Float64}(df_R1.range_nmi[idx_start:end])
wei_pay_N_R2 = Vector{Float64}(df_R2.payload_weight_N[idx_start:end])
ran_nmi_R2   = Vector{Float64}(df_R2.range_nmi[idx_start:end])
wei_pay_N_R3 = Vector{Float64}(df_R3.payload_weight_N[idx_start:end])
ran_nmi_R3   = Vector{Float64}(df_R3.range_nmi[idx_start:end])

#### Setup fuel
fuel_idx  = 32        
rho_fuel  = 789.0     
hvap_fuel = 918187.9

#### Run off design
out_R1_collect = Dict() #Collect range outputs
out_R2_collect = Dict()
out_R3_collect = Dict()
for (i,name_cur_range) in enumerate(ran_des_nmi)
    println("running $(name_cur_range)")
    # Load the baseline aircraft model as a starting guess
    ac = quickload_aircraft(par_path_base_prefix*"$(name_cur_range).jld2")
    
    # Design mission follows the baseline setup. Specify off-design missions requirements here
    range_off_des_nmi = [ran_nmi_R1[i], ran_nmi_R2[i], ran_nmi_R3[i]]
    wei_pay_off_des_N = [wei_pay_N_R1[i], wei_pay_N_R2[i], wei_pay_N_R3[i]]

    # Run the global local optimization process
    out = off_design_specified!(ac, fuel_idx, rho_fuel, hvap_fuel, range_off_des_nmi, wei_pay_off_des_N;
                                mod_ac_inplace=false, itermax=max_iter_sizing, save_model=true, 
                                save_dir = save_dir_actual, save_name = save_key*"_$(name_cur_range)_")

    # Collect the outputs
    if i == 1
        #
        out_R1_collect["design_range_nmi"] = [name_cur_range]
        out_R1_collect["range_nmi"] = [out.ran_nmi[1]]
        out_R1_collect["payload_weight_N"] = [out.wei_pay_N[1]]
        out_R1_collect["PFEI_JJ"] = [out.PFEI_JJ[1]]
        out_R1_collect["fuel_tank_frac"] = [out.fuel_tank_frac[1]]
        out_R1_collect["payload_frac"] = [out.payload_frac[1]]
        #
        out_R2_collect["design_range_nmi"] = [name_cur_range]
        out_R2_collect["range_nmi"] = [out.ran_nmi[2]]
        out_R2_collect["payload_weight_N"] = [out.wei_pay_N[2]]
        out_R2_collect["PFEI_JJ"] = [out.PFEI_JJ[2]]
        out_R2_collect["fuel_tank_frac"] = [out.fuel_tank_frac[2]]
        out_R2_collect["payload_frac"] = [out.payload_frac[2]]
        #
        out_R3_collect["design_range_nmi"] = [name_cur_range]
        out_R3_collect["range_nmi"] = [out.ran_nmi[3]]
        out_R3_collect["payload_weight_N"] = [out.wei_pay_N[3]]
        out_R3_collect["PFEI_JJ"] = [out.PFEI_JJ[3]]
        out_R3_collect["fuel_tank_frac"] = [out.fuel_tank_frac[3]]
        out_R3_collect["payload_frac"] = [out.payload_frac[3]]
    else
        #
        push!(out_R1_collect["design_range_nmi"], name_cur_range)
        push!(out_R1_collect["range_nmi"], out.ran_nmi[1])
        push!(out_R1_collect["payload_weight_N"], out.wei_pay_N[1])
        push!(out_R1_collect["PFEI_JJ"], out.PFEI_JJ[1])
        push!(out_R1_collect["fuel_tank_frac"], out.fuel_tank_frac[1])
        push!(out_R1_collect["payload_frac"], out.payload_frac[1])
        #
        push!(out_R2_collect["design_range_nmi"], name_cur_range)
        push!(out_R2_collect["range_nmi"], out.ran_nmi[2])
        push!(out_R2_collect["payload_weight_N"], out.wei_pay_N[2])
        push!(out_R2_collect["PFEI_JJ"], out.PFEI_JJ[2])
        push!(out_R2_collect["fuel_tank_frac"], out.fuel_tank_frac[2])
        push!(out_R2_collect["payload_frac"], out.payload_frac[2])
        #
        push!(out_R3_collect["design_range_nmi"], name_cur_range)
        push!(out_R3_collect["range_nmi"], out.ran_nmi[3])
        push!(out_R3_collect["payload_weight_N"], out.wei_pay_N[3])
        push!(out_R3_collect["PFEI_JJ"], out.PFEI_JJ[3])
        push!(out_R3_collect["fuel_tank_frac"], out.fuel_tank_frac[3])
        push!(out_R3_collect["payload_frac"], out.payload_frac[3])
    end
end

#### Print out the collect R1 R2 R3 data into csv a csv file
R1_df = DataFrame(; (Symbol(k) => v for (k,v) in out_R1_collect)...)
CSV.write(joinpath(save_dir_actual,save_key*"_R1.csv"), R1_df)
R2_df = DataFrame(; (Symbol(k) => v for (k,v) in out_R2_collect)...)
CSV.write(joinpath(save_dir_actual,save_key*"_R2.csv"), R2_df)
R3_df = DataFrame(; (Symbol(k) => v for (k,v) in out_R3_collect)...)
CSV.write(joinpath(save_dir_actual,save_key*"_R3.csv"), R3_df)