"""
This script compare the R1 R2 R3 missions performance between retrofit and specially sized aircraft
Only compare PFEI from the csv log file
"""

using CSV
using DataFrames
using Plots
using TASOPT
include(__TASOPTindices__)

#### Save folder
save_dir      = "../ModelProcessed"
save_name     = "Design_versus_retrofit"
# Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Data path
R1_Ret_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R1.csv")
R2_Ret_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R2.csv")
R3_Ret_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R3.csv")
R1_Des_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_OffDes_/Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_OffDes_R1.csv")
R2_Des_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_OffDes_/Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_OffDes_R2.csv")
R3_Des_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_OffDes_/Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_OffDes_R3.csv")

#### Read in the R1 R2 R3 requirements for retrofit
df_R1 = CSV.read(R1_Ret_csv_path, DataFrame)
df_R2 = CSV.read(R2_Ret_csv_path, DataFrame)
df_R3 = CSV.read(R3_Ret_csv_path, DataFrame)
ran_des_ret_nmi  = Vector{Float64}(df_R1.design_range_nmi)
wei_pay_ret_N_R1 = Vector{Float64}(df_R1.payload_weight_N)
ran_nmi_ret_R1   = Vector{Float64}(df_R1.range_nmi)
wei_pay_ret_N_R2 = Vector{Float64}(df_R2.payload_weight_N)
ran_nmi_ret_R2   = Vector{Float64}(df_R2.range_nmi)
wei_pay_ret_N_R3 = Vector{Float64}(df_R3.payload_weight_N)
ran_nmi_ret_R3   = Vector{Float64}(df_R3.range_nmi)
#
PFEI_ret_R1      = Vector{Float64}(df_R1.PFEI_JJ)
PFEI_ret_R2      = Vector{Float64}(df_R2.PFEI_JJ)
PFEI_ret_R3      = Vector{Float64}(df_R3.PFEI_JJ)

#### Read in the R1 R2 R3 requirements for design
df_R1 = CSV.read(R1_Des_csv_path, DataFrame)
df_R2 = CSV.read(R2_Des_csv_path, DataFrame)
df_R3 = CSV.read(R3_Des_csv_path, DataFrame)
ran_des_des_nmi  = Vector{Float64}(df_R1.design_range_nmi)
wei_pay_des_N_R1 = Vector{Float64}(df_R1.payload_weight_N)
ran_nmi_des_R1   = Vector{Float64}(df_R1.range_nmi)
wei_pay_des_N_R2 = Vector{Float64}(df_R2.payload_weight_N)
ran_nmi_des_R2   = Vector{Float64}(df_R2.range_nmi)
wei_pay_des_N_R3 = Vector{Float64}(df_R3.payload_weight_N)
ran_nmi_des_R3   = Vector{Float64}(df_R3.range_nmi)
#
PFEI_des_R1      = Vector{Float64}(df_R1.PFEI_JJ)
PFEI_des_R2      = Vector{Float64}(df_R2.PFEI_JJ)
PFEI_des_R3      = Vector{Float64}(df_R3.PFEI_JJ)


#### Plotting
markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5]
linestyles = [:solid, :dash, :dot, :dashdot]

p = plot(xlabel="Design Range (nmi)", ylabel="Off Design Range (nmi)", dpi=800)
plot!(p, ran_des_ret_nmi, ran_nmi_ret_R1, marker=markers[1], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R1")
plot!(p, ran_des_ret_nmi, ran_nmi_ret_R2, marker=markers[2], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R2")
plot!(p, ran_des_ret_nmi, ran_nmi_ret_R3, marker=markers[3], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R3")
plot!(p, ran_des_des_nmi, ran_nmi_des_R1, marker=markers[4], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R1")
plot!(p, ran_des_des_nmi, ran_nmi_des_R2, marker=markers[5], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R2")
plot!(p, ran_des_des_nmi, ran_nmi_des_R3, marker=markers[6], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R3")
savefig(p, joinpath(save_dir_sub, "Range_Requirements.png"))

p = plot(xlabel="Design Range (nmi)", ylabel="Payload Weight (N)", dpi=800)
plot!(p, ran_des_ret_nmi, wei_pay_ret_N_R1, marker=markers[1], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R1")
plot!(p, ran_des_ret_nmi, wei_pay_ret_N_R2, marker=markers[2], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R2")
plot!(p, ran_des_ret_nmi, wei_pay_ret_N_R3, marker=markers[3], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R3")
plot!(p, ran_des_des_nmi, wei_pay_des_N_R1, marker=markers[4], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R1")
plot!(p, ran_des_des_nmi, wei_pay_des_N_R2, marker=markers[5], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R2")
plot!(p, ran_des_des_nmi, wei_pay_des_N_R3, marker=markers[6], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R3")
savefig(p, joinpath(save_dir_sub, "Payload_Requirements.png"))

p = plot(xlabel="Off Design Range (nmi)", ylabel="PFEI (J/J)", dpi=800)
plot!(p, ran_nmi_ret_R1, PFEI_ret_R1, marker=markers[1], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R1")
plot!(p, ran_nmi_ret_R2, PFEI_ret_R2, marker=markers[2], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R2")
plot!(p, ran_nmi_des_R1, PFEI_des_R1, marker=markers[4], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R1")
plot!(p, ran_nmi_des_R2, PFEI_des_R2, marker=markers[5], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R2")
savefig(p, joinpath(save_dir_sub, "PFEI_R1R2.png"))

p = plot(xlabel="Off Design Range (nmi)", ylabel="PFEI (J/J)", dpi=800)
plot!(p, ran_nmi_ret_R3, PFEI_ret_R3, marker=markers[3], linestyle=linestyles[1], lw=2, markerstrokewidth=0, label="Retrofit R3")
plot!(p, ran_nmi_des_R3, PFEI_des_R3, marker=markers[6], linestyle=linestyles[3], lw=2, markerstrokewidth=0, label="Optimized R3")
savefig(p, joinpath(save_dir_sub, "PFEI_R3.png"))