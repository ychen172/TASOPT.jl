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

#### file keys and labels
read_dir = "../ModelSaved/"
keys = ["R1R2R3_Jet_NoACT_to_Eth_","Opti_Jet_NoACT_to_Eth_for_Eth_MatR1R2R3_","Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth_MatR1R2R3_"]
labels = ["Retrofit", "Optimized Retrofit", "Pure Optimization"]
ran_nmi_des_R1_lst = []
ran_nmi_des_R2_lst = []
ran_nmi_des_R3_lst = []
ran_nmi_R1_lst = []
ran_nmi_R2_lst = []
ran_nmi_R3_lst = []
wei_pay_N_R1_lst = []
wei_pay_N_R2_lst = []
wei_pay_N_R3_lst = []
PFEI_JJ_R1_lst = []
PFEI_JJ_R2_lst = []
PFEI_JJ_R3_lst = []
for key_cur in keys
    dir_cur = joinpath(read_dir,key_cur,key_cur)
    df_R1 = CSV.read(dir_cur*"R1.csv", DataFrame)
    df_R2 = CSV.read(dir_cur*"R2.csv", DataFrame)
    df_R3 = CSV.read(dir_cur*"R3.csv", DataFrame)
    push!(ran_nmi_des_R1_lst, Vector{Float64}(df_R1.design_range_nmi))
    push!(ran_nmi_des_R2_lst, Vector{Float64}(df_R2.design_range_nmi))
    push!(ran_nmi_des_R3_lst, Vector{Float64}(df_R3.design_range_nmi))
    push!(ran_nmi_R1_lst, Vector{Float64}(df_R1.range_nmi))
    push!(ran_nmi_R2_lst, Vector{Float64}(df_R2.range_nmi))
    push!(ran_nmi_R3_lst, Vector{Float64}(df_R3.range_nmi))
    push!(wei_pay_N_R1_lst, Vector{Float64}(df_R1.payload_weight_N))
    push!(wei_pay_N_R2_lst, Vector{Float64}(df_R2.payload_weight_N))
    push!(wei_pay_N_R3_lst, Vector{Float64}(df_R3.payload_weight_N))
    push!(PFEI_JJ_R1_lst, Vector{Float64}(df_R1.PFEI_JJ))
    push!(PFEI_JJ_R2_lst, Vector{Float64}(df_R2.PFEI_JJ))
    push!(PFEI_JJ_R3_lst, Vector{Float64}(df_R3.PFEI_JJ))
end

#### Plotting
markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :star5]
linestyles = [:solid, :dash, :dot, :dashdot]

p = plot(xlabel="Design Range (nmi)", ylabel="Retrofit Range (nmi)", dpi=800)
for i in eachindex(ran_nmi_des_R1_lst)
    plot!(p, ran_nmi_des_R1_lst[i], ran_nmi_R1_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i]*"_R1")
    plot!(p, ran_nmi_des_R2_lst[i], ran_nmi_R2_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i]*"_R2")
    plot!(p, ran_nmi_des_R3_lst[i], ran_nmi_R3_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i]*"_R3")
end
savefig(p, joinpath(save_dir_sub, "Design_Range_to_Retrofit_Range.png"))

p = plot(xlabel="Design Range (nmi)", ylabel="Payload Weight (N)", dpi=800)
for i in eachindex(ran_nmi_des_R1_lst)
    plot!(p, ran_nmi_des_R1_lst[i], wei_pay_N_R1_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i]*"_R1")
    plot!(p, ran_nmi_des_R2_lst[i], wei_pay_N_R2_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i]*"_R2")
    plot!(p, ran_nmi_des_R3_lst[i], wei_pay_N_R3_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i]*"_R3")
end
savefig(p, joinpath(save_dir_sub, "Design_Wpay_to_Retrofit_Wpay.png"))

p = plot(xlabel="Retrofit R1 (nmi)", ylabel="R1 PFEI (J/J)", dpi=800)
for i in eachindex(ran_nmi_R1_lst)
    plot!(p, ran_nmi_R1_lst[i], PFEI_JJ_R1_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i])
end
savefig(p, joinpath(save_dir_sub, "R1_PFEI.png"))

p = plot(xlabel="Retrofit R2 (nmi)", ylabel="R2 PFEI (J/J)", dpi=800)
for i in eachindex(ran_nmi_R2_lst)
    plot!(p, ran_nmi_R2_lst[i], PFEI_JJ_R2_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i])
end
savefig(p, joinpath(save_dir_sub, "R2_PFEI.png"))

p = plot(xlabel="Retrofit R3 (nmi)", ylabel="R3 PFEI (J/J)", dpi=800)
for i in eachindex(ran_nmi_R3_lst)
    plot!(p, ran_nmi_R3_lst[i], PFEI_JJ_R3_lst[i], marker=markers[i], linestyle=linestyles[i], lw=2, markerstrokewidth=0, label=labels[i])
end
savefig(p, joinpath(save_dir_sub, "R3_PFEI.png"))