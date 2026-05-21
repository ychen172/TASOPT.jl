"""
This script compare the off-design performance between cases in terms of PFEI
and flight range
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
using Plots

#### Setup IO
# I
case_keywords = ["off_designJet", "jetfuel_match_payload", "jetfuel_to_ethanolJet"]
case_names    = ["Base", "Match", "Retrofit"]
model_dir     = "ModelSaved"
des_ranges    = [3000] #design range to compare (nmi)
# O
save_dir      = "ModelProcessed"
save_name     = "compare" #sub_folder will be created
# Test conditions
offdes_ranges = float.(collect(0:100:8000)) # (nmi)

#### Save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
range_lst = [] #(nmi)
PFEI_lst = [] #(J/J)
massEmp_lst = [] #(Ton)
voluFuel_lst = [] #(m3)
voluFuelMax_lst = [] #(m3)

#### Extract data for each scenerio
for (i, keyword_cur) in enumerate(case_keywords)
    model_dir_sub = joinpath(model_dir, keyword_cur)
    range_lst_sub = [] #(nmi)
    PFEI_lst_sub = [] #(J/J)
    massEmp_lst_sub = [] #(Ton)
    voluFuel_lst_sub = [] #(m3)
    voluFuelMax_lst_sub = [] #(m3)
    # Extract data for each design case
    for (j, des_range_cur) in enumerate(des_ranges)
        model_dir_sub_sub = joinpath(model_dir_sub, keyword_cur*"$(round(Int,des_range_cur))")
        range_lst_sub_sub = [] #(nmi)
        PFEI_lst_sub_sub = [] #(J/J)
        massEmp_lst_sub_sub = [] #(Ton)
        voluFuel_lst_sub_sub = [] #(m3)
        voluFuelMax_lst_sub_sub = [] #(m3)
        ## extract for each off-design case
        for (k, offdes_range_cur) in enumerate(offdes_ranges)
            model_dir_sub_sub_sub = joinpath(model_dir_sub_sub, keyword_cur*"$(round(Int,des_range_cur))_$(round(Int,offdes_range_cur)).jld2")
            ## reading
            println("Attempt to read from $(model_dir_sub_sub_sub)")
            try
                ac_cur = quickload_aircraft(model_dir_sub_sub_sub)
                ## extract data
                range_cur = ac_cur.parm[imRange,2] / 1852.0 #(nmi)
                PFEI_cur = ac_cur.parm[imPFEI, 2] #(J/J)
                
                ## Check on masses
                massTO = ac_cur.parm[imWTO,2]/gee/1000.0 #Takeoff mass (Ton)
                massFuelTot = ac_cur.parm[imWfuel,2]/gee/1000.0 #Fuel mass (Ton)(Include reserved and burned)
                massPayload = ac_cur.parm[imWpay, 2]/gee/1000.0 #(Ton)
                massEmpty = massTO - massFuelTot - massPayload #(Ton) empty weight
                ## Check on volume
                rhoFuel = ac_cur.parg[igrhofuel] #kg/m3
                volFuel = massFuelTot * 1000.0 / rhoFuel #m3
                volFuelMax = ac_cur.parg[igWfmax] / gee / rhoFuel #m3 (The design mission fuel mass might be different from the maximum fuel mass with off-design fuel density)
                
                ## store
                push!(range_lst_sub_sub, range_cur)
                push!(PFEI_lst_sub_sub, PFEI_cur)
                push!(massEmp_lst_sub_sub, massEmpty)
                push!(voluFuel_lst_sub_sub, volFuel)
                push!(voluFuelMax_lst_sub_sub, volFuelMax)
            catch
                nothing
            end
        end
        push!(range_lst_sub,range_lst_sub_sub)
        push!(PFEI_lst_sub,PFEI_lst_sub_sub)
        push!(massEmp_lst_sub,massEmp_lst_sub_sub)
        push!(voluFuel_lst_sub,voluFuel_lst_sub_sub)
        push!(voluFuelMax_lst_sub,voluFuelMax_lst_sub_sub)
    end
    push!(range_lst,range_lst_sub)
    push!(PFEI_lst,PFEI_lst_sub)
    push!(massEmp_lst,massEmp_lst_sub)
    push!(voluFuel_lst,voluFuel_lst_sub)
    push!(voluFuelMax_lst,voluFuelMax_lst_sub)
end

#### Plotting
# PFEI
plot_PFEI = plot(xlabel="Range (nmi)", ylabel="PFEI (J/J)", dpi=800)
plot!(plot_PFEI, xlims=(100,1600),ylims=(0.6, 1.0))
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_PFEI, range_lst[i][j], PFEI_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
    end
end
savefig(plot_PFEI, joinpath(save_dir_sub, "PFEI_comparison.png"))

# Empty weight
plot_mEmpty = plot(xlabel="Range (nmi)", ylabel="Empty Mass (Ton)", dpi=800)
plot!(plot_mEmpty, ylims=(0.0, massEmp_lst[1][1][1]+30.0))
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_mEmpty, range_lst[i][j], massEmp_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
    end
end
savefig(plot_mEmpty, joinpath(save_dir_sub, "mass_empty_comparison.png"))

# Fuel volume
plot_volfuel = plot(xlabel="Range (nmi)", ylabel="Fuel Volume (m3)", dpi=800)
# plot!(plot_volfuel, ylims=(0.0, massEmp_lst[1][1][1]+30.0))
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_volfuel, range_lst[i][j], voluFuel_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
    end
end
# for (i, keyword_cur) in enumerate(case_keywords)
#     for (j, des_range_cur) in enumerate(des_ranges)
#         plot!(plot_volfuel, range_lst[i][j], voluFuelMax_lst[i][j], marker=:none, lw=0.5, label=label=case_names[i]*"_$(round(Int,des_ranges[j]))")
#     end
# end
savefig(plot_volfuel, joinpath(save_dir_sub, "volume_fuel_comparison.png"))