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

#### Extract data for each scenerio
for (i, keyword_cur) in enumerate(case_keywords)
    model_dir_sub = joinpath(model_dir, keyword_cur)
    range_lst_sub = [] #(nmi)
    PFEI_lst_sub = [] #(J/J)
    # Extract data for each design case
    for (j, des_range_cur) in enumerate(des_ranges)
        model_dir_sub_sub = joinpath(model_dir_sub, keyword_cur*"$(round(Int,des_range_cur))")
        range_lst_sub_sub = [] #(nmi)
        PFEI_lst_sub_sub = [] #(J/J)
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
                push!(range_lst_sub_sub, range_cur)
                push!(PFEI_lst_sub_sub, PFEI_cur)
            catch
                nothing
            end
        end
        push!(range_lst_sub,range_lst_sub_sub)
        push!(PFEI_lst_sub,PFEI_lst_sub_sub)
    end
    push!(range_lst,range_lst_sub)
    push!(PFEI_lst,PFEI_lst_sub)
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