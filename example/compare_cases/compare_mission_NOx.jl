"""
This script compare the EINOx and total NOx emissions across multiple missions
    Reading from csv files with *SweepOutPutTot.csv as extension
    Assume phases: [1:5]C1,C2,C3,C4,C5 -> [6:7]R1,R2 -> [8:12]D1,D2,D3,D4,D5
"""

using DataFrames, CSV
using Plots
using TASOPT

#### Setup IO
dir_out = joinpath(__TASOPTroot__,"../example/CombustorOutSaved/")
key_out = ["Opti_Jet_NoACT_V2_",       "Opti_Eth_NoACT_V2_"]
dir_cyc = joinpath(__TASOPTroot__,"../example/CombustorCycleSaved/")
key_cyc = ["Opti_Jet_NoACT_V2_CycIn_", "Opti_Eth_NoACT_V2_CycIn_"]
nam_read = ["Jet Fuel",                 "Ethanol"]
ranges   = collect(300:100:3000)
range_expose = 400
# Output
save_dir = joinpath(__TASOPTroot__,"../example/ModelProcessed/")
save_nam = "NOx_Jet_vs_Eth_No_ACT_Des"
save_dir = joinpath(save_dir, save_nam)
mkpath(save_dir)

#### Initialization
fields = (:ranges_nmi, :NOx_Cli_kg, :NOx_Cru_kg, :NOx_Des_kg, :NOx_Tot_kg)
results = [Dict(field => Vector{Any}(undef, length(ranges)) for field in fields) for _ in eachindex(nam_read)] #NumKey * Fields * numRange
fields_expose = (:phases, :Tt4_K_Comp, :Tt4_K_Exp, :Tt3_K, :mdot_air_kgs)
results_expose = [Dict(field => Vector{Any}(undef, 12) for field in fields_expose) for _ in eachindex(nam_read)] #NumKey * Fields_expose * numPhases

#### Extract Emissions Data
for idx_key in eachindex(nam_read)
    for (idx_cur,ran_cur) in enumerate(ranges)
        # Reading the CSV output files and cycle files
        sub_dir_out = joinpath(dir_out, key_out[idx_key], key_out[idx_key]*"$(ran_cur)SweepOutPutTot.csv")
        sub_dir_cyc = joinpath(dir_cyc, key_cyc[idx_key], key_cyc[idx_key]*"$(ran_cur).csv")
        out = CSV.read(sub_dir_out, DataFrame)
        cyc = CSV.read(sub_dir_cyc, DataFrame)
        # NOx emissions indices at each phase
        phases = cyc[:, Symbol("Phase")]
        time_cur = cyc[:, Symbol("Time[s]")] #[s]
        EINOx_cur = out[:, Symbol(" EINOx(g/kg)")] #[g/kg]
        mdot_fuel_cur = out[:, Symbol(" mdot_f(kg/s)")] #[kg/s]
        mdot_NOx_cur = EINOx_cur .* mdot_fuel_cur ./ 1000.0 #[kg/s]
        # Calculate mission NOx emissions at different phases and total
        NOx_Cli_cur = sum(0.5 .* (mdot_NOx_cur[1:5][1:end-1]  .+ mdot_NOx_cur[1:5][2:end])  .* diff(time_cur[1:5])) #[kg]
        NOx_Cru_cur = sum(0.5 .* (mdot_NOx_cur[6:7][1:end-1]  .+ mdot_NOx_cur[6:7][2:end])  .* diff(time_cur[6:7]))
        NOx_Des_cur = sum(0.5 .* (mdot_NOx_cur[8:12][1:end-1] .+ mdot_NOx_cur[8:12][2:end]) .* diff(time_cur[8:12]))
        NOx_Tot_cur = NOx_Cli_cur + NOx_Cru_cur + NOx_Des_cur
        results[idx_key][:NOx_Cli_kg][idx_cur] = NOx_Cli_cur
        results[idx_key][:NOx_Cru_kg][idx_cur] = NOx_Cru_cur
        results[idx_key][:NOx_Des_kg][idx_cur] = NOx_Des_cur
        results[idx_key][:NOx_Tot_kg][idx_cur] = NOx_Tot_cur
        results[idx_key][:ranges_nmi][idx_cur] = ran_cur
        #### Extract the combustor efficiency plot for a certain range case
        if Int(ran_cur) == Int(range_expose)
            results_expose[idx_key][:phases] .= phases
            results_expose[idx_key][:Tt4_K_Comp] .= out[:, Symbol(" T4(K)")] #[K]
            results_expose[idx_key][:Tt4_K_Exp] .= cyc[:, Symbol("Tt4[R]")] .* (5.0/9.0) #[K]
            results_expose[idx_key][:Tt3_K] .= cyc[:, Symbol("Tt3[R]")] .* (5.0/9.0) #[K]
            results_expose[idx_key][:mdot_air_kgs] .= cyc[:, Symbol("W3[lbm/s]")] .* 0.45359237 #[kg/s]
        end
    end
end


#### Compare the NOx mission emissions
linestyles = repeat([:solid, :dash, :dot, :dashdot, :dashdotdot],3)
linecolors = repeat([:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :black, :cyan,
                     :magenta, :teal, :navy, :maroon, :olive, :gold, :coral, :turquoise, :lime, :indigo], 3) 
markers = [:square, :circle, :diamond, :pentagon]

p = plot(xlabel="Design Range [nmi]", ylabel="NOx Emissions [kg]", dpi=800)
for idx_key in eachindex(nam_read)
    plot!(p, results[idx_key][:ranges_nmi], results[idx_key][:NOx_Tot_kg], marker=markers[idx_key], linecolors=linecolors[1], markercolor=linecolors[1],  linestyles=linestyles[1], lw=2, markerstrokewidth=0, label="$(nam_read[idx_key]) Total NOx")
    plot!(p, results[idx_key][:ranges_nmi], results[idx_key][:NOx_Cli_kg], marker=markers[idx_key], linecolors=linecolors[2], markercolor=linecolors[2],  linestyles=linestyles[2], lw=2, markerstrokewidth=0, label="$(nam_read[idx_key]) Climb NOx")
    plot!(p, results[idx_key][:ranges_nmi], results[idx_key][:NOx_Cru_kg], marker=markers[idx_key], linecolors=linecolors[3], markercolor=linecolors[3],  linestyles=linestyles[3], lw=2, markerstrokewidth=0, label="$(nam_read[idx_key]) Cruise NOx")
    plot!(p, results[idx_key][:ranges_nmi], results[idx_key][:NOx_Des_kg], marker=markers[idx_key], linecolors=linecolors[4], markercolor=linecolors[4],  linestyles=linestyles[4], lw=2, markerstrokewidth=0, label="$(nam_read[idx_key]) Descent NOx")
end
savefig(p, joinpath(save_dir, "NOx_Comp.png"))


#### Plot out the individiual exposed missions data
for idx_key in eachindex(nam_read)
    filename = replace(nam_read[idx_key], " " => "_")
    x = eachindex(results_expose[idx_key][:phases])
    # Tt4 plot
    global p = plot(lw=2, xticks=(collect(x), results_expose[idx_key][:phases]), xlabel="Phase", ylabel="Tt4 [K]", dpi=800)
    plot!(p, x, results_expose[idx_key][:Tt4_K_Comp],label="Computed",marker=markers[1])
    plot!(p, x, results_expose[idx_key][:Tt4_K_Exp],label="Expected",marker=markers[1])
    savefig(p, joinpath(save_dir, "Tt4_Comp_$(filename)_$(range_expose).png"))

    # Tt3 Plot
    global p = plot(lw=2, xticks=(collect(x), results_expose[idx_key][:phases]), xlabel="Phase", ylabel="Tt3 [K]", dpi=800)
    plot!(p, x, results_expose[idx_key][:Tt3_K],marker=markers[1])
    savefig(p, joinpath(save_dir, "Tt3_Comp_$(filename)_$(range_expose).png"))

    # Air flow rate
    global p = plot(lw=2, xticks=(collect(x), results_expose[idx_key][:phases]), xlabel="Phase", ylabel="Air Flow Rate [kg/s]", dpi=800)
    plot!(p, x, results_expose[idx_key][:mdot_air_kgs],marker=markers[1])
    savefig(p, joinpath(save_dir, "mdotAir_Comp_$(filename)_$(range_expose).png"))
end
