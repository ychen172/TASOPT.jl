"""
This script runs the OAG missions (4 bulk + 1 tail) for every seat capacity in each fuel case and compares:
    - Percentage per-mission PFEI penalty of the 2nd case (Ethanol) relative to the 1st (Jet), per selected seat capacity
    - Start-of-cruise L/D of the tail mission vs seat capacity, one curve per fuel case
    - Percentage OAG-weighted PFEI penalty of the 2nd case (Ethanol) relative to the 1st (Jet)
    - R1 payload vs R1 range across all seat capacities, one curve per fuel case (R1 = design mission, sized at max payload)
Models are read from merged "MI" folders (one seat layout per seat capacity)
"""

using TASOPT
include(__TASOPTindices__)
using Plots
using Glob
include(joinpath(@__DIR__,"../UtilityFunctions/objective_factory.jl"))
off_design_specified! = ObjectiveFactory.off_design_specified!
include(joinpath(@__DIR__,"../UtilityFunctions/Extract.jl"))
using .Extract: read_oag

#### Setup IO
model_dir    = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir     = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
save_name_mission_penalty = "OAG_Mission_PFEI_EthPenalty.png"
save_name_LD       = "OAG_Mission_LD.png"
save_name_penalty  = "OAG_Weighted_PFEI_EthPenalty.png"
save_name_R1       = "OAG_PRD_R1.png"
# Cases: merged folder, file pattern (capture = seat capacity), and the fuel each case was optimized with
caseFolders  = ["Opti_Jet_NoACT_OAG_MI_TailCorr_V4_",
                "Opti_Eth_NoACT_OAG_MI_TailCorr_V4_"]
caseKeys     = [r"^Opti_Jet_NoACT_OAG_Ml_.*_TailCorr_V4_(\d+)\.jld2$",
                r"^Opti_Eth_NoACT_OAG_Ml_.*_TailCorr_V4_(\d+)\.jld2$"]
caseNames    = ["Jet Fuel", "Ethanol"]
idx_fuel_case      = [24      , 32       ] #[24      , 32       , 322431  ]
rho_fuel_case_kgm3 = [817.0   , 789.0    ] #[817.0   , 789.0    , 805.368 ]
hvap_fuel_case_Jkg = [358694.0, 918187.9 ] #[358694.0, 918187.9 , 586408.0]
# Seat capacities shown in the per-mission plot (must exist in every case folder); every seat capacity is computed
seat_caps    = [140, 180, 220]
# OAG mission set and payload fractions
miss_dir     = joinpath(model_dir,"OffDesignMissions_50_300_300_Tail.csv")
pass_load_frac_off  = 0.825 # Off-design payload load factor, matches opt_from_multi_warm_starts_para_oag.jl
pass_load_frac_tail = 0.850 # The farthest (tail) mission payload fraction
constraints  = [:WPay,:MWTO,:VolFuel]
iter_max     = 150
# Plotting
plt_linecolor = [:red  , :purple, :teal, :blue, :orange, :green, :brown, :magenta, :gray40, :olive]                 # Per case (per seat in the penalty plot)
plt_linestyle = [:solid, :dash  , :dot , :dashdot, :dashdotdot, :solid, :dash, :dot, :dashdot, :dashdotdot]          # Per case
plt_marker    = [:diamond, :circle, :rect, :utriangle, :dtriangle, :star5, :hexagon, :pentagon, :xcross, :cross] # Per seat capacity

#### Create save directory
mkpath(save_dir)

#### Load OAG off-design mission data (ranges/weights per seat capacity)
miss_off_des = read_oag(miss_dir)

#### Run the 5 OAG missions for every seat capacity in each case
# results[j][sc] = (ran_nmi, PFEI_JJ, LD_tail, weighted_PFEI, R1_ran_nmi, R1_wpay_N) for case j at seat capacity sc
results = [Dict{Int,NamedTuple}() for _ in caseFolders]
for (j, caseFolder) in enumerate(caseFolders)
    model_dir_sub = joinpath(model_dir, caseFolder)
    for f in sort(glob([caseKeys[j]], model_dir_sub), by = f -> parse(Int, match(caseKeys[j], basename(f))[1]))
        sc = parse(Int, match(caseKeys[j], basename(f))[1])
        haskey(miss_off_des, sc) || continue #No OAG mission data for this seat capacity
        ac = quickload_aircraft(f)
        println("File, $(f), read successfully")

        # OAG missions: 4 bulk at the off-design load factor, the farthest (tail) one at the tail load factor
        ranges_off_nmi = miss_off_des[sc].ranges_nmi
        weights_off    = miss_off_des[sc].weights
        wei_pay_off_N  = fill(ac.parg[igWpaymax]*pass_load_frac_off, length(ranges_off_nmi))
        wei_pay_off_N[argmax(ranges_off_nmi)] = ac.parg[igWpaymax]*pass_load_frac_tail
        # Fly each mission on its own copy (mod_ac_inplace=true) so the flown mission-2 state is available for L/D
        ran_ok, pfei_ok = Float64[], Float64[]
        LD_tail = NaN # Start-of-cruise L/D of the tail mission (NaN if it failed)
        num, den = 0.0, 0.0 # OAG-weighted PFEI = sum(w*PFEI*payload*range) / sum(w*payload*range), over converged missions
        for k in eachindex(ranges_off_nmi)
            ac_k = deepcopy(ac)
            out = off_design_specified!(ac_k, idx_fuel_case[j], rho_fuel_case_kgm3[j], hvap_fuel_case_Jkg[j],
                                        [ranges_off_nmi[k]], [wei_pay_off_N[k]];
                                        mod_ac_inplace=true, itermax=iter_max, constraints=constraints, save_model=false)
            if isempty(out.ran_nmi) # off_design_specified! drops missions that fail
                @warn "OAG mission $(round(ranges_off_nmi[k], digits=1)) nmi failed for $(caseNames[j]) at $(sc) seats"
                continue
            end
            push!(ran_ok, out.ran_nmi[1])
            push!(pfei_ok, out.PFEI_JJ[1])
            if k == argmax(ranges_off_nmi)
                LD_tail = ac_k.para[iaCL,ipcruise1,2] / ac_k.para[iaCD,ipcruise1,2]
            end
            num += weights_off[k] * out.PFEI_JJ[1] * out.wei_pay_N[1] * out.ran_nmi[1]
            den += weights_off[k] * out.wei_pay_N[1] * out.ran_nmi[1]
        end
        results[j][sc] = (ran_nmi = ran_ok, PFEI_JJ = pfei_ok, LD_tail = LD_tail,
                          weighted_PFEI = den > 0.0 ? num/den : NaN,
                          R1_ran_nmi = ac.parm[imRange,1]/1852.0, R1_wpay_N = ac.parm[imWpay,1]) # Design mission = R1 (sized at max payload)
        println("Case $(caseNames[j]) at $(sc) seats: OAG-weighted PFEI = $(round(results[j][sc].weighted_PFEI, digits=4)) J/J")
    end
end

#### Plotting: per-mission PFEI penalty of case 2 (Ethanol) relative to case 1 (Jet), one curve per selected seat capacity
Range_nmi = [[results[j][sc].ran_nmi for sc in seat_caps] for j in eachindex(caseFolders)]
PFEI_JJ   = [[results[j][sc].PFEI_JJ for sc in seat_caps] for j in eachindex(caseFolders)]
p_mp = plot(xlabel="OAG Mission Range [nmi]", ylabel="$(caseNames[2]) PFEI Penalty [%]", dpi=800, legend=:best)
for (i, sc) in enumerate(seat_caps)
    # The penalty pairs missions one-to-one, so both cases must have kept the same missions
    (length(Range_nmi[1][i]) == length(Range_nmi[2][i]) && Range_nmi[1][i] ≈ Range_nmi[2][i]) ||
        error("$(caseNames[1]) and $(caseNames[2]) converged different OAG missions at $(sc) seats; cannot pair them")
    color = plt_linecolor[i]
    plot!(p_mp, Range_nmi[1][i], 100.0 .* (PFEI_JJ[2][i] .- PFEI_JJ[1][i]) ./ PFEI_JJ[1][i],
          color=color, lw=2, marker=plt_marker[i], mc=color, msc=color, label="$(sc) Seats")
end
savefig(p_mp, joinpath(save_dir, save_name_mission_penalty))
println("Saved per-mission PFEI penalty plot to $(joinpath(save_dir, save_name_mission_penalty))")

#### Plotting: tail-mission start-of-cruise L/D vs seat capacity, one curve per case
p_ld = plot(xlabel="Seat Capacity", ylabel="L/D Ratio at Cruise for the Tail Mission", dpi=800, legend=:best)
for (j, caseName) in enumerate(caseNames)
    scs = sort(collect(keys(results[j])))
    plot!(p_ld, scs, [results[j][sc].LD_tail for sc in scs],
          color=plt_linecolor[j], linestyle=plt_linestyle[j], lw=2,
          marker=:circle, mc=plt_linecolor[j], msc=plt_linecolor[j], label=caseName)
end
savefig(p_ld, joinpath(save_dir, save_name_LD))
println("Saved tail-mission L/D plot to $(joinpath(save_dir, save_name_LD))")

#### Plotting: percentage OAG-weighted PFEI penalty of case 2 (Ethanol) relative to case 1 (Jet), common seat capacities
scs_common = sort(collect(intersect(keys(results[1]), keys(results[2]))))
penalty_pct = [100.0 * (results[2][sc].weighted_PFEI - results[1][sc].weighted_PFEI) / results[1][sc].weighted_PFEI for sc in scs_common]
p_pen = plot(scs_common, penalty_pct, xlabel="Seat Capacity", ylabel="$(caseNames[2]) Mission-averaged PFEI Penalty [%]",
             dpi=800, legend=false, color=plt_linecolor[2], lw=2, marker=:circle, mc=plt_linecolor[2], msc=plt_linecolor[2])
hline!(p_pen, [0.0], color=:black, lw=1)
savefig(p_pen, joinpath(save_dir, save_name_penalty))
println("Saved OAG-weighted PFEI penalty plot to $(joinpath(save_dir, save_name_penalty))")
for (sc, pen) in zip(scs_common, penalty_pct)
    println("  $(sc) seats: $(caseNames[2]) penalty = $(round(pen, digits=2)) %")
end

#### Plotting: R1 payload vs R1 range across all seat capacities, one curve per case
p_R1 = plot(xlabel="R1 Range [nmi]", ylabel="R1 Payload Weight [Ton]", dpi=800, legend=:best)
for (j, caseName) in enumerate(caseNames)
    scs = sort(collect(keys(results[j])))
    plot!(p_R1, [results[j][sc].R1_ran_nmi for sc in scs], [results[j][sc].R1_wpay_N for sc in scs] ./ 9.81 ./ 1000.0,
          color=plt_linecolor[j], linestyle=plt_linestyle[j], lw=2,
          marker=:circle, mc=plt_linecolor[j], msc=plt_linecolor[j], label=caseName)
end
savefig(p_R1, joinpath(save_dir, save_name_R1))
println("Saved R1 payload vs range plot to $(joinpath(save_dir, save_name_R1))")
