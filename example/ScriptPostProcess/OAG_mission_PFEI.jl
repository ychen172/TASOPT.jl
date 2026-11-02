"""
This script compares the off-design PFEI of each OAG mission (4 bulk + 1 tail) across fuels,
for a few selected seat capacities.
    - One curve per (fuel case, seat capacity): PFEI vs OAG mission range
    - Models are read from merged "MI" folders (one seat layout per seat capacity)
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
save_name    = "OAG_Mission_PFEI.png"
# Cases: merged folder, file pattern (capture = seat capacity), and the fuel each case was optimized with
caseFolders  = ["Opti_Jet_NoACT_OAG_MI_TailCorr_V4_",
                "Opti_Eth_NoACT_OAG_MI_TailCorr_V4_"]
caseKeys     = [r"^Opti_Jet_NoACT_OAG_Ml_.*_TailCorr_V4_(\d+)\.jld2$",
                r"^Opti_Eth_NoACT_OAG_Ml_.*_TailCorr_V4_(\d+)\.jld2$"]
caseNames    = ["Jet Fuel", "Ethanol"]
idx_fuel_case      = [24      , 32       ] #[24      , 32       , 322431  ]
rho_fuel_case_kgm3 = [817.0   , 789.0    ] #[817.0   , 789.0    , 805.368 ]
hvap_fuel_case_Jkg = [358694.0, 918187.9 ] #[358694.0, 918187.9 , 586408.0]
# Seat capacities to compare (must exist in every case folder)
seat_caps    = [150, 200, 240]
# OAG mission set and payload fractions
miss_dir     = joinpath(model_dir,"OffDesignMissions_50_300_300_Tail.csv")
pass_load_frac_off  = 0.825 # Off-design payload load factor, matches opt_from_multi_warm_starts_para_oag.jl
pass_load_frac_tail = 0.850 # The farthest (tail) mission payload fraction
constraints  = [:WPay,:MWTO,:VolFuel]
iter_max     = 150
# Plotting
plt_linecolor = [:red  , :purple, :teal]   # Per case
plt_linestyle = [:solid, :dash  , :dot ]   # Per case
plt_marker    = [:diamond, :circle, :rect] # Per seat capacity

#### Create save directory
mkpath(save_dir)

#### Load OAG off-design mission data (ranges/weights per seat capacity)
miss_off_des = read_oag(miss_dir)

#### Run the 5 OAG missions for each (case, seat capacity)
Range_nmi = [[Float64[] for _ in seat_caps] for _ in caseFolders]
PFEI_JJ   = [[Float64[] for _ in seat_caps] for _ in caseFolders]
for (j, caseFolder) in enumerate(caseFolders)
    model_dir_sub = joinpath(model_dir, caseFolder)
    files = glob([caseKeys[j]], model_dir_sub)
    for (i, sc) in enumerate(seat_caps)
        # Find this seat capacity's model in the merged folder
        sc_files = filter(f -> parse(Int, match(caseKeys[j], basename(f))[1]) == sc, files)
        length(sc_files) == 1 || error("Expected exactly 1 model for $(sc) seats in $(model_dir_sub), found $(length(sc_files))")
        ac = quickload_aircraft(sc_files[1])
        println("File, $(sc_files[1]), read successfully")

        # OAG missions: 4 bulk at the off-design load factor, the farthest (tail) one at the tail load factor
        ranges_off_nmi = miss_off_des[sc].ranges_nmi
        wei_pay_off_N  = fill(ac.parg[igWpaymax]*pass_load_frac_off, length(ranges_off_nmi))
        wei_pay_off_N[argmax(ranges_off_nmi)] = ac.parg[igWpaymax]*pass_load_frac_tail
        out = off_design_specified!(ac, idx_fuel_case[j], rho_fuel_case_kgm3[j], hvap_fuel_case_Jkg[j],
                                    ranges_off_nmi, wei_pay_off_N;
                                    mod_ac_inplace=false, itermax=iter_max, constraints=constraints, save_model=false)

        # off_design_specified! silently drops missions that fail, so report any gap
        n_fail = length(ranges_off_nmi) - length(out.ran_nmi)
        n_fail > 0 && @warn "$(n_fail)/$(length(ranges_off_nmi)) OAG missions failed for $(caseNames[j]) at $(sc) seats"
        append!(Range_nmi[j][i], out.ran_nmi)
        append!(PFEI_JJ[j][i], out.PFEI_JJ)
        println("Case $(caseNames[j]) at $(sc) seats: PFEI = $(round.(out.PFEI_JJ, digits=4)) J/J")
    end
end

#### Plotting: PFEI vs OAG mission range, one curve per (case, seat capacity)
p = plot(xlabel="OAG Mission Range [nmi]", ylabel="Off-Design PFEI [J/J]", dpi=800, legend=:best)
for (j, caseName) in enumerate(caseNames), (i, sc) in enumerate(seat_caps)
    color = plt_linecolor[j]
    plot!(p, Range_nmi[j][i], PFEI_JJ[j][i],
          color=color, linestyle=plt_linestyle[j], lw=2,
          marker=plt_marker[i], mc=color, msc=color,
          label="$(caseName), $(sc) Seats")
end
savefig(p, joinpath(save_dir, save_name))
println("Saved PFEI vs range plot to $(joinpath(save_dir, save_name))")
