"""
Compare aircraft and engine performance between a model and a reference
engine parameters.
"""

using TASOPT
using Plots
using CSV, DataFrames
using Printf
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_engines/run_engine.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_offdesign/offdesign.jl"))
using .PRD: findR1R2R3

## IO Directory
miss_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/reference_engine_cycle_CFM5B/Target_Leap1B.csv")
save_key = "Leap1B_Max8_210S_M078"
save_dir = joinpath(__TASOPTroot__,"../example/ModelProcessed",save_key)
mkpath(save_dir)

## Load EEDB reference mission (engine comparison)
df_eedb       = CSV.read(miss_dir, DataFrame)
Fn_N          = Float64.(df_eedb[:, "Thrust (kN)"] .* 1000.0) #[N]
BPR_Ref       = df_eedb[:, "BPR"]
OPR_Ref       = df_eedb[:, "OPR"]
Wfuel_Ref_kgs = Float64.(df_eedb[:, "Wf[lbm/s]"] .* 0.453592)

## Setup engine inlet conditions (EEDB reference condition)
M0 = 0.0
P0 = 101320.0 #Pa
T0 = 288.2 #K
a0 = 340.2074661144284 #m/s
Fan_Dia_Ref_m = 1.76 #m

## Offdesign fuel parameters (Jet, matching the reference aircraft)
fuel_idx  = 24
rho_fuel  = 817.0    #kg/m3
hvap_fuel = 358694.0 #J/kg

## R1/R2/R3 search bounds -- found directly via findR1R2R3's own internal bisection
## (not the full off_design_PRD envelope sweep), so no need to pre-derive tighter bounds.
R_LB = 300.0  #nmi
R_UB = 5000.0 #nmi
epsRange = 1e-5
epsWpay  = 1e-8

## Candidates to compare -- each is a (label, path) pair for an already-sized/optimized aircraft
## model. Loaded directly with quickload_aircraft, no resizing.
model_cases = [
     ("Partial Optimization with Fixed PR-ratio",
     joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_CalibReproduced_Max8_/Opti_Jet_NoACT_V4_CalibReproduced_Max8_2621.jld2")),
     ("Full Global Optimization with Fixed LPC PR",
     joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Glo_FixLPC_/Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Glo_FixLPC_2621.jld2")),
     ("Full Local Optimization with Fixed LPC PR",
     joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Loc_FixLPC_/Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Loc_FixLPC_2621.jld2")),
     ("Full Global Optimization with Free PR-ratio",
     joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Glo_/Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Glo_2621.jld2")),
     ("Full Local Optimization with Free PR-ratio",
     joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Loc_/Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Loc_2621.jld2")),
]

## Reference 737MAX10 numbers
# ref_name = "B737Max10"
# ref_air  = Dict(
#     "MWTO_Ton"         => 89.8,
#     "WZF_R1_Ton"       => 72.6,
#     "R1_nmi"           => 2419.0,
#     "WZF_R2_Ton"       => 69.0,
#     "R2_nmi"           => 3094.0,
#     "WZF_R3_Ton"       => 48.6,
#     "R3_nmi"           => 4120.0,
#     "L_fuse_apu_end_m" => 43.0,
#     "R_fuse_m"         => 1.88,
#     "b_wing_m"         => 35.92,
#     "AR_wing"          => 10.16,
#     "Sweep_wing_deg"   => 25.03,
#     "S_vtail_m2"       => 31.62,
#     "AR_vtail"         => 1.94,
#     "b_htail_m"        => 14.35,
#     "x_box_wing_m"     => 20.5,
#     "x_box_vtail_m"    => 39.0,
#     "x_box_htail_m"    => 40.5,
#     "y_engine_m"       => 4.82,
#     "F_eng_static_kN" => 130.4,
#     "BPR_eng_static"  => 8.3,
#     "OPR_eng_static"  => 42.0,
#     "Dia_eng_m"        => 1.76,
# )

ref_name = "B737Max8"
ref_air  = Dict(
    "MWTO_Ton"         => 82.6,
    "WZF_R1_Ton"       => 65.952,
    "R1_nmi"           => 2621.4,
    "WZF_R2_Ton"       => 61.99,
    "R2_nmi"           => 3480.7,
    "WZF_R3_Ton"       => 44.67,
    "R3_nmi"           => 4427.5,
    "L_fuse_apu_end_m" => 39.47,
    "R_fuse_m"         => 1.88,
    "b_wing_m"         => 35.92,
    "AR_wing"          => 10.16,
    "Sweep_wing_deg"   => 25.03,
    "S_vtail_m2"       => 31.62,
    "AR_vtail"         => 1.94,
    "b_htail_m"        => 14.35,
    "x_box_wing_m"     => 18.09,
    "x_box_vtail_m"    => 33.45,
    "x_box_htail_m"    => 36.83,
    "y_engine_m"       => 4.82,
    "F_eng_static_kN" => 130.4,
    "BPR_eng_static"  => 8.3,
    "OPR_eng_static"  => 42.0,
    "Dia_eng_m"        => 1.76,
)

## Engine off-design results, one column per candidate (same layout/order as model_cases)
numTry = length(model_cases)
numFn  = length(Fn_N)
Fn_N_Calc      = fill(NaN, numFn, numTry)
Wfuel_Calc_kgs = fill(NaN, numFn, numTry)
BPR_Calc       = fill(NaN, numFn, numTry)
OPR_Calc       = fill(NaN, numFn, numTry)
DFan_Calc_m    = fill(NaN, numFn, numTry)

## Collect one aircraft-comparison row per candidate, and fill in its EEDB off-design sweep
rows = Vector{Dict{String,Any}}()
for (idx, (label, path)) in enumerate(model_cases)
    local ac_local
    try
        ac_local = quickload_aircraft(path)
    catch e
        @warn "Case \"$(label)\" failed to load from $(path), skipping" exception=e
        continue
    end
    @assert ac_local.is_sized[1] "Case \"$(label)\" aircraft model is not marked as sized"

    row = Dict{String,Any}("Case" => label)

    ## Design-point geometry/engine data
    row["MWTO_Ton"]       = ac_local.parg[igWMTO]/gee/1000.0
    row["L_fuse_apu_end_m"] = ac_local.fuselage.layout.x_cone_end
    row["R_fuse_m"]       = ac_local.fuselage.layout.radius
    row["b_wing_m"]       = ac_local.wing.layout.span
    row["AR_wing"]        = ac_local.wing.layout.AR
    row["Sweep_wing_deg"] = ac_local.wing.layout.sweep
    row["S_vtail_m2"]     = ac_local.vtail.layout.S
    row["AR_vtail"]       = ac_local.vtail.layout.AR
    row["b_htail_m"]      = ac_local.htail.layout.span
    row["x_box_wing_m"]   = ac_local.wing.layout.box_x
    row["x_box_vtail_m"]  = ac_local.vtail.layout.box_x
    row["x_box_htail_m"]  = ac_local.htail.layout.box_x
    row["y_engine_m"]     = ac_local.parg[igyeng]
    row["F_eng_static_kN"] = ac_local.pare[ieFe,ipstatic,1]/1000.0
    row["Dia_eng_m"]      = ac_local.parg[igdfan]
    row["BPR_eng_static"] = ac_local.pare[ieBPR,ipstatic,1]
    row["OPR_eng_static"] = ac_local.pare[ieOPR,ipstatic,1]
    row["pifan_cru"]      = ac_local.pare[iepif,ipcruise1,1]
    row["pilpc_cru"]      = ac_local.pare[iepilc,ipcruise1,1]
    row["pihpc_cru"]      = ac_local.pare[iepihc,ipcruise1,1]
    row["BPR_cru"]        = ac_local.pare[ieBPR,ipcruise1,1]
    row["Tt4_cru"]        = ac_local.pare[ieTt4,ipcruise1,1]
    row["OPR_cru"]        = ac_local.pare[iepilc,ipcruise1,1]*ac_local.pare[iepihc,ipcruise1,1]
    ## Technology level
    row["pib_cru"]        = ac_local.pare[iepib,ipcruise1,1]
    row["etapfan_cru"]    = ac_local.pare[ieepolf,ipcruise1,1]
    row["etaplpc_cru"]    = ac_local.pare[ieepollc,ipcruise1,1]
    row["etaphpc_cru"]    = ac_local.pare[ieepolhc,ipcruise1,1]
    row["etaphpt_cru"]    = ac_local.pare[ieepolht,ipcruise1,1]
    row["etaplpt_cru"]    = ac_local.pare[ieepollt,ipcruise1,1]
    row["etab_cru"]       = ac_local.pare[ieetab,ipcruise1,1]
    row["TmetalMax"]      = ac_local.parg[igTmetal]


    ## Design-point OEW (mission-invariant), used to convert R1/R2/R3 payload back to zero-fuel weight
    OEW_N = ac_local.parm[imWTO,1] - ac_local.parm[imWfuel,1] - ac_local.parm[imWpay,1]

    ## Find R1/R2/R3 directly via bisection (no full envelope sweep, no saved .jld2 output)
    for mode in (:R1, :R2, :R3)
        try
            out = findR1R2R3(mode, R_LB, R_UB, ac_local, fuel_idx, rho_fuel, hvap_fuel;
                             epsRange=epsRange, epsWpay=epsWpay, flg_save_ac=false)
            row["$(mode)_nmi"]    = out["range_nmi"][1]
            row["WZF_$(mode)_Ton"] = (OEW_N + out["payload_weight_N"][1])/gee/1000.0
        catch e
            @warn "Case \"$(label)\" failed to find $(mode)" exception=e
            row["$(mode)_nmi"]    = missing
            row["WZF_$(mode)_Ton"] = missing
        end
    end

    push!(rows, row)

    ## Run the engine through the EEDB off-design thrust sweep -- zero_offtake=true to match EEDB's
    ## no-customer-offtake measurement convention
    for (idx_Fn,Fn_N_cur) in enumerate(Fn_N)
        try
            res = RunEngine.runOffDes(ac_local, M0, P0, T0, a0, Fn_N_cur; zero_offtake=true)
            Fn_N_Calc[idx_Fn, idx]      = res.Fe        #[N]
            Wfuel_Calc_kgs[idx_Fn, idx] = res.mcore * res.ff #[kg/s]
            BPR_Calc[idx_Fn, idx]       = res.BPR
            OPR_Calc[idx_Fn, idx]       = res.OPR
            DFan_Calc_m[idx_Fn, idx]    = ac_local.parg[igdfan] #[m]
        catch e
            @warn "Case \"$(label)\" failed at Fn_N=$(Fn_N_cur), leaving as NaN" exception=e
        end
    end
end

########################################
## Part 1: Aircraft comparison 
########################################

## Reference row, prepended alongside the model candidates so it appears as its own
## column once the table is transposed below.
ref_row = Dict{String,Any}("Case" => "Reference ($(ref_name))")
for (k,v) in ref_air
    ref_row[k] = v
end
rows_with_ref = vcat([ref_row], rows)

## Assemble and display comparison table
all_cols = unique(vcat([collect(keys(r)) for r in rows_with_ref]...))
col_order = ["Case","MWTO_Ton","R1_nmi","WZF_R1_Ton","R2_nmi","WZF_R2_Ton","R3_nmi","WZF_R3_Ton",
             "b_wing_m","AR_wing","Sweep_wing_deg","S_vtail_m2","AR_vtail","b_htail_m",
             "L_fuse_apu_end_m","R_fuse_m","x_box_wing_m","x_box_vtail_m","x_box_htail_m",
             "y_engine_m","F_eng_static_kN","BPR_eng_static","OPR_eng_static",
             "pifan_cru","pilpc_cru","pihpc_cru","OPR_cru","BPR_cru","Tt4_cru",
             "pib_cru","etapfan_cru","etaplpc_cru","etaphpc_cru","etaphpt_cru","etaplpt_cru","etab_cru","TmetalMax",
             "Dia_eng_m"]
col_order = filter(c -> c in all_cols, col_order)

df = DataFrame([c => [get(r,c,missing) for r in rows_with_ref] for c in col_order])
println()
println("=== Aircraft comparison ($(nrow(df)) candidates incl. reference) ===")
show(df, allcols=true)
println()

## Transpose so each parameter is a row and each candidate (plus the reference) is its own column --
## easier to scan one parameter across all candidates at a glance.
param_names = filter(c -> c != "Case", col_order)
df_t = DataFrame(Parameter = param_names)
for i in 1:nrow(df)
    df_t[!, Symbol(df.Case[i])] = [df[i, p] for p in param_names]
end
println()
println("=== Aircraft comparison, transposed ===")
show(df_t, allcols=true)
println()

CSV.write(joinpath(save_dir,"Compare_Aircraft_Models.csv"), df_t)
println()
println("Saved aircraft comparison table to: ", joinpath(save_dir,"Compare_Aircraft_Models.csv"))

## Reference vs Model vs Difference%, one table per candidate.
## Uses BPR_eng_static/OPR_eng_static to match the reference row concept directly -- caveat above applies.
# diff_param_order = ["MWTO_Ton","WZF_R1_Ton","R1_nmi","WZF_R2_Ton","R2_nmi","WZF_R3_Ton","R3_nmi",
#                      "L_fuse_apu_end_m","R_fuse_m","b_wing_m","AR_wing","Sweep_wing_deg",
#                      "S_vtail_m2","AR_vtail","b_htail_m","x_box_wing_m","x_box_vtail_m","x_box_htail_m",
#                      "y_engine_m","F_eng_static_kN","BPR_eng_static","OPR_eng_static","Dia_eng_m"]
# for row in rows
#     println()
#     println("=== $(row["Case"]) vs $(ref_name) reference ===")
#     @printf("%-20s %12s %12s %12s\n","Parameter","Reference","Model","Difference")
#     for p in diff_param_order
#         ref_val = get(ref_air, p, missing)
#         mod_val = get(row, p, missing)
#         if ismissing(ref_val) || ismissing(mod_val)
#             @printf("%-20s %12s %12s %12s\n", p, string(ref_val), string(mod_val), "n/a")
#         else
#             diff_pct = 100*abs(mod_val-ref_val)/abs(ref_val)
#             @printf("%-20s %12.4g %12.4g %11.3f%%\n", p, ref_val, mod_val, diff_pct)
#         end
#     end
# end

########################################
## Part 2: Engine comparison (EEDB)
########################################

case_legend = [label for (label,path) in model_cases]
const MARKERS = [:rect, :circle, :diamond, :utriangle, :dtriangle, :cross, :xcross, :star4, :star5, :star6, :star7, :star8]

# Fuel flow rate
p = plot(xlabel="Thrust [kN]", ylabel="Fuel Flow [kg/s]", dpi=800, legend=:best)
for idx in 1:numTry
    plot!(p, Fn_N_Calc[:,idx] ./ 1000.0,  Wfuel_Calc_kgs[:,idx], marker=MARKERS[idx+1], label=case_legend[idx])
end
plot!(p, Fn_N ./ 1000.0,  Wfuel_Ref_kgs, marker=MARKERS[1], label="Reference")
savefig(p,joinpath(save_dir,"$(save_key)_Fuel_Flow.png"))

# Bypass ratio
p = plot(xlabel="Thrust [kN]", ylabel="Bypass Ratio", dpi=800, legend=:best)
for idx in 1:numTry
    plot!(p, Fn_N_Calc[:,idx] ./ 1000.0,  BPR_Calc[:,idx], marker=MARKERS[idx+1], label=case_legend[idx])
end
plot!(p, Fn_N ./ 1000.0,  BPR_Ref, marker=MARKERS[1], label="Reference")
savefig(p,joinpath(save_dir,"$(save_key)_Bypass_Ratio.png"))

# Overall Pressure Ratio
p = plot(xlabel="Thrust [kN]", ylabel="Overall Pressure Ratio", dpi=800, legend=:best)
for idx in 1:numTry
    plot!(p, Fn_N_Calc[:,idx] ./ 1000.0,  OPR_Calc[:,idx], marker=MARKERS[idx+1], label=case_legend[idx])
end
plot!(p, Fn_N ./ 1000.0,  OPR_Ref, marker=MARKERS[1], label="Reference")
savefig(p,joinpath(save_dir,"$(save_key)_Overall_Pressure_Ratio.png"))

# Fan Diameter
p = plot(xlabel="Thrust [kN]", ylabel="Fan Diameter [m]", dpi=800, legend=:best)
for idx in 1:numTry
    plot!(p, Fn_N_Calc[:,idx] ./ 1000.0,  DFan_Calc_m[:,idx], marker=MARKERS[idx+1], label=case_legend[idx])
end
plot!(p, Fn_N ./ 1000.0,  ones(length(Fn_N)) .* Fan_Dia_Ref_m, marker=MARKERS[1], label="Reference")
savefig(p,joinpath(save_dir,"$(save_key)_Fan_Diameter.png"))

println()
println("Saved engine comparison plots to: ", save_dir)
println("DONE")
