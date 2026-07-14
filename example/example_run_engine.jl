"""
This script intends to probe the sea-level static operation performance of engine model optimized for different aircraft design
"""

using TASOPT
include(__TASOPTindices__)
using CSV, DataFrames
include(joinpath(__TASOPTroot__,"../example/utilities_engines/run_engine.jl"))
using .RunEngine:runOffDes

#### IO
save_name = "Engine_Opti_Jet_NoACT_2900" #Save for engine off-design performance
save_dir = joinpath(__TASOPTroot__,"../example/ModelSaved")
ac_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_/Opti_Jet_NoACT_2900.jld2")
Fn_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/reference_engine_cycle_CFM5B/CFM5B_Finepts_NewLHV.csv")

#### Create a save directory
save_dir = joinpath(save_dir,save_name)
mkpath(save_dir)

#### Load aircraft model that has the sized engine model
ac = quickload_aircraft(ac_dir)

#### Load the requested thrust levels
df = CSV.read(Fn_dir, DataFrame)
Fn_run = df[:, "Thrust (kN)"] .* 1000.0 #[N]

#### Assume using the TASOPT sea-level static inlet conditions
M0 = 0.0
p0 = 101320.0 #Pa
T0 = 288.2 #K
a0 = 340.2074661144284 #m/s

#### Run each of the off-design thrust requirements
num_Fn = length(Fn_run)
Fn_kN = Fn_run./1000.0 #[kN]
flg_conv = Vector{Bool}(undef,num_Fn)
Pt3_psi = zeros(num_Fn)
Tt3_R = zeros(num_Fn)
Wf_lbms = zeros(num_Fn)
BPR = zeros(num_Fn)
W3_lbms = zeros(num_Fn)
T04_R = zeros(num_Fn)
OPR = zeros(num_Fn)
for (idx,Fn_cur) in enumerate(Fn_run)
    res_cur = runOffDes(ac,M0,p0,T0,a0,Fn_cur)
    println("For thrust $(Fn_cur/1000.0) kN, solution converged: $(res_cur.Lconv)")
    flg_conv[idx] = res_cur.Lconv
    Pt3_psi[idx] = res_cur.pt3 * 0.000145038 #[psi]
    Tt3_R[idx] = res_cur.Tt3 * 9.0/5.0 #[R]
    Wf_lbms[idx] = res_cur.mcore * res_cur.ff * 2.20462 #[lbm/s]
    BPR[idx]  = res_cur.Fe / res_cur.Fsp / res_cur.u0 / res_cur.mcore - 1.0
    W3_lbms[idx] = 2.20462 * (res_cur.Fe + res_cur.u0 * res_cur.mcore - BPR[idx] * (res_cur.u8 - res_cur.u0) * res_cur.mcore)/res_cur.u6 #[lbm/s]
    T04_R[idx] = res_cur.Tt4 * 9.0/5.0 #[R]
    OPR[idx] = res_cur.pilc * res_cur.pihc
end
W4_lbms = W3_lbms .+ Wf_lbms

#### Save data
df = DataFrame(
    "Converged"   => flg_conv,
    "Thrust (kN)" => Fn_kN,
    "Pt3[psi]"    => Pt3_psi,
    "Tt3[R]"      => Tt3_R,
    "W4[lbm/s]"   => W4_lbms,
    "W3[lbm/s]"   => W3_lbms,
    "Wf[lbm/s]"   => Wf_lbms,
    "T04[R]"      => T04_R,
    "BPR"         => BPR,
    "OPR"         => OPR,
)
CSV.write(joinpath(save_dir, save_name*".csv"), df)