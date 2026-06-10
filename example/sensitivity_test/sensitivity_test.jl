"""This script test the sensitivity for a group of parameters around optimized solution"""

using TASOPT
include(TASOPT.__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils","sensitivity.jl"))

modelPath = joinpath(@__DIR__,"../ModelSaved/acOptim_BatJet_CT/acOptim_BatJet_CT3000.jld2")
println(modelPath)
ac = quickload_aircraft(modelPath)
ac.htail.opt_sizing = TailSizing.CLmaxFwdCG
ac.htail.CL_max_fwd_CG = -0.7
ac.vtail.opt_sizing = TailSizing.OEI
ac.parg[igCLveout] = 0.5

input_params = [
    :(ac.wing.layout.ηs),                                     # Panel break eta location
    :(ac.parg[igxeng]),                                       # Engine axial location
    :(ac.parg[igyeng]),                                       # Engine vertical location
    :(ac.wing.layout.AR),                                     # Wing aspect ratio
    :(ac.para[iaCL, ipcruise1, 1]),                           # Cruise CL
    :(ac.wing.layout.sweep),                                  # Wing sweep
    :(ac.para[iaalt, ipcruise1, 1]),                          # Cruise altitude
    :(ac.wing.inboard.λ),                                     # Inboard taper ratio
    :(ac.wing.outboard.λ),                                    # Outboard taper ratio
    :(ac.wing.inboard.cross_section.thickness_to_chord),      # Inboard t/c
    :(ac.wing.outboard.cross_section.thickness_to_chord),     # Outboard t/c
    :(ac.para[iarcls, ipcruise1, 1]),                         # Break/root Cl ratio
    :(ac.para[iarclt, ipcruise1, 1]),                         # Tip/root Cl ratio
    :(ac.pare[ieTt4, ipcruise1, 1]),                          # Tt4 at cruise
    :(ac.pare[iepihc, ipcruise1, 1]),                         # HPC pressure ratio
    :(ac.pare[iepif, ipcruise1, 1]),                          # Fan pressure ratio
    :(ac.pare[iepilc, ipcruise1, 1]),                         # LPC pressure ratio
    :(ac.pare[ieBPR, ipcruise1, 1]),                          # Bypass ratio
    :(ac.htail.CL_max_fwd_CG),                                # CL max horizontal tail (probably skip)
    :(ac.parg[igCLveout])                                     # CL max vertical tail (Skip due to duplicated effect with yEng and also bounded by physical limit)
]

names_params = [
    "ηs",
    "xEng",
    "yEng",
    "AR",
    "CL_cruise",
    "sweep",
    "alt_cruise",
    "λ_inboard",
    "λ_outboard",
    "t/c_inboard",
    "t/c_outboard",
    "Cl_break_root",
    "Cl_tip_root",
    "Tt4_cruise",
    "pi_hc",
    "pi_f",
    "pi_lc",
    "BPR",
    "CL_h",
    "CL_v"
]
sensitivityVector = get_sensitivity(input_params; model_state=ac, eps=1e-5, optimizer=false, f_out_fn=nothing, diff_scheme=:central, metric=:impact)
println(sensitivityVector)