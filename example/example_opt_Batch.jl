using TASOPT, NLopt
include(joinpath(@__DIR__, "objective_factory.jl"))
using .ObjectiveFactory: OptHistory, make_obj
include(__TASOPTindices__)

#### IO Prepare
# Make a folder for saving optimized aircraft model
save_dir = "ModelSaved"
mkpath(save_dir)

# Make a save name for the optimized model
save_Name = "acOptimized_Cus2"

#### Loading an Baseline Aircraft Model
ac = read_aircraft_model("./customized/narrow_input.toml"; templatefile = "./customized/narrow_input.toml")
size_aircraft!(ac)

#### Modify the mission requirement for optimization (ToDo)

#### Get the objective function for optimization
hist_optim = OptHistory() #Optimization history
# Setup constraints
max_span = 35.814 #Maximum Span(m)
max_lenField = 2438.4 #Maximum Balanced Field Length(m)
min_TOCGamma = 0.015 #Minimum Top of Climb Flight Angle (rad)
max_Tt3 = 900.0 #Maximum Combustor Inlet Temperature (K)
max_TMetal = 1333.33 #Maximum Metal Tempeature (K)
max_DiaFan = 2.0 #Maximum Fan Diameter (m)
constraints_optim = [max_span, max_lenField, min_TOCGamma, max_Tt3, max_TMetal, max_DiaFan]
# Get objective function
obj = make_obj(ac, constraints_optim, hist_optim)

#### Setup the optimized parameters
""" These are the parameters to optimize
       x: [AR,    CL,  sweep(deg), altitude, λ_in,  λ_out, t/c_root, t/c_span, rcls, rclt,   Tt4,  π_hc,   π_f,   π_lc,   BPR]
            1.    2.    3.          4.         5.     6.     7.         8.       9.   10.     11.    12.    13.    14.     15.
""" 
lower   = [6.0,  0.45, 25.0,      10000.0,    0.65,  0.1,   0.125,    0.125,    0.9,  0.7,  1400.0, 10.0,  1.25,   2.99,   1.0]
upper   = [18.0, 0.75, 30.0,      20000.0,    0.85,  0.4,   0.15,     0.15,     1.3,  1.0,  1650.0, 15.0,  2.0,    3.01,  20.0]
init_dx = [0.5,  0.05,  0.1,       200.0,     0.01,  0.01,  0.01,     0.01,     0.01, 0.01, 100.0,   0.5,  0.05,   0.001  ,1.0]
initial = [
    ac.wing.layout.AR,                                     # 1.  Current aspect ratio
    ac.para[iaCL, ipcruise1, 1],                           # 2.  Current cruise CL
    ac.wing.layout.sweep,                                  # 3.  Current sweep angle (deg)
    ac.para[iaalt, ipcruise1, 1],                          # 4.  Current cruise altitude
    ac.wing.inboard.λ,                                     # 5.  Current inner taper ratio
    ac.wing.outboard.λ,                                    # 6.  Current outer taper ratio
    ac.wing.inboard.cross_section.thickness_to_chord,      # 7.  Current root t/c
    ac.wing.outboard.cross_section.thickness_to_chord,     # 8.  Current span t/c
    ac.para[iarcls, ipcruise1, 1],                         # 9.  Current rcls (Cl ratios)
    ac.para[iarclt, ipcruise1, 1],                         # 10. Current rclt
    ac.pare[ieTt4, ipcruise1, 1],                          # 11. Current Tt4
    ac.pare[iepihc, ipcruise1, 1],                         # 12. Current π_hc
    ac.pare[iepif, ipcruise1, 1],                          # 13. Current π_f
    ac.pare[iepilc, ipcruise1, 1],                         # 14. Current π_lc
    ac.pare[ieBPR, ipcruise1, 1],                          # 15. Current BPR
]
# Filter the initial conditions
for i in 1:length(initial)
    if initial[i]>upper[i] || initial[i]<lower[i]
        initial[i] = 0.5*(upper[i]+lower[i])
    end
end
# set tolerances
tol_rel = 1e-6
iters_max_opt = 1000

#### Setup and run optimization
# Set the optimizer
opt = NLopt.Opt(:LN_NELDERMEAD, length(initial))
opt.lower_bounds  = lower
opt.upper_bounds  = upper
opt.min_objective = obj
opt.initial_step  = init_dx
opt.ftol_rel      = tol_rel
opt.maxeval       = iters_max_opt

# Optimize
opt_time = @elapsed begin
    try
        (optf, optx, ret) = NLopt.optimize(opt, initial)
        quicksave_aircraft(ac,joinpath(save_dir, "$(save_Name).jld2"))
    catch e
        println("Optimization failed: $e")
        optf, optx, ret = Inf, initial, :FAILURE
    end
end


