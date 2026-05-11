module OptimizeRangeFuel
export MissionReq, BoundsOpt, ConstraintsOpt, optimize_rangefuel_fun!

using TASOPT, NLopt
include(joinpath(@__DIR__, "objective_factory.jl"))
using .ObjectiveFactory: OptHistory, make_obj, best_feasible
include(__TASOPTindices__)

Base.@kwdef mutable struct MissionReq
    range_des::Float64 = (3000.0 * 1852.0)  #Design flight range (m)
    idx_fuel::Int = 24 #Fuel Index: Jet Fuel(24), Ethanol(32)
    rho_fuel::Float64 = 817.0 #Fuel Density: Jet Fuel(817.0) (kg/m3), Ethanol(789.0) (kg/m3)
    hvap_fuel::Float64 = 358694.0 #Heat of Vaporization: Jet Fuel(358694.0) (J/kg), Ethanol(918187.9) (J/kg)
end

Base.@kwdef mutable struct BoundsOpt #Limites: [Lower Bound, Upper Bound, Initial dx]
    AR_lim::NTuple{3,Float64}             = (6.0,    18.0,   0.5)   #1  Wing aspect ratio
    CL_lim::NTuple{3,Float64}             = (0.45,   0.75,   0.05)  #2  Cruise CL
    sweep_deg_lim::NTuple{3,Float64}      = (25.0,   30.0,   0.1)   #3  Wing sweep angle (deg)
    alt_cruise_m_lim::NTuple{3,Float64}   = (10000.0,20000.0,200.0) #4  Cruise altitude (m)
    taper_in_lim::NTuple{3,Float64}       = (0.65,   0.85,   0.01)  #5  Inboard wing taper ratio
    taper_out_lim::NTuple{3,Float64}      = (0.10,   0.40,   0.01)  #6  Outboard wing taper ratio
    tc_root_lim::NTuple{3,Float64}        = (0.125,  0.15,   0.01)  #7  Inboard thickness-to-chord ratio
    tc_span_lim::NTuple{3,Float64}        = (0.125,  0.15,   0.01)  #8  Outboard thickness-to-chord ratio
    rcls_lim::NTuple{3,Float64}           = (0.90,   1.30,   0.01)  #9  Break/root Cl ratio at cruise
    rclt_lim::NTuple{3,Float64}           = (0.70,   1.00,   0.01)  #10 Tip/root Cl ratio at cruise
    Tt4_lim::NTuple{3,Float64}            = (1400.0, 1650.0, 100.0) #11 Tt4 at cruise (K)
    PR_hpc_lim::NTuple{3,Float64}         = (10.0,   15.0,   0.5)   #12 HPC pressure ratio at cruise
    PR_fan_lim::NTuple{3,Float64}         = (1.25,   2.0,    0.05)  #13 Fan pressure ratio at cruise
    PR_lpc_lim::NTuple{3,Float64}         = (2.999,  3.001,  0.0001)#14 LPC pressure ratio at cruise
    BPR_lim::NTuple{3,Float64}            = (1.0,    20.0,   1.0)   #15 Bypass ratio at cruise
end

Base.@kwdef mutable struct ConstraintsOpt #Constrained values for optimization
    span_max::Float64      = 35.814   # Maximum span (m)
    lenField_max::Float64  = 2400.0   # Maximum balanced field length (m)
    TOCGamma_min::Float64  = 0.015    # Minimum top-of-climb flight angle (rad)
    Tt3_max::Float64       = 900.0    # Maximum combustor inlet temperature (K)
    TMetal_max::Float64    = 1333.33  # Maximum metal temperature (K)
    DiaFan_max::Float64    = 2.0      # Maximum fan diameter (m)
end

"""
optimize_rangefuel_fun!(ac; mission_req=default, bounds_opt=default, constraints_opt=default
                        tol_rel=default, maxeval=default, optimizer=default)

Optimize an aircraft for a design flight range and fuel type
    Input aircraft model will be modified inplace to become the optimized solution
    The input aircraft model, which is already sized, will be use a initial guess
    Mission requirements, including range and fuel will be overwritting the model
    Optimization bounds will be limiting the 15 parameters to be tuned
    6 constraints will be ensured satisfied by the end of optimization
Inputs:
    ac: TASOPT.aircraft: Aircraft model to be optimized and used for initialization (Already sized)
    mission_req: MissionReq: Required range and fuel
    bounds_opt: BoundsOpt: 15 parameters to tunes and their absolute bounds [AR,CL,sweep(deg),altitude,λ_in,λ_out,t/c_root,t/c_span,rcls,rclt,Tt4,π_hc,π_f,π_lc,BPR]
    constraints_opt: ConstraintsOpt: 6 constraints for optimization [max_span,max_lenField,min_TOCGamma,max_Tt3,max_TMetal,max_DiaFan]
    tol_rel:Float64: Relative tolerance for the optmizer to converge
    iters_max_opt:Int: maximum iteration for the optimizers to run
    optimizer_type:Symbol: optimization method
Outputs:
    hist_optim: OptHistory: Optimization history [test_param,penalty,PFEI,violations]
    status: Symbol: status of the optimization [:FTOL_REACHED(Good), :MAXEVAL_REACHED(Okay), :FAILURE(Bad)]
Updates:
    ac
"""
function optimize_rangefuel_fun!(
    ac; 
    mission_req::MissionReq = MissionReq(),
    bounds_opt::BoundsOpt = BoundsOpt(),
    constraints_opt::ConstraintsOpt = ConstraintsOpt(),
    tol_rel::Float64=1e-6, iters_max_opt::Int=1000, 
    optimizer_type::Symbol=:LN_NELDERMEAD)

    #### Size check
    # Model
    ac.is_sized[1] || throw(ArgumentError("Input aircraft model has to be sized first"))
    tol_rel > 0 || throw(ArgumentError("Relative tolerance for optimization must be > 0"))
    iters_max_opt > 0 || throw(ArgumentError("maximum number of optimization interations must be > 0"))
    # Mission requirements
    mission_req.range_des > 0 || throw(ArgumentError("Design flight range > 0"))
    mission_req.rho_fuel > 0 || throw(ArgumentError("Fuel density > 0"))
    mission_req.hvap_fuel >= 0 || throw(ArgumentError("Fuel heat of vaporization >= 0"))
    # Constraints
    constraints_opt.span_max > 0 || throw(ArgumentError("span_max > 0")) # Maximum span (m)
    constraints_opt.lenField_max > 0 || throw(ArgumentError("lenField_max > 0")) # Maximum balanced field length (m)
    constraints_opt.TOCGamma_min > 0 || throw(ArgumentError("TOCGamma_min > 0")) # Minimum top-of-climb flight angle (rad)
    constraints_opt.Tt3_max > 0 || throw(ArgumentError("Tt3_max > 0")) # Maximum combustor inlet temperature (K)
    constraints_opt.TMetal_max > 0 || throw(ArgumentError("TMetal_max > 0")) # Maximum metal temperature (K)
    constraints_opt.DiaFan_max > 0 || throw(ArgumentError("DiaFan_max > 0")) # Maximum fan diameter (m)

    #### Modify the mission requirement for optimization
    ac.parm[imRange,:] .= mission_req.range_des #Design flight range (m)
    ac.options.ifuel = mission_req.idx_fuel #Fuel Index
    ac.parg[igrhofuel] = mission_req.rho_fuel #Fuel Density (kg/m3)
    ac.pare[iehvap, :, :] .= mission_req.hvap_fuel #Heat of Vaporization (J/kg)
    ac.pare[iehvapcombustor, :, :] .= mission_req.hvap_fuel

    #### Create an objective function
    # Initialize a history of optimization
    hist_optim = OptHistory() #Optimization history
    # Setup the constraints
    constraints_optim = [constraints_opt.span_max, constraints_opt.lenField_max, constraints_opt.TOCGamma_min,
    constraints_opt.Tt3_max, constraints_opt.TMetal_max, constraints_opt.DiaFan_max]
    # Get an objective function
    obj = make_obj(ac, constraints_optim, hist_optim) #This obj modify the aircraft model in place

    #### Setup the optimized parameters
    # Construct the bounds
    numBounds = fieldcount(BoundsOpt)
    lower   = Vector{Float64}(undef, numBounds)
    upper   = Vector{Float64}(undef, numBounds)
    init_dx = Vector{Float64}(undef, numBounds)
    for (idx, fNam) in enumerate(fieldnames(BoundsOpt))
        lower_cur, upper_cur, init_dx_cur = getfield(bounds_opt, fNam)
        # Check error
        lower_cur < upper_cur || throw(ArgumentError("Optimization parameters: $(fNam): lower bound must be < upper bound"))
        init_dx_cur > 0  || throw(ArgumentError("Optimization initial step for: $(fNam): must be > 0"))
        # Form vector
        lower[idx]   = lower_cur
        upper[idx]   = upper_cur
        init_dx[idx] = init_dx_cur
    end
    # Use the initial conditions from aircraft model (Reused previous optimized solution)
    initial = [
    ac.wing.layout.AR,                                     #1  Wing aspect ratio
    ac.para[iaCL, ipcruise1, 1],                           #2  Cruise CL
    ac.wing.layout.sweep,                                  #3  Wing sweep angle (deg)
    ac.para[iaalt, ipcruise1, 1],                          #4  Cruise altitude (m)
    ac.wing.inboard.λ,                                     #5  Inboard wing taper ratio
    ac.wing.outboard.λ,                                    #6  Outboard wing taper ratio
    ac.wing.inboard.cross_section.thickness_to_chord,      #7  Inboard thickness-to-chord ratio
    ac.wing.outboard.cross_section.thickness_to_chord,     #8  Outboard thickness-to-chord ratio
    ac.para[iarcls, ipcruise1, 1],                         #9  Break/root Cl ratio at cruise
    ac.para[iarclt, ipcruise1, 1],                         #10 Tip/root Cl ratio at cruise
    ac.pare[ieTt4, ipcruise1, 1],                          #11 Tt4 at cruise (K)
    ac.pare[iepihc, ipcruise1, 1],                         #12 HPC pressure ratio at cruise
    ac.pare[iepif, ipcruise1, 1],                          #13 Fan pressure ratio at cruise
    ac.pare[iepilc, ipcruise1, 1],                         #14 LPC pressure ratio at cruise
    ac.pare[ieBPR, ipcruise1, 1],                          #15 Bypass ratio at cruise
    ]
    # Filter the out of bound initial conditions
    initial_saved = initial #save an unclamped initial to retract back the aircraft model state
    initial = clamp.(initial, lower, upper)

    #### Setup the optimizer
    opt = NLopt.Opt(optimizer_type, numBounds)
    opt.lower_bounds  = lower
    opt.upper_bounds  = upper
    opt.min_objective = obj
    opt.initial_step  = init_dx
    opt.ftol_rel      = tol_rel
    opt.maxeval       = iters_max_opt

    #### Run the optimization
    try
        (_, _, status) = NLopt.optimize(opt, initial)
    catch errOpt
        println("Current O=optimization failed: $(errOpt)")
        (_, _, status) = (Inf, initial, :FAILURE)
    end

    #### Prepare an updated aircraft model state
    if status != :FAILURE
        bestSol = best_feasible(hist_optim)
        if isnothing(bestSol)
            println("No constraint satisfying solutions found, return to starting state")
            obj(initial_saved, Float64[]) #This likely wont be feasible either
        else
            println("Identified and setup the current best solution with optimization reached by $(status)")
            obj(bestSol.test_param, Float64[])
        end
    else
        println("Optimization process failed, return to starting state")
        obj(initial_saved, Float64[]) #This likely wont be feasible either
    end

    return status, hist_optim
end

end # module OptimizeRangeFuel