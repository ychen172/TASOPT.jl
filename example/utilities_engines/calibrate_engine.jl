module CaliEng

using TASOPT
using NLopt
include(__TASOPTindices__)
import ..RunEngine
import ..ObjectiveFactory

"""
    make_obj(ac, parameters::AbstractVector{<:ObjectiveFactory.Parameter}, Fn_N::AbstractVector{<:Real},
            WFuel_kgs::AbstractVector{<:Union{Real,Missing}}, OPR::AbstractVector{<:Union{Real,Missing}}, BPR::AbstractVector{<:Union{Real,Missing}};
            M0::Real=0.0, P0::Real=101320.0, T0::Real=288.2, a0::Real=340.21, dia_fan_m::Real=-1.0, pen_wei::AbstractVector{<:Real}=[1.0/3.0,1.0/3.0,1.0/3.0],
            print_every::Int=10, max_iter_sizing::Int=150, pen_failed_sizing::Float64=10000.0, pen_failed_engine::Float64=300.0)

`make_obj` constructs and returns an objective function for engine calibration that will be fed to optimizer

!!!! details "🔃 Inputs and Outputs"
    **Inputs:**
        - `ac`: Baseline aircraft model for sizing (Passed by copy)
        - `parameters`: Parameters for optimization. (SI unit) (Unaltered by obj! call)
        - `Fn_N`: Thrust for engine mission (N)
        - `WFuel_kgs`: Fuel flow rate for engine mission (kg/s) (Same length as Fn_N, but accept Missing)
        - `OPR`: Overall pressure ratio for engine mission (Same length as Fn_N, but accept Missing)
        - `BPR`: Bypass ratio for engine mission (Same length as Fn_N, but accept Missing)
        - `M0`: Inlet engine Mach number
        - `P0`: Inlet engine pressure (Pa)
        - `T0`: Inlet engine temperature (K)
        - `a0`: Inlet speed of sound (m/s)
        - `dia_fan_m`: Reference fan diameter to match (m) set negative to disable this criterion
        - `pen_wei`: [Weight_mdotFuel, Weight_OPR, Weight_BPR]
        - `print_every`: Every number of iterations before print out for the current optimized state (Performed by obj! call)
        - `max_iter_sizing`: Maximum number of iteration used for the sizing loop
        - `pen_failed_sizing`: Penalty value for failed aircraft sizing process
        - `pen_failed_engine`: Penalty value for failed engine operation point
    
    **Outputs:**
        - `obj!`::a function: Objective function for optimizer. Only return a penalty value and modify the optimization history. 
                  A copy of aircraft model is changed within obj! call. It does not affect the model, initially inputted into make_obj function. 
                  The initially inputted model is used to refresh the aircraft model used by obj! before every call for deterministic behavior.
        - `hist`::OptHistory: Optimization history for only succeed sizing case, containing test_param, penalty, PFEI, violations. returned by reference and will be modified by obj! call.

!!! note "Behavior"
        - Penalty is defined for the deviation from the deviation of fuel flow rate, OPR, and BPR
        - Call to `obj!` will modify the returned `hist` by reference. However, the `ac` changed within obj! call does not affect the ac for make_obj input or the ac outside of the function call. (passed and reset by copy)
        - In the outer loop, if parameters length got changed, then the number of parameters expected by obj! here is also altered correspondingly. Beware.
        - All inputs are in SI unit
"""
function make_obj(ac, parameters::AbstractVector{<:ObjectiveFactory.Parameter}, Fn_N::AbstractVector{<:Real},
                  WFuel_kgs::AbstractVector{<:Union{Real,Missing}}, OPR::AbstractVector{<:Union{Real,Missing}}, BPR::AbstractVector{<:Union{Real,Missing}};
                  M0::Real=0.0, P0::Real=101320.0, T0::Real=288.2, a0::Real=340.21, dia_fan_m::Real=-1.0, pen_wei::AbstractVector{<:Real}=[1.0/3.0,1.0/3.0,1.0/3.0],
                  print_every::Int=10, max_iter_sizing::Int=150, pen_failed_sizing::Float64=10000.0, pen_failed_engine::Float64=300.0)
    #### Size check
    isempty(parameters)   && throw(ArgumentError("`parameters` cannot be empty. Provide at least one optimization variable."))
    ((length(Fn_N) >= 1)  && (length(WFuel_kgs) >= 1) && (length(OPR) >= 1) && (length(BPR) >= 1)) || throw(ArgumentError("more than one mission point needs to be provided"))
    ((length(Fn_N) == length(WFuel_kgs)) && (length(WFuel_kgs) == length(OPR)) && (length(OPR) == length(BPR))) || throw(ArgumentError("Thrust, fuel flow, OPR, and BPR need to provided for each mission point"))
    length(pen_wei) == 3  || throw(ArgumentError("three weighting needs to be provided for fuel flow rate, OPR, and BPR"))
    print_every >= 1      || throw(ArgumentError("`print_every` must be ≥ 1"))
    max_iter_sizing >= 10 || throw(ArgumentError("`max_iter_sizing` must be ≥ 10"))
    pen_failed_sizing > 0 || throw(ArgumentError("`pen_failed_sizing` must be > 0"))
    
    #### Construct an optimization history to be returned by reference for successive update
    hist = ObjectiveFactory.OptHistory() #History with empty entries

    #### make a backup
    ac_bak = deepcopy(ac)
    
    #### make a iteration counter
    iterCount = 0
    num_WFuel = count(.!ismissing.(WFuel_kgs))
    num_OPR = count(.!ismissing.(OPR))
    num_BPR = count(.!ismissing.(BPR))
    #### Construct the objective function
    function obj!(x, grad)
        #### Always reset back to the same initial aircraft model for the current optimization attempt
        ac = deepcopy(ac_bak)

        #### Update current function call count
        iterCount += 1
        println("Current Objective Function Call # $(iterCount)")

        #### Sanity check
        (grad === nothing || isempty(grad)) || throw(ArgumentError("Gradient requested, but obj! has no gradient implementation. Use derivative-free NLopt algorithms (NM/DIRECT)."))
        (length(x) == length(parameters)) || throw(DimensionMismatch("#Param in by optimizer $(length(x)) not equal to the #Param expected $(length(parameters))"))

        #### Insert the parameters into the aircraft model 
        for (i_param, param_cur) in enumerate(parameters)
            ObjectiveFactory.setNestedProp_fromExpr!(ac, x[i_param]; field_path=param_cur.field_path, index=param_cur.index)
        end

        #### Test sizing the current case
        flgSuccessed = true
        penalty = 0.0
        try
            TASOPT.size_aircraft!(ac, iter=max_iter_sizing, printiter=false)
            if dia_fan_m > 0.0
                penalty += 100.0*abs((ac.parg[igdfan]-dia_fan_m)/dia_fan_m)
            end
        catch e
            if e isa InterruptException #Unless user interruption
                rethrow()
            end
            flgSuccessed = false #Continue but mark current as failed sizing attempt
        end

        #### Primary penalty value
        penalty += flgSuccessed ? 0.0 : pen_failed_sizing #(J/J)
        
        #### Run engine mission
        for (idxFn,Fn_N_cur) in enumerate(Fn_N)
            try
                res = RunEngine.runOffDes(ac, M0, P0, T0, a0, Fn_N_cur)
                penalty += !res.Lconv ? pen_failed_engine : 0.0
                if !ismissing(WFuel_kgs[idxFn])
                    penalty += 100*(abs(res.mcore*res.ff-WFuel_kgs[idxFn])/WFuel_kgs[idxFn])*(pen_wei[1]/num_WFuel)
                end
                if !ismissing(OPR[idxFn])
                    penalty += 100*(abs(res.OPR-OPR[idxFn])/OPR[idxFn])*(pen_wei[2]/num_OPR)
                end
                if !ismissing(BPR[idxFn])
                    penalty += 100*(abs(res.BPR-BPR[idxFn])/BPR[idxFn])*(pen_wei[3]/num_BPR)
                end
            catch e
                if e isa InterruptException #Unless user interruption
                    rethrow()
                end
                penalty += pen_failed_engine
            end
        end

        #### Store history only if sizing were successful
        if flgSuccessed
            # Update the Optimization History
            push!(hist.test_param, copy(x))
            push!(hist.penalty, penalty)
            push!(hist.PFEI, 0.0)
            push!(hist.violations, Vector{ObjectiveFactory.Constraint}())
        end

        #### In-flight Printout
        ObjectiveFactory.InFlightPrintOutBest(hist, print_every, iterCount)
        
        return penalty
    end
    return (; obj!, hist)
end

"""
    optimize_match_EEDB!(ac,parameters::AbstractVector{<:ObjectiveFactory.Parameter},
                         Fn_N::AbstractVector{<:Real},WFuel_kgs::AbstractVector{<:Union{Real,Missing}}, 
                         OPR::AbstractVector{<:Union{Real,Missing}}, BPR::AbstractVector{<:Union{Real,Missing}};
                         M0::Real=0.0, P0::Real=101320.0, T0::Real=288.2, a0::Real=340.21, dia_fan_m::Real=-1.0, 
                         pen_wei::AbstractVector{<:Real}=[1.0/3.0,1.0/3.0,1.0/3.0],
                         print_every::Int=10, max_iter_sizing::Int=150, 
                         pen_failed_sizing::Float64=100.0, pen_failed_engine::Float64=3.0,
                         ftol_rel::Float64=1e-6,
                         max_iter_optim::Int=1000,
                         optimizer_type::Symbol=:LN_NELDERMEAD)

`optimize_match_EEDB!` Single round optimization such that the engine performance can match well with the reference

    ***Inputs:***
        - `ac`: Baseline aircraft model for sizing (Passed by copy)
        - `parameters`: Parameters for optimization. (SI unit) (Unaltered by obj! call)
        - `Fn_N`: Thrust for engine mission (N)
        - `WFuel_kgs`: Fuel flow rate for engine mission (kg/s) (Same length as Fn_N, but accept Missing)
        - `OPR`: Overall pressure ratio for engine mission (Same length as Fn_N, but accept Missing)
        - `BPR`: Bypass ratio for engine mission (Same length as Fn_N, but accept Missing)
        - `M0`: Inlet engine Mach number
        - `P0`: Inlet engine pressure (Pa)
        - `T0`: Inlet engine temperature (K)
        - `a0`: Inlet speed of sound (m/s)
        - `dia_fan_m`: Fan diameter to match (m) Set negative to disable this criterion
        - `pen_wei`: [Weight_mdotFuel, Weight_OPR, Weight_BPR]
        - `print_every`: Every number of iterations before print out for the current optimized state (Performed by obj! call)
        - `max_iter_sizing`: Maximum number of iteration used for the sizing loop
        - `pen_failed_sizing`: Penalty value for failed aircraft sizing process
        - `pen_failed_engine`: Penalty value for failed engine operation point
        - `ftol_rel`: relative tolerance to converge in optimization
        - `max_iter_optim`: maixmum optimization iterations (At least 3)
        - `optimizer_type`: optimization method: Ex. Local Search: :LN_NELDERMEAD or Global Search: :GN_CRS2_LM, :GN_DIRECT. (Only non-gradient base supported)
    
    ***Outputs***
        - status: Status of optimization: 
          Success cases: (:RECOVERED_FEASIBLE_FROM_HISTORY, :SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED, :MAXEVAL_REACHED, :MAXTIME_REACHED)
          Fail cases: (:NO_FEASIBLE_SOLUTION, :FAILURE, :INVALID_ARGS, :OUT_OF_MEMORY, :ROUNDOFF_LIMITED, :FORCED_STOP)
        - hist: Optimization history: {test_param::Vector{Vector{T}}, penalty::Vector{T}, PFEI::Vector{T}, violations::Vector{Vector{Constraint}}}
        - bestSol: Best solution from the history: either nothing or one element from OptHistory variable
    
    ***Updated***
        - parameters: The value `:val` of the parameters are updated to the best feasible solution found if optimization succeeded. `:val` is still the initial value if optimization failed.

    ***Behavior***
        - Has to sized the aircraft beforehead
        - parameters and mission requirements needs to be in SI unit except for (deg used by sweep angle)
        - parameters and mission requirements needs to be constructed using special constructor but not manual
        - Only support non-gradient base optimization method
        - Always unchange the mission requirements
        - Update the value inside parameters, status, and history of optimization

"""
function optimize_match_EEDB!(ac,parameters::AbstractVector{<:ObjectiveFactory.Parameter},
                              Fn_N::AbstractVector{<:Real},WFuel_kgs::AbstractVector{<:Union{Real,Missing}}, 
                              OPR::AbstractVector{<:Union{Real,Missing}}, BPR::AbstractVector{<:Union{Real,Missing}};
                              M0::Real=0.0, P0::Real=101320.0, T0::Real=288.2, a0::Real=340.21, dia_fan_m::Real=-1.0,
                              pen_wei::AbstractVector{<:Real}=[1.0/3.0,1.0/3.0,1.0/3.0],
                              print_every::Int=10, max_iter_sizing::Int=150, 
                              pen_failed_sizing::Float64=100.0, pen_failed_engine::Float64=3.0,
                              ftol_rel::Float64=1e-6,
                              max_iter_optim::Int=1000,
                              optimizer_type::Symbol=:LN_NELDERMEAD)
    #### make a copy of the input aircraft model state to ensure model is unaltered by function call for successive deterministic behavior
    ac = deepcopy(ac)
    
    #### Size check (Only check parameters that are not checked by sub functions)
    ftol_rel > 0       || throw(ArgumentError("Relative tolerance for optimization must be > 0"))
    max_iter_optim > 2 || throw(ArgumentError("maximum number of optimization interations must be > 2"))
    
    #### Create an objective function and optimization history container
    (; obj!, hist) = make_obj(ac, parameters, Fn_N, WFuel_kgs, OPR, BPR; M0=M0, P0=P0, T0=T0, a0=a0, dia_fan_m=dia_fan_m, pen_wei=pen_wei,
                              print_every=print_every, max_iter_sizing=max_iter_sizing, pen_failed_sizing=pen_failed_sizing, pen_failed_engine=pen_failed_engine)
    
    #### Setup the optimized parameters
    # Reshape the optimization parameters
    upper_bounds = Float64.(getfield.(parameters, :bon_up)) #Boardcasting extraction by copy
    lower_bounds = Float64.(getfield.(parameters, :bon_lo))
    initial_step = Float64.(getfield.(parameters, :d_val))
    initial_gues = Float64.(getfield.(parameters, :val))
    number_param = length(parameters) #Number of parameters use to intiailize optimizer
    # Chek again the bounds and initial guess (In case the bounds and intial values are updated without the use of the constructor)
    any(upper_bounds .<= lower_bounds) && error("Some upper bounds are below or equal to lower bounds")
    any(initial_step .<= 0) && error("Some initial step sizes are non-positive")
    if any((initial_gues .< lower_bounds) .| (initial_gues .> upper_bounds))
        @warn "Some initial guesses are out of bounds. Clamping them to within bounds."
        initial_gues = clamp.(initial_gues, lower_bounds, upper_bounds)
    end
    
    #### Optimization
    status = :FAILURE #Initialization
    try
        # Setup the optimizer
        opt               = NLopt.Opt(optimizer_type, number_param)
        opt.lower_bounds  = lower_bounds
        opt.upper_bounds  = upper_bounds
        opt.min_objective = obj!
        opt.initial_step  = initial_step
        opt.ftol_rel      = ftol_rel
        opt.maxeval       = max_iter_optim
        # Runing
        (_, _, status)    = NLopt.optimize(opt, initial_gues)
    catch e
        e isa InterruptException && rethrow()
        @error "Current optimization failed $(typeof(e)): $(e)"
    end

    #### Post-process the optimization results
    # Extract the best feasible solution
    bestSol = ObjectiveFactory.best_feasible(hist) #Name tuple of `index`, `test_param`, `penalty`, `PFEI`
    if isnothing(bestSol)
        println("No constraint satisfying solutions found, \nOptimize_par is left at its initial state. \nAircraft model is insid the current function call may change but NOT changing the input model to the function due to passed by copy")
        if status in ObjectiveFactory.success_statuses #Ex. FTOL is reached but no feasible solution were found.
            status = :NO_FEASIBLE_SOLUTION
        end
    else
        println("Identified and setup the current best solution.")
        if status in ObjectiveFactory.failure_statuses #Ex. Maximum memory reached, but some feasible constraint satisfying solution have beed tested.
            @warn "Using feasible solution from history despite optimizer termination status = $(status)"
            status = :RECOVERED_FEASIBLE_FROM_HISTORY
        end
        ObjectiveFactory.setproperty!.(parameters, :val, bestSol.test_param) #This set the test parameters from the best run back into the optimization parameter container (Use as mild initial guess for next case)
    end

    return status, hist, bestSol
end

"""
    UpdAcEngMod!(ac_ref, x; tol::Float64=1e-6, maxIter::Int=150)

`UpdAcEngMod!` sets the engine-cycle design variables on `ac_ref` and re-sizes its engine at the
cruise design point, iterating the cruise cycle sizing against the takeoff-rotation cooling-flow
sizing until TSFC settles.

    ***Inputs:***
        - `ac_ref`: Aircraft model to update (mutated in place)
        - `x::AbstractVector{<:Real}`: `[BPR, pif, pilc, pihc, Tt4]` candidate cycle design variables (SI units)
        - `tol::Float64`: Relative tolerance on TSFC for the design/cooling_sizing coupling loop to be considered converged
        - `maxIter::Int`: Maximum number of design/cooling_sizing alternations (same convention as `size_aircraft!`'s own weight-closure loop)

    ***Outputs***
        - `ac_ref`: The same aircraft model, mutated in place and returned for convenience

    ***Behavior***
        - Overwrites `ac_ref.pare[ieBPR/iepif/iepilc/iepihc/ieTt4, ipcruise1, 1]` with `x`
        - Alternates `"design"` (cruise cycle sizing, `ipcruise1`) and `"cooling_sizing"` (Tmetal-driven
          cooling flow at `iprotate`) engine calcs -- the two are coupled through the shared design-reference
          values (`pifD` etc.) and the cooling-flow ratio (`epsrow`), so a single pass of each is NOT enough
          for the two to be mutually consistent
        - Uses `initializes_engine=true` on every call, so each evaluation starts from a clean bootstrap
          rather than warm-starting off whatever was previously in `ac_ref`
        - Does not throw on non-convergence of the underlying engine calcs -- caller is responsible for
          checking/catching (see `make_obj_engine_opt`)
"""
function UpdAcEngMod!(ac_ref, x; tol::Float64=1e-6, maxIter::Int=150)
    BPR,pif,pilc,pihc,Tt4 = x
    ac_ref.pare[ieBPR,ipcruise1,1] = BPR
    ac_ref.pare[iepif,ipcruise1,1] = pif
    ac_ref.pare[iepilc,ipcruise1,1] = pilc
    ac_ref.pare[iepihc,ipcruise1,1] = pihc
    ac_ref.pare[ieTt4,ipcruise1,1] = Tt4

    TSFC_prev = Inf
    for _ in 1:maxIter
        ac_ref.engine.enginecalc!(ac_ref, "design", 1, ipcruise1, true, 1)
        ac_ref.engine.enginecalc!(ac_ref, "cooling_sizing", 1, iprotate, true, 1)
        TSFC_cur = ac_ref.pare[ieTSFC, ipcruise1, 1]
        if abs(TSFC_cur - TSFC_prev) < tol*max(abs(TSFC_prev), 1e-30)
            break
        end
        TSFC_prev = TSFC_cur
    end
    return ac_ref
end

"""
    make_obj_engine_opt(ac, printEvery::Int64)

`make_obj_engine_opt` constructs and returns an objective function for the engine-cycle-only
inner-loop optimization: minimize TSFC at a single fixed flight condition/thrust (the design point
already stored on `ac`), searching over `[BPR, pif, pilc, pihc, Tt4]` via `UpdAcEngMod!`.

    ***Inputs:***
        - `ac`: Baseline aircraft model, providing the fixed flight condition/thrust/technology
          parameters for the design point (passed by copy, never mutated)
        - `printEvery::Int64`: Print the current best solution every this many objective calls
        - `tol_coupling::Float64`: Forwarded to `UpdAcEngMod!` as its design/cooling_sizing
          coupling-loop tolerance -- should be tighter than whatever NLopt tolerance the caller's
          own optimizer uses (e.g. `engine_opt`'s `ftol`), since that search can't reliably resolve
          improvements finer than the noise floor left by an under-converged coupling loop

    ***Outputs***
        - `obj!`::a function: Objective function `(x, grad) -> penalty` for the optimizer. A fresh
          copy of `ac` is used and updated via `UpdAcEngMod!` on every call; `ac` itself is never
          modified.
        - `histPara`::Vector{Vector{Float64}}: Parameter vectors `x` for every successful evaluation,
          returned by reference and updated by `obj!` calls
        - `histPenl`::Vector{Float64}: Penalty (TSFC/gee) for every successful evaluation, same order
          as `histPara`, returned by reference and updated by `obj!` calls

    ***Behavior***
        - Penalty is `TSFC/gee` (kg/s/N); if the engine calc fails to converge or a non-positive TSFC
          results, a fixed fallback penalty (1.78e-3) is returned instead and the evaluation is NOT
          recorded in `histPara`/`histPenl`
        - All inputs/outputs are in SI units
"""
function make_obj_engine_opt(ac,printEvery::Int64;tol_coupling::Float64=1e-8)
    ac_used = deepcopy(ac)
    histPara_engine_opt = Vector{Vector{Float64}}()
    histPenl_engine_opt = Vector{Float64}()
    count = 0
    function obj_engine_opt!(x, grad)
        count += 1
        ac_ref = deepcopy(ac_used)
        penal=0.0
        try
            UpdAcEngMod!(ac_ref, x; tol=tol_coupling)
            penal = ac_ref.pare[ieTSFC, ipcruise1, 1] / gee #(kg/s/N) ~1.78e-5
            if penal<=0.0
                penal = 1.78e-3
            else
                push!(histPara_engine_opt,copy(x))
                push!(histPenl_engine_opt,penal)
            end
        catch e
            e isa InterruptException && rethrow()
            penal = 1.78e-3
        end
        if (mod(count,printEvery)==0)&&(length(histPenl_engine_opt)>0)
            idxMin = argmin(histPenl_engine_opt)
            println("EngRun#$(count): Current best for obj_engine_opt: Penalty: $(histPenl_engine_opt[idxMin]) with Parameters: $(histPara_engine_opt[idxMin])")
        else
            println("EngRun#$(count): No Feasible Sol Yet")
        end
        return penal
    end
    return (;obj! = obj_engine_opt!,histPara = histPara_engine_opt, histPenl = histPenl_engine_opt)
end

"""
    engine_opt(ac; ini::Vector{Float64}, upBon::Vector{Float64}, loBon::Vector{Float64},
              printEvery::Int64, ftol::Float64=1e-6, tol_coupling::Float64=1e-8, maxIter::Int=1000, optTyp::Symbol=:LN_NELDERMEAD)

`engine_opt` runs the engine-cycle-only inner-loop optimization (see `make_obj_engine_opt`) to find
the minimum-TSFC engine cycle `[BPR, pif, pilc, pihc, Tt4]` at `ac`'s fixed design point.

    ***Inputs:***
        - `ac`: Baseline aircraft model for sizing (passed by copy)
        - `ini::Vector{Float64}`: Initial guess `[BPR, pif, pilc, pihc, Tt4]` (SI units)
        - `upBon::Vector{Float64}`: Upper bounds, same order as `ini`
        - `loBon::Vector{Float64}`: Lower bounds, same order as `ini`
        - `printEvery::Int64`: Print the current best solution every this many objective calls
        - `ftol::Float64`: Relative tolerance for this function's own NLopt optimizer to converge
        - `tol_coupling::Float64`: Forwarded to `make_obj_engine_opt`/`UpdAcEngMod!` as the design/
          cooling_sizing coupling-loop tolerance -- should be tighter than `ftol`, since this NLopt
          search can't reliably resolve improvements finer than the noise floor left by an
          under-converged coupling loop
        - `maxIter::Int`: Maximum number of optimization iterations
        - `optTyp::Symbol`: NLopt algorithm symbol (only non-gradient-based algorithms supported,
          e.g. `:LN_NELDERMEAD` local, `:GN_CRS2_LM`/`:GN_DIRECT` global)

    ***Outputs***
        - `bestSol::Vector{Float64}`: `[BPR, pif, pilc, pihc, Tt4]` of the minimum-penalty evaluation
          found, or a vector of `NaN` if the optimizer status wasn't a success or no evaluation
          succeeded
        - `status::Symbol`: `:NO_FEASIBLE_SOLUTION` if the optimizer reported success but no
          evaluation succeeded; otherwise NLopt's own termination status
        - `histPara::Vector{Vector{Float64}}`: All successful evaluations' parameter vectors
        - `histPenl::Vector{Float64}`: All successful evaluations' penalties, same order as `histPara`

    ***Behavior***
        - Errors immediately if `ini` is outside `[loBon,upBon]`
        - `ac` is never mutated -- a fresh copy is made internally
        - Only supports non-gradient-based NLopt algorithms
"""
function engine_opt(ac;
                    ini::Vector{Float64},upBon::Vector{Float64},loBon::Vector{Float64},
                    printEvery::Int64,ftol::Float64=1e-6,tol_coupling::Float64=1e-8,maxIter::Int=1000,optTyp::Symbol=:LN_NELDERMEAD)
    any((ini .< loBon) .| (ini .> upBon)) && error("Initial guess `ini` is outside the bounds [loBon,upBon]: ini=$(ini), loBon=$(loBon), upBon=$(upBon)")
    ac_used = deepcopy(ac)
    (; obj!, histPara, histPenl) = make_obj_engine_opt(ac_used, printEvery; tol_coupling=tol_coupling)
    status = :FAILURE
    try
        opt               = NLopt.Opt(optTyp, length(ini))
        opt.lower_bounds  = loBon
        opt.upper_bounds  = upBon
        opt.min_objective = obj!
        opt.initial_step  = (upBon .- loBon)*0.1
        opt.ftol_rel      = ftol
        opt.maxeval       = maxIter
        # Runing
        (_, _, status)    = NLopt.optimize(opt, ini)
    catch e
        e isa InterruptException && rethrow()
        @error "Current optimization failed $(typeof(e)): $(e)"
    end

    #### Post-process the optimization results
    bestSol = ini.*NaN
    if (status in ObjectiveFactory.success_statuses) && length(histPenl)>0
        idxMin = argmin(histPenl)
        bestSol = histPara[idxMin]
    else
        status = :NO_FEASIBLE_SOLUTION
    end
    return bestSol, status, histPara, histPenl
end

"""
    UpdAcTecLvl!(ac_ref, x::Vector{Float64}, ini_eng::Vector{Float64},
                upBon_eng::Vector{Float64}, loBon_eng::Vector{Float64}; printEvery::Int64=10,
                ftol_eng::Float64=1e-7, tol_coupling::Float64=1e-8, maxIter::Int=1000, optTyp::Symbol=:LN_NELDERMEAD)

`UpdAcTecLvl!` sets the 8 engine technology parameters on `ac_ref`, then runs `engine_opt` to find
the minimum-TSFC engine cycle for that technology level (see `make_obj_engine_opt`/`engine_opt`),
and re-materializes that optimal cycle onto `ac_ref` via `UpdAcEngMod!`.

    ***Inputs:***
        - `ac_ref`: Aircraft model to update (mutated in place)
        - `x::Vector{Float64}`: `[pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal]` candidate technology parameters (SI units)
        - `ini_eng::Vector{Float64}`: Initial guess `[BPR,pif,pilc,pihc,Tt4]` for the inner `engine_opt` cycle-design search
        - `upBon_eng::Vector{Float64}`: Upper bounds for the inner search, same order as `ini_eng`
        - `loBon_eng::Vector{Float64}`: Lower bounds for the inner search, same order as `ini_eng`
        - `printEvery::Int64`: Forwarded to `engine_opt` -- print the current best solution every this many objective calls
        - `ftol_eng::Float64`: Forwarded to `engine_opt` as its NLopt relative tolerance (the middle of the
          three nested loops -- should be looser than `tol_coupling` but tighter than whatever tolerance
          the caller's own outer optimizer uses, e.g. `tech_opt`'s `ftol_tec`)
        - `tol_coupling::Float64`: Forwarded to `UpdAcEngMod!` as its design/cooling_sizing coupling-loop
          tolerance (the innermost of the three nested loops -- should be tighter than `ftol_eng`, since
          `engine_opt`'s own NLopt search can't reliably resolve improvements finer than the noise floor
          left by an under-converged coupling loop)
        - `maxIter::Int`: Forwarded to `engine_opt` as its maximum optimization iterations
        - `optTyp::Symbol`: Forwarded to `engine_opt` as its NLopt algorithm symbol

    ***Outputs***
        - `ac_ref`: The same aircraft model, mutated in place and returned for convenience
        - `flgSizSuc::Bool`: `true` only if the inner `engine_opt` search found a feasible cycle AND
          `UpdAcEngMod!` successfully re-converged it onto `ac_ref`; `false` otherwise
        - `bestSol::Vector{Float64}`: `[BPR,pif,pilc,pihc,Tt4]` found by `engine_opt`, or a vector of
          `NaN` if that search itself failed

    ***Behavior***
        - Broadcasts the 7 `pare`-based technology parameters across every mission point/mission
          (`pare[ie...,:,:] .= ...`) rather than just `ipcruise1` -- these are constant engine-hardware
          characteristics read directly at whatever point is being processed (`tfcalc!` has no
          `ipcruise1` fallback), so setting only one point would leave e.g. `iprotate`'s `cooling_sizing`
          call reading stale values
        - `Tmetal` (`parg[igTmetal]`) is a single aircraft-wide scalar, set once (no per-point broadcast needed)
        - Never throws except `InterruptException` -- any other failure (in either `engine_opt` or the
          `UpdAcEngMod!` re-solve) is caught and reported via `flgSizSuc=false`, matching this file's
          established "always same return shape, never throw" convention
"""
function UpdAcTecLvl!(ac_ref,x::Vector{Float64},ini_eng::Vector{Float64},
                      upBon_eng::Vector{Float64},loBon_eng::Vector{Float64};printEvery::Int64=10,
                      ftol_eng::Float64=1e-7,tol_coupling::Float64=1e-8,maxIter::Int=1000,optTyp::Symbol=:LN_NELDERMEAD)
    # Unpack
    length(x) == 8 || throw(DimensionMismatch("`x` must have exactly 8 elements [pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal], got $(length(x))"))
    pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal = x
    # Update aircraft with technology parameters
    ac_ref.pare[iepib,   :, :] .= pib
    ac_ref.pare[ieepolf, :, :] .= epolf
    ac_ref.pare[ieepollc,:, :] .= epollc
    ac_ref.pare[ieepolhc,:, :] .= epolhc
    ac_ref.pare[ieepolht,:, :] .= epolht
    ac_ref.pare[ieepollt,:, :] .= epollt
    ac_ref.pare[ieetab,  :, :] .= etab
    ac_ref.parg[igTmetal]       = Tmetal
    # Optimize engine at the design point
    flgSizSuc = true
    bestSol = fill(NaN, length(ini_eng))
    try
        bestSol, _, _, _ = engine_opt(ac_ref;ini=ini_eng,upBon=upBon_eng,loBon=loBon_eng,
                                              printEvery=printEvery,ftol=ftol_eng,tol_coupling=tol_coupling,maxIter=maxIter,optTyp=optTyp)
    catch e
        e isa InterruptException && rethrow()
        flgSizSuc = false
    end
    if any(isnan,bestSol)
        flgSizSuc = false
    elseif flgSizSuc
        try
            # Update aircraft with update engine cycle
            UpdAcEngMod!(ac_ref, bestSol; tol=tol_coupling, maxIter=150) #Update the engine model
        catch e
            e isa InterruptException && rethrow()
            flgSizSuc = false
        end
    end
    return ac_ref, flgSizSuc, bestSol
end

"""
    make_obj_tech_cali(ac, ini_eng::Vector{Float64}, upBon_eng::Vector{Float64}, loBon_eng::Vector{Float64};
                       printEvery::Int64, ftol_eng::Float64=1e-7, tol_coupling::Float64=1e-8, maxIter::Int=1000, optTyp::Symbol=:LN_NELDERMEAD,
                       Fn_N::Vector{Float64}, WFuel_kgs_ref::Vector{<:Union{Missing,Float64}},
                       OPR_ref::Vector{<:Union{Missing,Float64}}, BPR_ref::Vector{<:Union{Missing,Float64}},
                       DFan_m_ref::Float64=-1.0, M0::Float64=0.0, P0::Float64=101320.0, T0::Float64=288.2, a0::Float64=340.21)

`make_obj_tech_cali` constructs and returns the outer-loop objective function for engine technology
calibration: for a candidate set of 8 technology parameters, find the minimum-TSFC engine cycle for
that technology (via `UpdAcTecLvl!`), then fly it through a set of EEDB reference thrust levels
(via `RunEngine.runOffDes`) and penalize deviation from the reference fuel flow/OPR/BPR (and,
optionally, fan diameter).

    ***Inputs:***
        - `ac`: Baseline aircraft model, providing the fixed cruise design point (flight condition,
          thrust) that `UpdAcTecLvl!`'s inner cycle-design search holds fixed (passed by copy, never mutated)
        - `ini_eng::Vector{Float64}`: Initial guess `[BPR,pif,pilc,pihc,Tt4]` for the inner `engine_opt` search
        - `upBon_eng::Vector{Float64}`: Upper bounds for the inner search, same order as `ini_eng`
        - `loBon_eng::Vector{Float64}`: Lower bounds for the inner search, same order as `ini_eng`
        - `printEvery::Int64`: Print the current best solution every this many outer objective calls
          (also forwarded to `UpdAcTecLvl!`/`engine_opt` for their own inner-loop printouts)
        - `ftol_eng::Float64`: Forwarded to `UpdAcTecLvl!` as the inner `engine_opt` NLopt relative tolerance
        - `tol_coupling::Float64`: Forwarded to `UpdAcTecLvl!` as the innermost design/cooling_sizing
          coupling-loop tolerance (should be tighter than `ftol_eng`, which should in turn be tighter
          than whatever tolerance the caller's own outer optimizer uses, e.g. `tech_opt`'s `ftol_tec`)
        - `maxIter::Int`: Forwarded to `UpdAcTecLvl!` as the inner `engine_opt` maximum optimization iterations
        - `optTyp::Symbol`: Forwarded to `UpdAcTecLvl!` as the inner `engine_opt` NLopt algorithm symbol
        - `Fn_N::Vector{Float64}`: EEDB reference thrust levels to fly off-design (N)
        - `WFuel_kgs_ref::Vector{<:Union{Missing,Float64}}`: Reference fuel flow rate at each `Fn_N` point
          (kg/s), same length as `Fn_N`; `missing` entries skip that quantity's penalty at that point
        - `OPR_ref::Vector{<:Union{Missing,Float64}}`: Reference overall pressure ratio at each `Fn_N` point
        - `BPR_ref::Vector{<:Union{Missing,Float64}}`: Reference bypass ratio at each `Fn_N` point
        - `DFan_m_ref::Float64`: Reference fan diameter to match (m); set negative to disable this criterion
        - `M0`, `P0`, `T0`, `a0`: Inlet flight condition (Mach, Pa, K, m/s) used for every `Fn_N` off-design point

    ***Outputs***
        - `obj!`::a function: Objective function `(x, grad) -> penalty` for the outer optimizer. `x` is
          `[pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal]`. `ac` itself is never modified.
        - `histTechPara`::Vector{Vector{Float64}}: Technology parameter vectors `x` for every fully-feasible
          evaluation (sized AND every `Fn_N` point converged), returned by reference and updated by `obj!` calls
        - `histEngPara`::Vector{Vector{Float64}}: The corresponding `[BPR,pif,pilc,pihc,Tt4]` found by the
          inner loop for each fully-feasible evaluation, same order as `histTechPara`
        - `histPenl`::Vector{Float64}: Penalty for every fully-feasible evaluation, same order as `histTechPara`

    ***Behavior***
        - Penalty is ordered so that failure modes rank worse the earlier they occur: a fully-failed
          engine/tech sizing (`flgSizSuc=false`) is penalized at `1000*(length(Fn_N)+1)` -- guaranteed
          strictly worse than every `Fn_N` point independently failing to converge (`1000` each) -- so
          the optimizer is never incentivized to prefer failing sizing outright over getting partway
          through the EEDB off-design sweep
        - A candidate's evaluation continues through all of `Fn_N` even after one point fails, so the
          optimizer still sees a graded signal (fewer failed points = lower penalty) rather than a flat wall
        - History (`histTechPara`/`histEngPara`/`histPenl`) is only recorded when sizing succeeded AND every
          `Fn_N` point converged -- a candidate that "flies" everywhere but is badly mismatched against the
          reference values is NOT excluded from history by that mismatch alone, only genuine convergence failure is
        - All inputs/outputs are in SI units
"""
function make_obj_tech_cali(ac,ini_eng::Vector{Float64},upBon_eng::Vector{Float64},loBon_eng::Vector{Float64};
                            printEvery::Int64,ftol_eng::Float64=1e-7,tol_coupling::Float64=1e-8,maxIter::Int=1000,optTyp::Symbol=:LN_NELDERMEAD,
                            Fn_N::Vector{Float64},WFuel_kgs_ref::Vector{<:Union{Missing,Float64}},OPR_ref::Vector{<:Union{Missing,Float64}},BPR_ref::Vector{<:Union{Missing,Float64}},DFan_m_ref::Float64=-1.0,
                            M0::Float64=0.0,P0::Float64=101320.0,T0::Float64=288.2,a0::Float64=340.21)
    ac_used = deepcopy(ac)
    histPara_tech_cali = Vector{Vector{Float64}}()
    histPara_eng_opt =  Vector{Vector{Float64}}()
    histPenl_tech_cali = Vector{Float64}()
    num_WFuel = count(.!ismissing.(WFuel_kgs_ref))
    num_OPR = count(.!ismissing.(OPR_ref))
    num_BPR = count(.!ismissing.(BPR_ref))
    count = 0
    function obj_tech_cali!(x, grad)
        count += 1
        ac_ref = deepcopy(ac_used)
        penal=0.0
        try
            # Optimize the engine with the tech parameters
            ac_ref,flgSizSuc,bestEngSol = UpdAcTecLvl!(ac_ref,x,ini_eng,upBon_eng,loBon_eng;printEvery=printEvery,
                                             ftol_eng=ftol_eng,tol_coupling=tol_coupling,maxIter=maxIter,optTyp=optTyp)
            if flgSizSuc
                if DFan_m_ref>0.0
                    penal += 100.0*abs((ac_ref.parg[igdfan]-DFan_m_ref)/DFan_m_ref) #Percentage deviation
                end
                # Perform EEDB off-design
                flgOffDesSuc = true
                for (idxFn,Fn_N_cur) in enumerate(Fn_N)
                    try
                        res = RunEngine.runOffDes(ac_ref, M0, P0, T0, a0, Fn_N_cur)
                        penal += !res.Lconv ? 1000 : 0.0 #Represents 10 times of the deviation
                        if !ismissing(WFuel_kgs_ref[idxFn])
                            penal += 100*(abs(res.mcore*res.ff-WFuel_kgs_ref[idxFn])/WFuel_kgs_ref[idxFn])*(0.3333/num_WFuel)
                        end
                        if !ismissing(OPR_ref[idxFn])
                            penal += 100*(abs(res.OPR-OPR_ref[idxFn])/OPR_ref[idxFn])*(0.3333/num_OPR)
                        end
                        if !ismissing(BPR_ref[idxFn])
                            penal += 100*(abs(res.BPR-BPR_ref[idxFn])/BPR_ref[idxFn])*(0.3333/num_BPR)
                        end
                    catch e
                        flgOffDesSuc = false
                        if e isa InterruptException #Unless user interruption
                            rethrow()
                        end
                        penal += 1000.0 #Represents 10 times of the deviation
                    end
                end
                if flgOffDesSuc
                    # Record the feasible solution
                    push!(histPara_eng_opt,bestEngSol)
                    push!(histPara_tech_cali,x)
                    push!(histPenl_tech_cali,penal)
                end
            else
                penal += 1000.0 * (length(Fn_N) + 1) #Represents 10 times of the deviation
            end
        catch e
            (e isa InterruptException) && rethrow()
            penal += 1000.0 * (length(Fn_N) + 1) #Represents 10 times of the deviation
        end
        # Print
        if (mod(count,printEvery)==0)&&(length(histPenl_tech_cali)>0)
            idxMin = argmin(histPenl_tech_cali)
            println("TechRun#$(count): Current best for obj_tech_cali: Penalty: $(histPenl_tech_cali[idxMin]) with Tech Parameters: $(histPara_tech_cali[idxMin]) and Eng Parameters: $(histPara_eng_opt[idxMin])")
        else
            println("TechRun#$(count): No Feasible Sol Yet")
        end
        return penal
    end
    return (;obj! = obj_tech_cali!, histTechPara = histPara_tech_cali, histEngPara = histPara_eng_opt, histPenl = histPenl_tech_cali)
end

function tech_opt(ac;ini_tec::Vector{Float64},ini_eng::Vector{Float64},
                  upBon_tec::Vector{Float64},upBon_eng::Vector{Float64},
                  loBon_tec::Vector{Float64},loBon_eng::Vector{Float64},
                  Fn_N::Vector{Float64},WFuel_kgs_ref::Vector{<:Union{Missing,Float64}},OPR_ref::Vector{<:Union{Missing,Float64}},BPR_ref::Vector{<:Union{Missing,Float64}},DFan_m_ref::Float64=-1.0,
                  M0::Float64=0.0,P0::Float64=101320.0,T0::Float64=288.2,a0::Float64=340.21,
                  printEvery::Int64=10,ftol_tec::Float64=1e-6,ftol_eng::Float64=1e-7,tol_coupling::Float64=1e-8,
                  maxIter::Int=1000,optTyp::Symbol=:LN_NELDERMEAD)
    # Size check
    length(ini_tec) == 8 || error("`ini_tec` must have exactly 8 elements [pib,epolf,epollc,epolhc,epolht,epollt,etab,Tmetal], got $(length(ini_tec))")
    length(ini_eng) == 5 || error("`ini_eng` must have exactly 5 elements [BPR,pif,pilc,pihc,Tt4], got $(length(ini_eng))")
    any((ini_tec .< loBon_tec) .| (ini_tec .> upBon_tec)) && error("Initial guess `ini_tec` is outside the bounds [loBon_tec,upBon_tec]: ini_tec=$(ini_tec), loBon_tec=$(loBon_tec), upBon_tec=$(upBon_tec)")
    any((ini_eng .< loBon_eng) .| (ini_eng .> upBon_eng)) && error("Initial guess `ini_eng` is outside the bounds [loBon_eng,upBon_eng]: ini_eng=$(ini_eng), loBon_eng=$(loBon_eng), upBon_eng=$(upBon_eng)")
    ac_used = deepcopy(ac)
    # setup objective
    (;obj!,histTechPara,histEngPara,histPenl)=
    make_obj_tech_cali(ac_used,ini_eng,upBon_eng,loBon_eng;
                       printEvery=printEvery,ftol_eng=ftol_eng,tol_coupling=tol_coupling,maxIter=maxIter,optTyp=optTyp,
                       Fn_N=Fn_N,WFuel_kgs_ref=WFuel_kgs_ref,OPR_ref=OPR_ref,BPR_ref=BPR_ref,DFan_m_ref=DFan_m_ref,
                       M0=M0,P0=P0,T0=T0,a0=a0)
    # optimize
    status = :FAILURE
    try
        opt               = NLopt.Opt(optTyp, length(ini_tec))
        opt.lower_bounds  = loBon_tec
        opt.upper_bounds  = upBon_tec
        opt.min_objective = obj!
        opt.initial_step  = (upBon_tec .- loBon_tec)*0.1
        opt.ftol_rel      = ftol_tec
        opt.maxeval       = maxIter
        # Runing
        (_, _, status)    = NLopt.optimize(opt, ini_tec)
    catch e
        e isa InterruptException && rethrow()
        @error "Current optimization failed $(typeof(e)): $(e)"
    end

    #### Post-process the optimization results
    bestSol_tec = ini_tec.*NaN
    bestSol_eng = ini_eng.*NaN
    if (status in ObjectiveFactory.success_statuses) && length(histPenl)>0
        idxMin = argmin(histPenl)
        bestSol_tec = histTechPara[idxMin]
        bestSol_eng = histEngPara[idxMin]
    else
        status = :NO_FEASIBLE_SOLUTION
    end
    return  status, bestSol_tec, bestSol_eng, histTechPara, histEngPara, histPenl
end

end #CaliEng