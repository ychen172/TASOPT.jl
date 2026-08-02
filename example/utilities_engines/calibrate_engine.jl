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


end #CaliEng