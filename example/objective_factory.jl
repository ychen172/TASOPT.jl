module ObjectiveFactory
export OptHistory, make_obj, best_feasible
using TASOPT
using Printf
include(__TASOPTindices__)
"""
make_obj() creates an objective function that store an aircraft model with setup mission information and store
constraints for the output aircraft parameters. The objective function can then be tested with different design variables.
"""
struct Constraint
    name::String
    value::Float64
    limit::Float64
    penalty::Float64
end

mutable struct OptHistory
    test_param::Vector{Vector{Float64}}
    penalty::Vector{Float64}
    PFEI::Vector{Float64}
    violations::Vector{Vector{Constraint}}
end

OptHistory() = OptHistory(
    Vector{Vector{Float64}}(),
    Vector{Float64}(),
    Vector{Float64}(),
    Vector{Vector{Constraint}}()
)

function best_feasible(hist_all)
    idx = findall(v -> isempty(v), hist_all.violations)
    isempty(idx) && return nothing
    i = idx[argmin(hist_all.penalty[idx])]
    return (index=i, test_param=hist_all.test_param[i], penalty=hist_all.penalty[i], PFEI=hist_all.PFEI[i])
end

function make_obj(ac::TASOPT.aircraft, constraints::Vector{Float64}, hist::OptHistory, 
                  penal_scale::Vector{Float64}=[25.0, 1.0, 1.0, 5.0, 5.0, 1.0]; print_every::Int=10, iter_size_loop::Int=150)
    """
    All inputs in SI units (Exceptions: sweep(Deg), )
    ac: Aircraft model with the modified mission parameters (Ex. target range, and target payload, or fuel type)
    constraints: [max_span, max_lenField, min_TOCGamma, max_Tt3, max_TMetal, max_DiaFan] #There are also other sanity checks will be automatically included
                    1          2             3            4          5           6
    penal_scale: [max_span, max_lenField, min_TOCGamma, max_Tt3, max_TMetal, max_DiaFan] #scalings for the constraints
    hist: Store the optimization history

    print_every: iterations
    iter_size_loop: for the sizing loop in each iteration
    """
    ####Size check
    @assert length(constraints) == 6
    @assert length(penal_scale) == length(constraints)
    @assert print_every >= 1
    @assert iter_size_loop >= 1

    ####Construct objective function
    function obj(x, grad)
        """
        x: [AR, CL, sweep, altitude, λ_in, λ_out, t/c_root, t/c_span, rcls, rclt, Tt4, π_hc, π_f, π_lc, BPR]
             1.  2.   3.      4.      5.     6.      7.         8.      9.   10.   11.   12.  13.  14.  15.
        """
        #### Check if gradients provided:
        if !isempty(grad)
            # for derivative-free algorithms this is usually empty anyway
            # for gradient-based ones, you must compute and assign grad[:] properly
            fill!(grad, 0.0)
        end

        ####Size check
        @assert length(x) == 15

        ####Overwrite the design parameters into the aircraft model
        # Update wing parameters
        ac.wing.layout.AR    = x[1]                                # Aspect Ratio 
        ac.wing.layout.sweep = x[3]                                # Sweep angle (deg)
        ac.wing.inboard.λ    = x[5]                                # Inner panel taper ratio
        ac.wing.outboard.λ   = x[6]                                # Outer panel taper ratio
        ac.wing.inboard.cross_section.thickness_to_chord  = x[7]   # Root thickness-to-chord ratio
        ac.wing.outboard.cross_section.thickness_to_chord = x[8]   # Spanbreak thickness-to-chord ratio

        # Update flight condition parameters
        ac.para[iaCL, ipclimb1+1:ipdescentn-1, 1] .= x[2]       # Cruise lift coefficient
        ac.para[iaalt, ipclimbn:ipcruise1, 1]     .= x[4]       # Cruise altitude
        ac.para[iarcls, ipclimb2:ipdescent4, 1]   .= x[9]       # Break/root CL ratio = cls/clo
        ac.para[iarclt, ipclimb2:ipdescent4, 1]   .= x[10]      # Tip/root CL ratio = clt/clo

        # Update engine parameters
        ac.pare[ieTt4, ipcruise1:ipcruise2, 1] .= x[11]         # Turbine inlet temperature
        ac.pare[iepihc, ipcruise1, 1]           = x[12]         # High pressure compressor pressure ratio
        ac.pare[iepif, ipcruise1, 1]            = x[13]         # Fan pressure ratio
        ac.pare[iepilc, ipcruise1, 1]           = x[14]         # Low pressure compressor pressure ratio (fixed)
        ac.pare[ieBPR, ipcruise1, 1]            = x[15]         # Bypass ratio

        #### Test size the current case
        flgSuccessed = true
        try
            TASOPT.size_aircraft!(ac, iter=iter_size_loop, printiter=false)
        catch e
            if e isa InterruptException
                rethrow()
            end
            flgSuccessed = false
        end

        ####Collect the primary cost function
        if flgSuccessed
            main_penalty = ac.parm[imPFEI, 1] #J/J
        else
            main_penalty = 1e10
        end
        
        ####Additional penalty because of the violation of the specified constraints
        add_penalty = 0.0
        violated_constraints = Constraint[] #Vector of constraints

        # Maximum span (Max) #1
        max_span = constraints[1]
        cur_span = ac.wing.span
        if cur_span > max_span
            err_frac = cur_span/max_span - 1.0 #fractional error
            penalty = penal_scale[1] * ac.parg[igWpay] * err_frac^2
            add_penalty += penalty
            push!(violated_constraints, Constraint("Wing Span", cur_span, max_span, penalty))
        end

        # Balanced field length (Max) #2
        max_lenField = constraints[2]
        cur_lenField = ac.parm[imlBF,1]
        if cur_lenField > max_lenField
            err_frac = cur_lenField/max_lenField - 1.0
            penalty = penal_scale[2] * ac.parg[igWpay] * err_frac^2
            add_penalty += penalty
            push!(violated_constraints, Constraint("Balanced Field Length", cur_lenField, max_lenField, penalty))
        end

        # Top of Climb Flight Angle (Min) #3
        min_TOCGamma = constraints[3]
        cur_TOCGamma = ac.para[iagamV, ipclimbn, 1]
        if cur_TOCGamma < min_TOCGamma
            err_frac = 1.0 - cur_TOCGamma/min_TOCGamma
            penalty = penal_scale[3] * ac.parg[igWpay] * err_frac^2
            add_penalty += penalty
            push!(violated_constraints, Constraint("Top of Climb Flight Angle", cur_TOCGamma, min_TOCGamma, penalty))
        end

        # Combustor Inlet Temperature (Max) #4
        max_Tt3 = constraints[4]
        cur_Tt3 = maximum(ac.pare[ieTt3, :, 1])
        if cur_Tt3 > max_Tt3
            err_frac = cur_Tt3/max_Tt3 - 1.0
            penalty = penal_scale[4] * ac.parg[igWpay] * err_frac^2
            add_penalty += penalty
            push!(violated_constraints, Constraint("Combustor Inlet Temperature", cur_Tt3, max_Tt3, penalty))
        end

        # Metal Temperature (Max) #5
        max_TMetal = constraints[5]
        cur_TMetal = maximum(ac.pare[ieTmet1, :, 1])
        if cur_TMetal > max_TMetal
            err_frac = cur_TMetal/max_TMetal - 1.0
            penalty = penal_scale[5] * ac.parg[igWpay] * err_frac^2
            add_penalty += penalty
            push!(violated_constraints, Constraint("Metal Temperature", cur_TMetal, max_TMetal, penalty))
        end

        # Fan Diameter (Max) #6
        max_DiaFan = constraints[6]
        cur_DiaFan = ac.parg[igdfan]
        if cur_DiaFan > max_DiaFan
            err_frac = cur_DiaFan/max_DiaFan - 1.0
            penalty = penal_scale[6] * ac.parg[igWpay] * err_frac^2
            add_penalty += penalty
            push!(violated_constraints, Constraint("Fan Diameter", cur_DiaFan, max_DiaFan, penalty))
        end

        ####Additional penalty from sanity check

        # Takeoff Weight (Max) #1
        max_WTO = ac.parg[igWMTO]
        cur_WTO = ac.parm[imWTO, 1]
        penScale_WTO = 10.0
        if cur_WTO > max_WTO+1.0
            err_frac = cur_WTO/max_WTO - 1.0
            penalty = penScale_WTO * ac.parg[igWpay] * err_frac^2
            add_penalty += penalty
            push!(violated_constraints, Constraint("Takeoff Weight (San)", cur_WTO, max_WTO, penalty))
        end

        # Fuel Volume/Weight (Max) #2
        max_WFuel = ac.parg[igWfmax] #by volume limitation
        cur_WFuel = ac.parg[igWfuel] #by energy limitation
        penScale_WFuel = 10.0
        if cur_WFuel > max_WFuel
            err_frac = cur_WFuel/max_WFuel - 1.0
            penalty = penScale_WFuel * ac.parg[igWpay] * err_frac^2
            add_penalty += penalty
            push!(violated_constraints, Constraint("Fuel Weight (San)", cur_WFuel, max_WFuel, penalty))
        end

        ####Final Penalty
        total_penalty = main_penalty + add_penalty

        ####Update the Optimization History
        push!(hist.test_param, copy(x))
        push!(hist.penalty, total_penalty)
        push!(hist.PFEI, main_penalty)
        push!(hist.violations, copy(violated_constraints))

        ####In-flight Printout
        cur_iter_optim = length(hist.penalty)
        if cur_iter_optim == 1 || cur_iter_optim % print_every == 0
            println()
            @printf("Iter %4d | total_penalty = %.4f | PFEI = %.4f\n", cur_iter_optim, total_penalty, main_penalty)
            println()
            @printf("Violates: ")
            for cur_violates in violated_constraints
                @printf("| %10s", cur_violates.name*"⚠️")
            end
            println()
        end
        return total_penalty
    end
    return obj
end

end # module ObjectiveFactory