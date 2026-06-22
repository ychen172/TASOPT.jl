module PRD
export off_design_PRD, off_design_specified, off_design_R1R2

using TASOPT
include(__TASOPTindices__)

"""
off_design_PRD(ac::TASOPT.aircraft, idxFuel::Int64, rhoFuel::Float64, hvap_fuel::Float64, ranges::Vector{Float64}; 
               epsWpay::Float64 = 1e-8, epsBuff::Float64 = 0.0, save_dir::String = "ModelSaved", save_name::String = "OffDesign",  flg_save_ac::Bool = true)

This function determine the maximum payload range envelope with a prescribe list of possible flight range
    Assume fixed CL for off-design
    Assume starting from design point engine parameters
    In terms of flight model to be saved. Beware the following variable will not be same at design point anymore. (idxFuel, rhoFuel, hvap_fuel)
Inputs:
    ac: TASOPT aircraft model: has to be sized at the design mission (1st mission)
    idxFuel: int: the fuel to use
    rhoFuel: float: fuel density (kg/m3)
    hvap_fuel: float: Heat of vaporization of the fuel (J/kg)
    ranges: vector{float}: A list of potential off-design ranges to test [nmi]
    epsWpay: float: fractional search range for convergence
    epsBuff: float: small fractional buffer given to the constraint for roundoff error (Set to 0 for strictly constraint satisfying)
    save_dir: String: name of the save directory (No need if no saving)
    save_name: String: name for the saved model (save_name*string(round(Int,ran_cur))*".jld2") (No need if no saving)
    flg_save_ac: bool: true then the off-design ac models will be saved
Outpus:
    output: Dict: ["payload_weight_N": Vector{Float64} , "range_nmi": Vector{Float64}, "PFEI_JJ": Vector{Float64}, "fuel_tank_frac": Vector{Float64}, "payload_frac": Vector{Float64}] 
            If all ranges are not feasible, each element will have length 0 (Empty element)
"""
function off_design_PRD(ac::TASOPT.aircraft, idxFuel::Int64, rhoFuel::Float64, hvap_fuel::Float64, ranges::Vector{Float64}; 
    epsWpay::Float64 = 1e-8, epsBuff::Float64 = 0.0, save_dir::String = "ModelSaved", save_name::String = "OffDesign",  flg_save_ac::Bool = false)

    #### Check on sizing
    ac.is_sized[1] || error("Aircraft model needs to be sized before runing offdesign.")
    length(ranges) > 0 || error("Need to provide a list of ranges to test")

    #### Duplicate the design mission to the off-design mission for tuning
    ac = deepcopy(ac)
    parm = cat(ac.parm[:,1], ac.parm[:,1], dims=2)
    para = cat(ac.para[:,:,1], ac.para[:,:,1], dims=3)
    pare = cat(ac.pare[:,:,1], ac.pare[:,:,1], dims=3)
    ac = aircraft(ac.name, ac.description, ac.options, ac.parg, parm, para, pare, [true], 
                  ac.fuselage, ac.fuse_tank, ac.wing, ac.htail, ac.vtail, ac.engine, ac.landing_gear)
    for HX in ac.engine.heat_exchangers
        HX.HXgas_mission = cat(HX.HXgas_mission[:,1], HX.HXgas_mission[:,1], dims=2)
    end

    #### Design parameter to check the validity of the solution
    weight_TO_max = ac.parg[igWMTO] #Maximum takeoff mass (N)
    vol_fuel_max = ac.parg[igVfmax] #Fuel volume maximum m3
    weight_empty = ac.parm[imWTO,1] - ac.parm[imWfuel,1] - ac.parm[imWpay, 1] #Empty weight (N)
    weight_one_passen = ac.parm[imWperpax, 1] #Weight of one passenger (N)
    weight_payload_max = ac.parg[igWpaymax] #Maximum weight of payload including cargo (N)
    # Add some buffer to the limit
    weight_TO_max *= (1.0 + epsBuff)
    vol_fuel_max *= (1.0 + epsBuff)
    weight_payload_max *= (1.0 + epsBuff)

    #### Parameters to collect
    payloads_feasible = Vector{Float64}() #N
    ranges_feasible = Vector{Float64}() #nmi
    PFEI_list = Vector{Float64}()
    fuel_tank_frac = Vector{Float64}() #fraction of fuel tank used by volume
    payload_frac = Vector{Float64}() #fraction of payload capacity taken
    
    #### Initialization
    range_max_found = Inf #above which skip the test
    weight_payload_LB_ref = weight_one_passen #N Bound the problem
    weight_payload_UB_ref = weight_payload_max #N
    weight_payload_ini_ref = 0.5*(weight_payload_LB_ref+weight_payload_UB_ref)
    
    #### Sweep through the test range
    for range_cur in ranges
        #### Perform initial feasibility test
        # Check for feasible range
        range_cur < range_max_found || continue

        # Set the range
        ac.parm[imRange,2] = range_cur * 1852.0 # m

        # Set the fuel
        ac.options.ifuel = idxFuel
        ac.parg[igrhofuel] = rhoFuel
        ac.pare[iehvap, :, :] .= hvap_fuel #Both design and off-design mission change
        ac.pare[iehvapcombustor, :, :] .= hvap_fuel #As initialize engine try to overwrite this with design mission

        # Set the test payload
        ac.para[iaalt,ipcruise1,2] = ac.para[iaalt,ipcruise1,1] #Reset aerodynamics requirement
        ac.para[iaCL,ipcruise1,2] = ac.para[iaCL,ipcruise1,1] #Reset aerodynamics requirement
        ac.parm[imWpay,2] = weight_payload_LB_ref
        
        # Run the sizing
        try
            fly_mission!(ac, 2; itermax = 150, initializes_engine = true, opt_prescribed_cruise_parameter = "CL")
            weight_TO = weight_empty + ac.parm[imWfuel,2] + ac.parm[imWpay,2] #N
            vol_fuel = ac.parm[imVfuel, 2] #m3
            if weight_TO > weight_TO_max || vol_fuel > vol_fuel_max || weight_TO <= 0.0 || vol_fuel <= 0.0
                range_max_found = range_cur
                continue
            end
        catch err
            println(err)
            range_max_found = range_cur
            continue
        end

        # Record the minimum feasible result
        payload_good = weight_payload_LB_ref

        #### Solve for the actual maximum payload at this point
        weight_payload_LB = weight_payload_LB_ref
        weight_payload_UB = weight_payload_UB_ref
        weight_payload_cur = weight_payload_ini_ref #N
        eps_payload_cur = abs((weight_payload_UB-weight_payload_LB)/(weight_payload_UB_ref-weight_payload_LB_ref))
        count_iter = 1
        flg_good = false
        while (eps_payload_cur>epsWpay) && (count_iter<1001)
            # Set the test payload
            ac.para[iaalt,ipcruise1,2] = ac.para[iaalt,ipcruise1,1] #Reset aerodynamics requirement
            ac.para[iaCL,ipcruise1,2] = ac.para[iaCL,ipcruise1,1] #Reset aerodynamics requirement
            ac.parm[imWpay,2] = weight_payload_cur #N

            # Run the sizing
            try
                fly_mission!(ac, 2; itermax = 150, initializes_engine = true, opt_prescribed_cruise_parameter = "CL")
                weight_TO = weight_empty + ac.parm[imWfuel,2] + ac.parm[imWpay,2] #N
                vol_fuel = ac.parm[imVfuel, 2] #m3
                if weight_TO > weight_TO_max || vol_fuel > vol_fuel_max || weight_TO <= 0.0 || vol_fuel <= 0.0
                    flg_good = false
                else
                    flg_good = true
                end
            catch err
                println(err)
                flg_good = false
            end

            # Judge
            if flg_good
                # Record
                payload_good = weight_payload_cur
                # Update bound
                weight_payload_LB = weight_payload_cur
            else
                weight_payload_UB = weight_payload_cur
            end
            weight_payload_cur = 0.5*(weight_payload_LB + weight_payload_UB)
            eps_payload_cur = abs((weight_payload_UB-weight_payload_LB)/(weight_payload_UB_ref-weight_payload_LB_ref))

            # update counter
            count_iter += 1
        end

        #### Rerun the working case
        ac.para[iaalt,ipcruise1,2] = ac.para[iaalt,ipcruise1,1] #Reset aerodynamics requirement
        ac.para[iaCL,ipcruise1,2] = ac.para[iaCL,ipcruise1,1] #Reset aerodynamics requirement
        payload_good = min(payload_good, weight_payload_max)  #N
        ac.parm[imWpay,2] = payload_good
        try
            fly_mission!(ac, 2; itermax = 150, initializes_engine = true, opt_prescribed_cruise_parameter = "CL")
            weight_TO = weight_empty + ac.parm[imWfuel,2] + ac.parm[imWpay,2] #N
            vol_fuel = ac.parm[imVfuel, 2] #m3
            if weight_TO > weight_TO_max || vol_fuel > vol_fuel_max || weight_TO <= 0.0 || vol_fuel <= 0.0
                println("The final rerun did not work: Wpay: $(payload_good) N at range $(range_cur) nmi with weight_TO/max 
                $(weight_TO/weight_TO_max) and vol_fuel/max $(vol_fuel/vol_fuel_max)")
                continue
            else
                println("For off_design at range $(range_cur) nmi, find good payload $(payload_good) N with PFEI $(ac.parm[imPFEI, 2]) J/J")
                ## Logging the current good result
                push!(payloads_feasible, ac.parm[imWpay,2])
                push!(ranges_feasible, range_cur) #nmi
                push!(PFEI_list, ac.parm[imPFEI, 2])
                push!(fuel_tank_frac, vol_fuel/(vol_fuel_max/(1.0 + epsBuff)))
                push!(payload_frac, ac.parm[imWpay,2]/(weight_payload_max/(1.0 + epsBuff)))
                ## Save the current off-design model
                if flg_save_ac
                    mkpath(save_dir)
                    quicksave_aircraft(ac, joinpath(save_dir, "$(save_name*string(round(Int,range_cur))).jld2"))
                end
            end
        catch err
            println(err)
            println("The final rerun did not work: Wpay: $(payload_good) N at range $(range_cur) nmi")
            continue
        end
    end
    #### prepare output
    output = Dict(
        "payload_weight_N" => payloads_feasible,
        "range_nmi" => ranges_feasible,
        "PFEI_JJ" => PFEI_list,
        "fuel_tank_frac" => fuel_tank_frac,
        "payload_frac" => payload_frac
    )
    return output
end

"""
    findR1R2R3(mode::Symbol, R_LB::Float64, R_UB::Float64, ac, idxFuel::Int, rhoFuel::Float64, hvap_fuel::Float64; 
               epsRange::Float64 = 1e-8, epsWpay::Float64 = 1e-8, minBoundSpan::Float64 = 1.0, eps_buff_pay_vol_Range::Float64 = 1e-6, eps_buff_weight_vol_Weight::Float64 = 0.0,
               flg_save_ac::Bool = false, save_name::String = "Dummy", save_dir::String = "ModelSaved")

`findR1R2R3` finds R1, R2, R3 of a aircraft model

    **Inputs**
        - mode: Select for :R1, :R2, or :R3 to find
        - R_LB: [nmi] Lower bound for the search, needs to be at least feasible for one passenger. If R1 smaller or equal to this, then this will be the final output range too. 
        - R_UB: [nmi] upper bound. Can be large where solution where to not feasible
                R_LB and R_UB both need to be bigger than zero and R_UB needs to be minBoundSpan nmi larger than R_LB
        - ac: Aircraft model. The first mission needs to be sized. Will NOT be altered by the current function.
        - idxFuel: fuel index(type) for this mission
        - rhoFuel: (kg/m3) fuel density for this mission
        - hvap_fuel: (J/kg) fuel heat of vaporization
        - epsRange: Fractional range span WRT initial below which bisection search is complete (Outer loop for range)
        - epsWpay: Fractional payload span WRT initial below which bisection search is complete (Inner loop for weight)
        - minBoundSpan: [nmi] Minimum starting span or gap between the upper and lower range bounds
        - eps_buff_pay_vol_Range: (Outer loop) Buffer so that R1 has payload fraction > 1-this and R2 has fuel volume < 1-this. Necessary as a feasible case always have fuel volume and payload weight strictly within constraint.
        - eps_buff_weight_vol_Weight: (Inner loop) Buffer so that feasible case allow WTO fraction < 1+this and fuel volume < 1+this, NOT Necessary. 0.0 allows strictly constraint satisfying
        - flg_save_ac: To save the R mission found
        - save_name
        - save_dir: will be created if not exist

    **Outputs**
        - out_dict::Dict: (Assured not empty or else will error out) 
                    for the R mission "payload_weight_N", "range_nmi", "PFEI_JJ", "fuel_tank_frac", "payload_frac"
"""
function findR1R2R3(mode::Symbol, R_LB::Float64, R_UB::Float64, ac, idxFuel::Int, rhoFuel::Float64, hvap_fuel::Float64;
                    epsRange::Float64 = 1e-8, epsWpay::Float64 = 1e-8, minBoundSpan::Float64 = 1.0, eps_buff_pay_vol_Range::Float64 = 1e-6, eps_buff_weight_vol_Weight::Float64 = 0.0,
                    flg_save_ac::Bool = false, save_name::String = "Dummy", save_dir::String = "ModelSaved")
    #### Size check
    ((R_LB > 0.0) && (R_UB > 0.0)) || throw(ArgumentError("Upper $(R_UB) and lower $(R_LB) range bounds have to have larger than 0.0"))
    ((R_UB - R_LB)>minBoundSpan) || throw(ArgumentError("R_UB - R_LB $(R_UB - R_LB) has to be bigger than $(minBoundSpan) nmi"))
    test_out = off_design_PRD(ac, idxFuel, rhoFuel, hvap_fuel, [R_LB]; flg_save_ac = false, epsWpay = epsWpay, epsBuff = eps_buff_weight_vol_Weight)
    (length(test_out["payload_frac"]) > 0) || throw(ArgumentError("Input lower range bound $(R_LB) does not lead to any feasible solution"))

    #### Initialization
    gap_range_ref = R_UB-R_LB #Reference search range span for convergence calculation
    frac_range_cur = 1.0 #Current fractional search range
    R_guess = 0.5*(R_UB+R_LB) #Initial guess [nmi]
    R_out = R_LB #final solution
    count_iter = 0
    while  (frac_range_cur>epsRange) && (count_iter<1000)
        # Test the new guess range
        out_dict = off_design_PRD(ac, idxFuel, rhoFuel, hvap_fuel, [R_guess]; flg_save_ac = false, epsWpay = epsWpay, epsBuff = eps_buff_weight_vol_Weight)

        # Judge
        cond_upper_bound = nothing
        if mode == :R1
            #If case simulation failed (flight range too far to be feasible) or payload is not at max capacity, that is an upper bound. If R1<MinRange, UB will be as short as possible.
            cond_upper_bound = (length(out_dict["payload_frac"]) <= 0) || (out_dict["payload_frac"][1] <= 1.0-eps_buff_pay_vol_Range)
        elseif mode == :R2
            #If case simulation failed (flight range too far to be feasible) or fuel tank is at max capacity, that is an upper bound.
            cond_upper_bound = (length(out_dict["payload_frac"]) <= 0) || (out_dict["fuel_tank_frac"][1] >= 1.0-eps_buff_pay_vol_Range)
        elseif mode == :R3
            #If case simulation failed (flight range too far to be feasible), that is an upper bound.
            cond_upper_bound = (length(out_dict["payload_frac"]) <= 0)
        else
            error("Unknown Mode of Calculation $(mode)")
        end
        
        # Set new guess
        if cond_upper_bound
            R_UB = R_guess
        else
            R_LB = R_guess
            R_out = R_guess
        end
        
        # Recompute convergence criteria
        frac_range_cur = (R_UB-R_LB)/gap_range_ref
        R_guess = 0.5*(R_UB+R_LB)
        count_iter += 1
    end
    (count_iter < 1000) || @warn "Maximum iterations reached for range search loop(Unlikely due to bisection nature(something is off))"
    
    #### Rerun the feasible case and output
    out_dict = off_design_PRD(ac, idxFuel, rhoFuel, hvap_fuel, [R_out]; 
                              epsWpay = epsWpay, epsBuff = eps_buff_weight_vol_Weight,
                              flg_save_ac = flg_save_ac, save_name = save_name, save_dir = save_dir)
    (length(out_dict["payload_frac"]) > 0) || throw(ErrorException("Determined range $(R_out) does not lead to any feasible solution"))
    return out_dict
end

"""
off_design_specified(ac::TASOPT.aircraft, idxFuel::Int64, rhoFuel::Float64, hvap_fuel::Float64, ranges::Vector{Float64}, weights_payload::Vector{Float64};
                     save_dir::String = "ModelSaved", save_name::String = "OffDesign")

This function take a design-point model and run a series of off-design missions
    Specify one fuel type for off-design
    Specify a series of payload and range pair (No saving if mission failed)
    Assume fixed CL for off-design
    Assume starting from design point engine parameters
Inputs:
    ac: TASOPT aircraft model: has to be sized
    idxFuel: int: the fuel to be use
    rhoFuel: float: fuel density (kg/m3)
    hvap_fuel: float: Heat of vaporization of the fuel (J/kg)
    ranges: vector{float}: A list of off-design ranges to test [nmi]
    weights_payload: Vector{Float64}: A list of off-design payload to test [N]
    save_dir: String: name of the save directory
    save_name: String: name for the saved model (save_name*string(round(Int,ran_cur))*".jld2")
Outpus:
    output: Dict: ["payload_weight_N": Vector{Float64} , "range_nmi": Vector{Float64}, "PFEI_JJ": Vector{Float64}]
"""
function off_design_specified(ac::TASOPT.aircraft, idxFuel::Int64, rhoFuel::Float64, hvap_fuel::Float64, ranges::Vector{Float64}, weights_payload::Vector{Float64};
                        save_dir::String = "ModelSaved", save_name::String = "OffDesign")

    #### Check on sizing
    ac.is_sized[1] || error("Aircraft model needs to be sized before runing offdesign.")
    length(ranges) > 0 || error("Need to provide a list of ranges to test")
    length(ranges) == length(weights_payload) || error("Need to provide pairs of ranges and payloads")

    #### Duplicate the design mission to the off-design mission for tuning
    ac = deepcopy(ac)
    parm = cat(ac.parm[:,1], ac.parm[:,1], dims=2)
    para = cat(ac.para[:,:,1], ac.para[:,:,1], dims=3)
    pare = cat(ac.pare[:,:,1], ac.pare[:,:,1], dims=3)
    ac = aircraft(ac.name, ac.description, ac.options, ac.parg, parm, para, pare, [true], 
                  ac.fuselage, ac.fuse_tank, ac.wing, ac.htail, ac.vtail, ac.engine, ac.landing_gear)
    for HX in ac.engine.heat_exchangers
        HX.HXgas_mission = cat(HX.HXgas_mission[:,1], HX.HXgas_mission[:,1], dims=2)
    end

    #### Design parameter to check the validity of the solution
    weight_TO_max = ac.parg[igWMTO] #Maximum takeoff mass (N)
    vol_fuel_max = ac.parg[igWfmax] / gee / ac.parg[igrhofuel] #Fuel volume maximum (m3)
    weight_empty = ac.parm[imWTO,1] - ac.parm[imWfuel,1] - ac.parm[imWpay, 1] #Empty weight (N)
    weight_payload_max = ac.parg[igWpaymax] #Maximum weight of payload including cargo (N). A +1 N margin will be given for numerical error.
    # Add some buffer to the limit
    weight_TO_max *= 1.0001 
    vol_fuel_max *= 1.0001
    weight_payload_max *= 1.0001

    #### Parameters to collect
    payloads_feasible = Vector{Float64}() # (N)
    ranges_feasible = Vector{Float64}() # (nmi)
    PFEI_list = Vector{Float64}() #(J/J)
        
    #### Sweep through the test range
    for (i, range_cur) in enumerate(ranges)
        #### Setup the case to test
        # Set the range
        ac.parm[imRange,2] = range_cur * 1852.0 # m
        
        # Set the payload
        ac.parm[imWpay,2] = weights_payload[i]

        # Set the fuel
        ac.options.ifuel = idxFuel
        ac.parg[igrhofuel] = rhoFuel
        ac.pare[iehvap, :, :] .= hvap_fuel #Both design and off-design mission change
        ac.pare[iehvapcombustor, :, :] .= hvap_fuel #As initialize engine try to overwrite this with design mission

        # Reset flight requirement
        ac.para[iaalt,ipcruise1,2] = ac.para[iaalt,ipcruise1,1] #Reset aerodynamics requirement
        ac.para[iaCL,ipcruise1,2] = ac.para[iaCL,ipcruise1,1] #Reset aerodynamics requirement
        
        #### Attempt the off-design mission
        try
            fly_mission!(ac, 2; itermax = 100, initializes_engine = true, opt_prescribed_cruise_parameter = "CL")
            weight_TO = weight_empty + ac.parm[imWfuel,2] + ac.parm[imWpay,2] #N
            vol_fuel = ac.parm[imWfuel,2] / gee / ac.parg[igrhofuel] #m3
            if weight_TO > weight_TO_max || vol_fuel > vol_fuel_max || ac.parm[imWpay,2] > weight_payload_max ||
               weight_TO < 0.0 || vol_fuel < 0.0 || ac.parm[imWpay,2] < 0.0
                println("Off-design failed constraints: Wpay: $(ac.parm[imWpay,2]) N at range $(range_cur) nmi 
                with WTO/max $(weight_TO/weight_TO_max); Vf/max $(vol_fuel/vol_fuel_max); Wpay/max $(ac.parm[imWpay,2]/weight_payload_max)")
                continue
            else
                println("Off-design succeeded: Wpay: $(ac.parm[imWpay,2]) N at range $(range_cur) nmi")
                ## Logging the result
                push!(payloads_feasible, ac.parm[imWpay,2]) #(N)
                push!(ranges_feasible, range_cur) #(nmi)
                push!(PFEI_list, ac.parm[imPFEI, 2]) #(J/J)

                ## Save the current off-design model
                mkpath(save_dir)
                quicksave_aircraft(ac, joinpath(save_dir, "$(save_name*string(round(Int,range_cur))).jld2"))
            end
        catch err
            println(err)
            println("Off-design failed to run: Wpay: $(ac.parm[imWpay,2]) N at range $(range_cur) nmi")
            continue
        end
    end

    #### prepare output
    output = Dict(
        "payload_weight_N" => payloads_feasible,
        "range_nmi" => ranges_feasible,
        "PFEI_JJ" => PFEI_list
    )

    return output
end

end #end module PRD