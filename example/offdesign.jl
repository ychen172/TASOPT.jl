module PRD
export off_design_PRD

using TASOPT
include(__TASOPTindices__)

"""
off_design_PRD(ac::TASOPT.aircraft, idxFuel::Int64, rhoFuel::Float64, hvap_fuel::Float64, ranges::Vector{Float64}; epsWpay::Float64 = 1e-8)

This function determine the maximum payload range envelope with a prescribe list of possible flight range
    Assume fixed CL for off-design
    Assume starting from design point engine parameters
Inputs:
    ac: TASOPT aircraft model: has to be sized
    idxFuel: int: the fuel to be use
    rhoFuel: float: fuel density (kg/m3)
    hvap_fuel: float: Heat of vaporization of the fuel (J/kg)
    ranges: Vector{Float64}: A list of potential off-design ranges to test [m]
Outpus:
    output: Dict: ["weight_pay_list": [N] , "range_list": [m], "PFEI_list": [J/J]]
"""
function off_design_PRD(ac::TASOPT.aircraft, idxFuel::Int64, rhoFuel::Float64, hvap_fuel::Float64, ranges::Vector{Float64}; 
    epsWpay::Float64 = 1e-8)

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
    vol_fuel_max = ac.parg[igWfmax] / gee / ac.parg[igrhofuel] #Fuel volume maximum m3
    weight_empty = ac.parm[imWTO,1] - ac.parm[imWfuel,1] - ac.parm[imWpay, 1] #Empty weight (N)
    weight_one_passen = ac.parm[imWperpax, 1] #Weight of one passenger (N)
    weight_payload_max = ac.parg[igWpaymax] #Maximum weight of payload including cargo (N)

    #### Parameters to collect
    payloads_feasible = Vector{Float64}() #N
    ranges_feasible = Vector{Float64}()
    PFEI_list = Vector{Float64}()
    
    #### Initialization
    range_max_found = 1e10 #above which skip the test
    weight_payload_LB_ref = weight_one_passen #N Bound the problem
    weight_payload_UB_ref = weight_payload_max #N
    weight_payload_ini_ref = 0.5*(weight_payload_LB_ref+weight_payload_UB_ref)
    
    #### Sweep through the test range
    for range_cur in ranges
        #### Perform initial feasibility test
        # Check for feasible range
        range_cur < range_max_found || continue

        # Set the range
        ac.parm[imRange,2] = range_cur

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
            fly_mission!(ac, 2; itermax = 100, initializes_engine = true, opt_prescribed_cruise_parameter = "CL")
            weight_TO = weight_empty + ac.parm[imWfuel,2] + ac.parm[imWpay,2] #N
            vol_fuel = ac.parm[imWfuel,2] / gee / ac.parg[igrhofuel] #m3
            if weight_TO > weight_TO_max || vol_fuel > vol_fuel_max || weight_TO < 0.0 || vol_fuel < 0.0
                range_max_found = range_cur
                continue
            end
        catch
            range_max_found = range_cur
            continue
        end

        # Record the minimum result
        payload_good = weight_payload_LB_ref
        ranges_good = range_cur
        PFEI_good = ac.parm[imPFEI, 2]

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
                fly_mission!(ac, 2; itermax = 100, initializes_engine = true, opt_prescribed_cruise_parameter = "CL")
                weight_TO = weight_empty + ac.parm[imWfuel,2] + ac.parm[imWpay,2] #N
                vol_fuel = ac.parm[imWfuel,2] / gee / ac.parg[igrhofuel] #m3
                if weight_TO > weight_TO_max || vol_fuel > vol_fuel_max || weight_TO < 0.0 || vol_fuel < 0.0
                    flg_good = false
                else
                    flg_good = true
                end
            catch
                flg_good = false
            end

            # Judge
            if flg_good
                # Record
                payload_good = weight_payload_cur
                PFEI_good = ac.parm[imPFEI, 2]
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

        #### logging the current good value
        push!(payloads_feasible, payload_good)
        push!(ranges_feasible, ranges_good)
        push!(PFEI_list, PFEI_good)
    end
    #### prepare output
    output = Dict(
        "weight_pay_list" => payloads_feasible,
        "range_list" => ranges_feasible,
        "PFEI_list" => PFEI_list
    )
    return output
end

end #end module PRD