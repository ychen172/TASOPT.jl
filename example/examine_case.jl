using CSV, DataFrames
using TASOPT, NLopt
include(joinpath(@__DIR__, "optimize_rangefuel.jl"))
using .OptimizeRangeFuel: BoundsOpt, ConstraintsOpt

"""
example_design!(ac, save_dir, save_name; bounds_opt, constraints_opt)

Post-process the design point information by outputing the parameters into csv and return design parameters too
Inputs:
    ac: TASOPT.aircraft: Sized aircraft model for post-processings
    save_dir: string: Directory for saving model
    save_name: string: Name of to save the csv file(Keywork only, exclude the csv part)
    flg_save: Bool: if true extracted design values will be saved
    bounds_opt: BoundsOpt: the boundaries for optimization, 15 parameters, [AR,CL,sweep(deg),altitude,λ_in,λ_out,t/c_root,t/c_span,rcls,rclt,Tt4,π_hc,π_f,π_lc,BPR]
    constraints_opt: ConstraintsOpt: the constraints for optimization, 6 parameters, [max_span,max_lenField,min_TOCGamma,max_Tt3,max_TMetal,max_DiaFan]
Saves(Activated if flg_save):
    DesignParameters.csv: Key design parameters (Always)
    StickPlot.png: Aicraft size plot (Always)
    LimitingParameters.csv: Constraints plots (if constraints_opt provided)
    BoundaryParameters.csv: Design parameters plots (if bounds_opt provided)
Outputs:
    designParam: name tuple: extracted design parameters
"""
function example_design!(ac::TASOPT.aircraft, save_dir::AbstractString, save_name::AbstractString;
    flg_save::Bool=true, bounds_opt::Union{Nothing,BoundsOpt}=nothing, constraints_opt::Union{Nothing,ConstraintsOpt}=nothing)
    
    #### Setup a saving directory 
    out_dir = joinpath(save_dir, save_name)
    mkpath(out_dir)
    ac.is_sized[1] || throw(ArgumentError("Input aircraft model has to be sized"))

    #### Extract design parameters
    # Mission Parameters
    massPayload = ac.parm[imWpay, 1]/gee/1000.0 #(Ton)
    numberPassengers = ac.parm[imWpay, 1]/ac.parm[imWperpax, 1] #Number
    range = ac.parm[imRange, 1]/1852.0 #nmi
    PFEI = ac.parm[imPFEI, 1] #(J/J)
    # Input Optimized Parameter
    AR = ac.wing.layout.AR #Wing aspect ratio
    CL_cruise = ac.para[iaCL, ipcruise1, 1] #Cruise CL
    sweep = ac.wing.layout.sweep #Wing sweep angle (deg)
    alt_cruise = ac.para[iaalt, ipcruise1, 1]* 3.280839895 #Cruise altitude (ft)
    taper_wing_in = ac.wing.inboard.λ
    taper_wing_out = ac.wing.outboard.λ
    thick_to_chord_in = ac.wing.inboard.cross_section.thickness_to_chord
    thick_to_chord_out = ac.wing.outboard.cross_section.thickness_to_chord
    break_root_cl_ratio_cruise = ac.para[iarcls, ipcruise1, 1]
    tip_root_cl_ratio_cruise = ac.para[iarclt, ipcruise1, 1]
    Tt4_cruise = ac.pare[ieTt4, ipcruise1, 1] #(K)
    PR_hpc_cruise = ac.pare[iepihc, ipcruise1, 1] #HPC pressure ratio
    PR_fan_cruise = ac.pare[iepif, ipcruise1, 1] #Fan
    PR_lpc_cruise = ac.pare[iepilc, ipcruise1, 1] #LPC
    BPR_cruise = ac.pare[ieBPR, ipcruise1, 1] #bypass ratio at cruise
    # Parameters to Flight Efficiencies
    LD_cruise = 0.5 * (ac.para[iaCL, ipcruise1, 1]/ac.para[iaCD, ipcruise1, 1] + 
                            ac.para[iaCL, ipcruise2, 1]/ac.para[iaCD, ipcruise2, 1]) #Averaged cruise lift-to-drag ratio
    LHV_cruise = 0.5 * (ac.pare[iehfuel, ipcruise1, 1] + ac.pare[iehfuel, ipcruise2, 1]) #Averaged cruise heating value (J/kg) (Include vaporization already)
    TSFC_cruise = 0.5 * (ac.pare[ieTSFC, ipcruise1, 1] + ac.pare[ieTSFC, ipcruise2,1]) / gee #Averaged cruise thrust specfic heat consumption (kg/s/N)
    TSEC_cruise = TSFC_cruise*LHV_cruise #Averaged cruise thrust specific energy consumption (J/s/N)
    vel_cruise = 0.5 * (cos(ac.para[iagamV, ipcruise1,1]) * ac.pare[ieu0, ipcruise1,1] + 
                        cos(ac.para[iagamV, ipcruise2,1]) * ac.pare[ieu0, ipcruise2, 1]) #Averaged cruise horizontal velocity (m/s)
    massTO = ac.parm[imWTO,1]/gee/1000.0 #Takeoff mass (Ton)
    massFuelTot = ac.parm[imWfuel,1]/gee/1000.0 #Fuel mass (Ton)
    fracFuelReserved_ = ac.parg[igfreserve]/(ac.parg[igfreserve]+1.0) #fraction reserved fuel (mFuelRes/mFuelTotal)
    massFuelReserved = massFuelTot*fracFuelReserved_ # Fuel mass reserved (Ton)
    massFuelBurned = massFuelTot*(1.0 - fracFuelReserved_) #Fuel mass burned (Ton)
    massEmpty = massTO - massFuelBurned - massFuelReserved - massPayload #By definition, the empty weight (Ton)
    weightRatio_ = massTO/(massTO-massFuelBurned) #Initial weight / Final weight
    rangeBreguet = ((vel_cruise * LD_cruise)/(gee * TSFC_cruise)) * log(weightRatio_) / 1852.0 #Estimated range using Breguet range equation
    # Parameters to engine and combustor performance
    numberEngines = ac.parg[igneng]
    thrustOneEngine_takeoff = ac.pare[ieFe, iptakeoff, 1] / 1000.0 #Takeoff thrust [kN] 
    thrustOneEngine_climb1 = ac.pare[ieFe, ipclimb1, 1] / 1000.0 #Takeoff thrust [kN]
    thrustOneEngine_cruise1 = ac.pare[ieFe, ipcruise1, 1] / 1000.0 #Cruise thrust [kN]
    # Constrained parameters
        # (massTO) (massFuelTot)
    diaFan = ac.parg[igdfan] #Fan diameter (m)
    TMetalMax = maximum(ac.pare[ieTmet1, :, 1]) #Maximum metal temperature (K)
    Tt3Max = maximum(ac.pare[ieTt3, :, 1]) #Maximum compressor outlet temperature (K)
    gamTOC = ac.para[iagamV, ipclimbn, 1]*180/pi #Top of climb flight angle (deg)
    lenFieldBalanced = ac.parm[imlBF, 1] #balanced field length (m)
    spanWing = ac.wing.span #Wing span (m)

    # Create a named tuple for these design parameters
    designParam = (;
    massPayload,
    numberPassengers,
    range,
    PFEI,
    AR,
    CL_cruise,
    sweep,
    alt_cruise,
    taper_wing_in,
    taper_wing_out,
    thick_to_chord_in,
    thick_to_chord_out,
    break_root_cl_ratio_cruise,
    tip_root_cl_ratio_cruise,
    Tt4_cruise,
    PR_hpc_cruise,
    PR_fan_cruise,
    PR_lpc_cruise,
    BPR_cruise,
    LD_cruise,
    LHV_cruise,
    TSFC_cruise,
    TSEC_cruise,
    vel_cruise,
    massTO,
    massFuelTot,
    massFuelReserved,
    massFuelBurned,
    massEmpty,
    rangeBreguet,
    numberEngines,
    thrustOneEngine_takeoff,
    thrustOneEngine_climb1,
    thrustOneEngine_cruise1,
    diaFan,
    TMetalMax,
    Tt3Max,
    gamTOC,
    lenFieldBalanced,
    spanWing
    )
    # Quick exit if not saving anything
    !flg_save && return designParam

    #### Output the design parameters
    desVar_names = [
    "Payload mass (Ton)",
    "Number of passengers",
    "Flight range (nmi)",
    "PFEI (J/J)",
    "Wing aspect ratio",
    "Cruise CL",
    "Wing sweep angle (deg)",
    "Cruise altitude (ft)",
    "Inboard wing taper ratio",
    "Outboard wing taper ratio",
    "Inboard thickness-to-chord ratio",
    "Outboard thickness-to-chord ratio",
    "Break/root CL ratio at cruise",
    "Tip/root CL ratio at cruise",
    "Tt4 at cruise (K)",
    "HPC pressure ratio at cruise",
    "Fan pressure ratio at cruise",
    "LPC pressure ratio at cruise",
    "Bypass ratio at cruise",
    "Averaged cruise lift-to-drag ratio",
    "Averaged cruise heating value (J/kg)",
    "Averaged cruise TSFC (kg/s/N)",
    "Averaged cruise TSEC (J/s/N)",
    "Averaged cruise horizontal velocity (m/s)",
    "Takeoff mass (Ton)",
    "Total fuel mass (Ton)",
    "Reserved fuel mass (Ton)",
    "Burned fuel mass (Ton)",
    "Empty mass (Ton)",
    "Breguet flight range (nmi)",
    "Number of engines",
    "Thrust at takeoff (one engine) (kN)",
    "Thrust at climb1 (one engine) (kN)",
    "Thrust at cruise1 (one engine) (kN)",
    "Fan diameter (m)",
    "Maximum metal temperature (K)",
    "Maximum compressor exit temperature Tt3 (K)",
    "Top of climb flight angle (deg)",
    "Balanced field length (m)",
    "Wing span (m)"
    ]

    desVar_values = [
    massPayload,
    numberPassengers,
    range,
    PFEI,
    AR,
    CL_cruise,
    sweep,
    alt_cruise,
    taper_wing_in,
    taper_wing_out,
    thick_to_chord_in,
    thick_to_chord_out,
    break_root_cl_ratio_cruise,
    tip_root_cl_ratio_cruise,
    Tt4_cruise,
    PR_hpc_cruise,
    PR_fan_cruise,
    PR_lpc_cruise,
    BPR_cruise,
    LD_cruise,
    LHV_cruise,
    TSFC_cruise,
    TSEC_cruise,
    vel_cruise,
    massTO,
    massFuelTot,
    massFuelReserved,
    massFuelBurned,
    massEmpty,
    rangeBreguet,
    numberEngines,
    thrustOneEngine_takeoff,
    thrustOneEngine_climb1,
    thrustOneEngine_cruise1,
    diaFan,
    TMetalMax,
    Tt3Max,
    gamTOC,
    lenFieldBalanced,
    spanWing
    ]
    df_desVar = DataFrame(
    Symbol("Variables") => desVar_names,
    Symbol("Values") => desVar_values
    )
    # write out the parameters
    CSV.write(joinpath(save_dir,"DesPara_$(save_name).csv"), df_desVar)

    #### Output the aircraft stick plots
    pic = TASOPT.stickfig(ac)
    savefig(pic, joinpath(save_dir,"StickPlot_$(save_name).png"))

    #### Output constraints
    

    return designParam
end