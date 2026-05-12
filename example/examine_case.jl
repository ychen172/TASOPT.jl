# This is an example file to examine an aircraft model

# 1) Load TASOPT
using TASOPT
using DataFrames, CSV
using Plots
include(__TASOPTindices__)

# 2) Include input file for desired aircraft/
#  load the target model
load_dir = "ModelSaved"
load_name = "acOptimized_Cus2"
ac = quickload_aircraft(joinpath(load_dir,"$(load_name).jld2")) # simply a synonym to read_aircraft_model()
save_dir = "ModelProcessed"
save_dir = joinpath(save_dir,load_name)
mkpath(save_dir)

# 3) Extract constrained parameters
wingSpan_Cur = ac.wing.span #(m)
wingSpan_Max = ac.wing.layout.max_span #(m)
println("Wing Span (m): Cur: $wingSpan_Cur, Max: $wingSpan_Max")

fieldLengthBalanced_Cur = ac.parm[imlBF, 1] #(m)
fieldLengthBalanced_Max = 2.4e3 #(m)
println("Balanced Field Length (m): Cur: $fieldLengthBalanced_Cur, Max: $fieldLengthBalanced_Max")

TOCFlightAngle_Cur = ac.para[iagamV, ipclimbn, 1]*180/pi #(deg)
TOCFlightAngle_Min = ac.parg[iggtocmin]*180/pi #(deg)
println("Top of Climb Flight Angle (deg): Cur: $TOCFlightAngle_Cur, Min: $TOCFlightAngle_Min")

Tt3_Cur = maximum(ac.pare[ieTt3, :, 1]) #(K)
Tt3_Max = 900.0
println("Compressor Exit Temperature (K): Cur: $Tt3_Cur, Max: $Tt3_Max")

TMetal_Cur = maximum(ac.pare[ieTmet1, :, 1]) #(K)
TMetal_Max = 1333.33 #(K)
println("Metal Temperature (K): Cur: $TMetal_Cur, Max: $TMetal_Max")

DiaFan_Cur = ac.parg[igdfan] #(m)
DiaFan_Max = 2.0 #(m)
println("Fan Diameter (m): Cur: $DiaFan_Cur, Max: $DiaFan_Max")

FuelMass_Cur = ac.parg[igWfuel]/9.81/1000. #(Ton)
FuelMass_Max = ac.parg[igWfmax]/9.81/1000. #(Ton)
println("Fuel Mass (Ton): Cur: $FuelMass_Cur, Max: $FuelMass_Max")

TakeoffMass_Cur = ac.parm[imWTO, 1]/9.81/1000. #(Ton)
TakeoffMass_Max = ac.parg[igWMTO]/9.81/1000. #(Ton)
println("Takeoff Mass (Ton): Cur: $TakeoffMass_Cur, Max: $TakeoffMass_Max")

# 3) Collect and save design parameters
#Mission Parameters
massPayload = ac.parm[imWpay, 1]/gee/1000.0 #(Ton)
numberPassengers = ac.parm[imWpay, 1]/ac.parm[imWperpax, 1] #Number
range = ac.parm[imRange, 1]/1852.0 #nmi
PFEI = ac.parm[imPFEI, 1] #(J/J)
#Input Optimized Parameter
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
#Parameters to Flight Efficiencies
LD_cruise = 0.5 * (ac.para[iaCL, ipcruise1, 1]/ac.para[iaCD, ipcruise1, 1] + 
                        ac.para[iaCL, ipcruise2, 1]/ac.para[iaCD, ipcruise2, 1]) #Averaged cruise lift-to-drag ratio
LHV_cruise = 0.5 * (ac.pare[iehfuel, ipcruise1, 1] + ac.pare[iehfuel, ipcruise2, 1]) #Averaged cruise heating value (J/kg) (Include vaporization already)
TSFC_cruise = 0.5 * (ac.pare[ieTSFC, ipcruise1, 1] + ac.pare[ieTSFC, ipcruise2,1]) / gee #Averaged cruise thrust specfic heat consumption (kg/s/N)
TSEC_cruise = TSFC_cruise*LHV_cruise #Averaged cruise thrust specific energy consumption (J/s/N)
vel_cruise = 0.5 * (cos(ac.para[iagamV, ipcruise1,1]) * ac.pare[ieu0, ipcruise1,1] + 
                    cos(ac.para[iagamV, ipcruise2,1]) * ac.pare[ieu0, ipcruise2, 1]) #Averaged cruise horizontal velocity (m/s)
massTO = ac.parm[imWTO,1]/gee/1000.0 #Takeoff mass (Ton)
massFuel_ = ac.parm[imWfuel,1]/gee/1000.0 #Fuel mass (Ton)
fracFuelReserved_ = ac.parg[igfreserve]/(ac.parg[igfreserve]+1.0) #fraction reserved fuel (mFuelRes/mFuelTotal)
massFuelReserved = massFuel_*fracFuelReserved_ # Fuel mass reserved (Ton)
massFuelBurned = massFuel_*(1.0 - fracFuelReserved_) #Fuel mass burned (Ton)
massEmpty = massTO - massFuelBurned - massFuelReserved - massPayload #By definition, the empty weight (Ton)
weightRatio_ = massTO/(massTO-massFuelBurned) #Initial weight / Final weight
rangeBreguet = ((vel_cruise * LD_cruise)/(gee * TSFC_cruise)) * log(weightRatio_) / 1852.0 #Estimated range using Breguet range equation
#Parameters to engine and combustor performance
numberEngines = ac.parg[igneng]
thrustOneEngine_takeoff = ac.pare[ieFe, iptakeoff, 1] / 1000.0 #Takeoff thrust [kN] 
thrustOneEngine_climb1 = ac.pare[ieFe, ipclimb1, 1] / 1000.0 #Takeoff thrust [kN]
thrustOneEngine_cruise1 = ac.pare[ieFe, ipcruise1, 1] / 1000.0 #Cruise thrust [kN]

var_names = [
    "Payload mass (Ton)",
    "Number of passenger",
    "Flight range (nmi)",
    "PFEI (J/J)",
    "Wing span (m)",
    "Balanced field length (m)",
    "Top of climb flight angle (deg)",
    "Compressor exit temperature (K)",
    "Metal temperature (K)",
    "Fan diameter (m)",
    "Fuel mass (Ton)",
    "Takeoff mass (Ton)",
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
    "Lift-to-drag ratio at cruise",
    "Heating value at cruise (from liquid) (J/kg)",
    "Thrust specific fuel consumption at cruise (kg/s/N)",
    "Thrust specific energy consumption at cruise (J/s/N)",
    "Horizontal velocity at cruise (m/s)",
    "Takeoff mass (Ton)",
    "Fuel mass reserved (Ton)",
    "Fuel mass burned (Ton)",
    "Empty mass (Ton)",
    "Breguet flight range (nmi)",
    "Number of engines",
    "Thrust at start takeoff (one engine) (kN)",
    "Thrust at start climb (one engine) (kN)",
    "Thrust at start cruise (one engine) (kN)"
]

val = [
    massPayload,
    numberPassengers,
    range,
    PFEI,
    wingSpan_Cur,
    fieldLengthBalanced_Cur,
    TOCFlightAngle_Cur,
    Tt3_Cur,
    TMetal_Cur,
    DiaFan_Cur,
    FuelMass_Cur,
    TakeoffMass_Cur,
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
    massFuelReserved,
    massFuelBurned,
    massEmpty,
    rangeBreguet,
    numberEngines,
    thrustOneEngine_takeoff,
    thrustOneEngine_climb1,
    thrustOneEngine_cruise1
]

df = DataFrame(
    Symbol("Variables") => var_names,
    Symbol("Values") => val
)

CSV.write(joinpath(save_dir,"DesignParameters.csv"), df)

# 4) Save the limiting parameters
var_names = [
    "Wing Span (m)",
    "Balanced Field Length (m)",
    "Top of Climb Flight Angle (deg)",
    "Compressor Exit Temperature (K)",
    "Metal Temperature (K)",
    "Fan Diameter (m)",
    "Fuel Mass (Ton)",
    "Takeoff Mass (Ton)"
]

var_cur = [
    wingSpan_Cur,
    fieldLengthBalanced_Cur,
    TOCFlightAngle_Cur,
    Tt3_Cur,
    TMetal_Cur,
    DiaFan_Cur,
    FuelMass_Cur,
    TakeoffMass_Cur
]

var_min = [
    NaN,
    NaN,
    TOCFlightAngle_Min,
    NaN,
    NaN,
    NaN,
    NaN,
    NaN
]

var_max = [
    wingSpan_Max,
    fieldLengthBalanced_Max,
    NaN,
    Tt3_Max,
    TMetal_Max,
    DiaFan_Max,
    FuelMass_Max,
    TakeoffMass_Max
]

df = DataFrame(
    Symbol("Variables") => var_names,
    Symbol("Current Values") => var_cur,
    Symbol("Lower Limits")  => var_min,
    Symbol("Upper Limites")  => var_max,
)
CSV.write(joinpath(save_dir,"LimitingParameters.csv"), df)

# 5) Plot figures
pic = TASOPT.stickfig(ac)
savefig(pic, joinpath(save_dir,"StickPlot.png"))