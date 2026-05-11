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
Range = ac.parm[imRange, 1]/1852. #nmi
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
LHV_cruise = 0.5 * (ac.pare[iehfuel, ipcruise1, 1] + ac.pare[iehfuel, ipcruise2, 1]) #Averaged cruise heating value (J/kg)
TSFC_cruise = 0.5 * (ac.pare[ieTSFC, ipcruise1, 1] + ac.pare[ieTSFC, ipcruise2,1]) / gee #Averaged cruise thrust specfic heat consumption (kg/s/N)
TSEC_cruise = TSFC_cruise*LHV_cruise #Averaged cruise thrust specific energy consumption (J/s/N)
vel_cruise = 0.5 * (cos(ac.para[iagamV, ipcruise1,1]) * ac.pare[ieu0, ipcruise1,1] + 
                    cos(ac.para[iagamV, ipcruise2,1]) * ac.pare[ieu0, ipcruise2, 1]) #Averaged cruise horizontal velocity (m/s)
massTO = ac.parm[imWTO,1]/gee/1000.0 #Takeoff mass (Ton)
massFuel = ac.parm[imWfuel]/gee/1000.0 #Fuel mass (Ton)
massEmpty = massTO-massFuel-massPayload #By definition, the empty weight (Ton)


var_names = [
    "Payload Mass (Ton)",
    "Number of passenger",
    "Flight Range (nmi)",
    "PFEI (J/J)",
    "Wing Span (m)",
    "Balanced Field Length (m)",
    "Top of Climb Flight Angle (deg)",
    "Compressor Exit Temperature (K)",
    "Metal Temperature (K)",
    "Fan Diameter (m)",
    "Fuel Mass (Ton)",
    "Takeoff Mass (Ton)",
    "Wing Aspect Ratio",
    "Cruise CL",
    "Wing Sweep Angle (deg)",
    "Cruise Altitude (ft)",
    "Inboard Wing Taper Ratio",
    "Outboard Wing Taper Ratio",
    "Inboard Thickness-to-Chord Ratio",
    "Outboard Thickness-to-Chord Ratio",
    "Break/Root CL Ratio at Cruise",
    "Tip/Root CL Ratio at Cruise",
    "Tt4 at Cruise (K)",
    "HPC Pressure Ratio at Cruise",
    "Fan Pressure Ratio at Cruise",
    "LPC Pressure Ratio at Cruise",
    "Bypass Ratio at Cruise",
    "Lift-to-Drag Ratio at Cruise",
    "Heating Value at Cruise (J/kg)",
    "Thrust Specific Fuel Consumption at Cruise (kg/s/N)",
    "Thrust Specific Energy Consumption at Cruise (J/s/N)",
    "Horizontal Velocity at Cruise (m/s)"
]

val = [
    massPayload,
    numberPassengers,
    Range,
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
    vel_cruise
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