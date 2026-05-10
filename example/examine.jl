# This is an example file to examine an aircraft model

# 1) Load TASOPT
using TASOPT
using DataFrames, CSV
using Plots
include(__TASOPTindices__)

# 2) Include input file for desired aircraft/
#  load the target model
load_dir = "ModelSaved"
load_name = "acOptimized"
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
massPayload = ac.parm[imWpay, 1]/9.81/1000. #(Ton)
numberPassengers = ac.parm[imWpay, 1]/ac.parm[imWperpax, 1] #Number
Range = ac.parm[imRange, 1]/1852. #nmi
PFEI = ac.parm[imPFEI, 1] #(J/J)

Output = (; 
    Symbol("massPayload (Ton)") => [massPayload],
    Symbol("numberPassengers") => [numberPassengers],
    Symbol("Range (nmi)") => [Range],
    Symbol("PFEI (J/J)") => [PFEI]
)
CSV.write(joinpath(save_dir,"DesignParameters.csv"), Output; writeheader=true)

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