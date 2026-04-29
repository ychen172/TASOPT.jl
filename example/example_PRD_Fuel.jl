"""
An example script to plot a payload-range comparison between jet fuel and ethanol
fleet of missions
"""

# 1. Import modules
using TASOPT
using Plots
# import indices for calling parameters

# Load default model
ac = load_default_model() #Use default model for payload-range diagram
size_aircraft!(ac)

#Run off-design on jet fuel
idxFuel = 25
rhoFuel = 817.0
LHVaporFuel = 358694.0 #J/kg
mPay_Lst, Ranges_Lst, PFEIs_Lst = PayloadRangeFuel(ac, idxFuel, rhoFuel, LHVaporFuel)
println(Ranges_Lst,"\n")
println(mPay_Lst,"\n")
println(PFEIs_Lst,"\n")


idxFuel = 32
rhoFuel = 789.0
LHVaporFuel = 918187.9 #J/kg
mPay_Lst, Ranges_Lst, PFEIs_Lst = PayloadRangeFuel(ac, idxFuel, rhoFuel, LHVaporFuel)
println(Ranges_Lst,"\n")
println(mPay_Lst,"\n")
println(PFEIs_Lst,"\n")