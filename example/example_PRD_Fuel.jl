"""
An example script to plot a payload-range comparison between jet fuel and ethanol
fleet of missions
"""

# 1. Import modules
using TASOPT
using Plots
using CSV
# Get directory
savedir = "Movie/"
if !isdir(savedir)
    mkdir(savedir)
end


# Load default model
ac = load_default_model() #Use default model for payload-range diagram
size_aircraft!(ac)

#Run off-design on jet fuel
idxFuel = 25
rhoFuel = 817.0
LHVaporFuel = 358694.0 #J/kg
mPay_Jet, Ranges_Jet, PFEIs_Jet = PayloadRangeFuel(ac, idxFuel, rhoFuel, LHVaporFuel, Rpts = 40, Ppts = 41)

Output = (; 
    Symbol("mPay (Ton)") => mPay_Jet,
    Symbol("Ranges (nmi)") => Ranges_Jet,
    Symbol("PFEIs") => PFEIs_Jet
)
CSV.write("$(savedir)PayRanDataJetFuel.csv", Output; writeheader=true)

idxFuel = 32
rhoFuel = 789.0
LHVaporFuel = 918187.9 #J/kg
mPay_Eth, Ranges_Eth, PFEIs_Eth = PayloadRangeFuel(ac, idxFuel, rhoFuel, LHVaporFuel, Rpts = 40, Ppts = 41)

Output = (; 
    Symbol("mPay (Ton)") => mPay_Eth,
    Symbol("Ranges (nmi)") => Ranges_Eth,
    Symbol("PFEIs") => PFEIs_Eth
)
CSV.write("$(savedir)PayRanDataEthanol.csv", Output; writeheader=true)