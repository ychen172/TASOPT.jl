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
include("../src/data_structs/index.inc")

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

n = min(length(Ranges_Eth), length(Ranges_Jet))
Ranges_Com = Ranges_Eth[1:n-2] #[nmi] #Drop the last two to avoid 0 range ending point
mPay_Com = min.(mPay_Jet[1:n-2], mPay_Eth[1:n-2])  #[Ton] smaller payload of the two
mPay_Com .*= 0.9
"""
#Run off-design on jet fuel (Specified common ranges and payloads)
idxFuel = 25
rhoFuel = 817.0
LHVaporFuel = 358694.0 #J/kg
mPay_Jet_Com, Ranges_Jet_Com, PFEIs_Jet_Com, EneTO_Jet, EneCR_Jet, EneDE_Jet = PayloadRangeSpecified(ac, idxFuel, rhoFuel, LHVaporFuel, mPay_Com, Ranges_Com)

Output = (; 
    Symbol("mPay (Ton)") => mPay_Jet_Com,
    Symbol("Ranges (nmi)") => Ranges_Jet_Com,
    Symbol("PFEIs") => PFEIs_Jet_Com,
    Symbol("EnergyTO (J)") => EneTO_Jet,
    Symbol("EnergyCR (J)") => EneCR_Jet,
    Symbol("EnergyDE (J)") => EneDE_Jet
)
CSV.write("$(savedir)PayRanDataJetFuel_Com.csv", Output; writeheader=true)

#Run off-design on ethanol (Specified common ranges and payloads)
idxFuel = 32
rhoFuel = 789.0
LHVaporFuel = 918187.9 #J/kg
mPay_Eth_Com, Ranges_Eth_Com, PFEIs_Eth_Com, EneTO_Eth, EneCR_Eth, EneDE_Eth = PayloadRangeSpecified(ac, idxFuel, rhoFuel, LHVaporFuel, mPay_Com, Ranges_Com)

Output = (; 
    Symbol("mPay (Ton)") => mPay_Eth_Com,
    Symbol("Ranges (nmi)") => Ranges_Eth_Com,
    Symbol("PFEIs") => PFEIs_Eth_Com,
    Symbol("EnergyTO (J)") => EneTO_Eth,
    Symbol("EnergyCR (J)") => EneCR_Eth,
    Symbol("EnergyDE (J)") => EneDE_Eth
)
CSV.write("$(savedir)PayRanDataEthanol_Com.csv", Output; writeheader=true)
"""

idxFuelPri = 25
idxFuelSec = 32
rhoFuelPri = 817.0
rhoFuelSec = 789.0
hVapFuelPri = 358694.0
hVapFuelSec = 918187.9
flgPhaseSwitch = zeros(size(ac.pare, 2))
flgPhaseSwitch[:] .= 1.0
#Run dual-fuel
mPay_Dua_Com, Ranges_Dua_Com, PFEIs_Dua_Com, EneTO_Dua, EneCR_Dua, EneDE_Dua = PayloadRangeSpecDual(ac, idxFuelPri, rhoFuelPri, hVapFuelPri, mPay_Com, Ranges_Com, idxFuelSec, rhoFuelSec, hVapFuelSec, flgPhaseSwitch)
Output = (; 
    Symbol("mPay (Ton)") => mPay_Dua_Com,
    Symbol("Ranges (nmi)") => Ranges_Dua_Com,
    Symbol("PFEIs") => PFEIs_Dua_Com,
    Symbol("EnergyTO (J)") => EneTO_Dua,
    Symbol("EnergyCR (J)") => EneCR_Dua,
    Symbol("EnergyDE (J)") => EneDE_Dua
)
CSV.write("$(savedir)PayRanDataDual_Com.csv", Output; writeheader=true)