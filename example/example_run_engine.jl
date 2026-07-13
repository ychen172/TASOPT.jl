"""
This script intends to probe the sea-level static operation performance of engine model optimized for different aircraft design
"""

using TASOPT
include(__TASOPTindices__)

####Load aircraft model to for engine test
read_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_/Opti_Jet_NoACT_2900.jld2")
ac = quickload_aircraft(read_dir)

####Extract parameters
gee = TASOPT.gee
# sea-level static inlet conditions
M0 = 0.0
p0 = 101320.0 #Pa
T0 = 288.2 #K
a0 = 340.2074661144284 #m/s
# for scaling
pref = 101320.0 #Pa
Tref = 288.2 #K
# for BLI
Phiinl = 0.0
Kinl = 0.0
eng_has_BLI_cores = false #For 0 speed, having core or not does not matter
# about the engine design (constant pressure ratio)
pid = ac.pare[iepid,ipcruise1,1] #diffuser pressure ratio (0.998)
pib = ac.pare[iepib,ipcruise1,1] #burner (0.94)
pifn = ac.pare[iepifn,ipcruise1,1] #fan nozzle (0.98)
pitn = ac.pare[iepitn,ipcruise1,1] #turbine nozzle (0.989)
epsl = ac.pare[ieepsl,ipcruise1,1] #low spool shaft power loss fraction 1%
epsh = ac.pare[ieepsh,ipcruise1,1] #high spool shaft power loss fraction 2.2%
etab = ac.pare[ieetab,ipcruise1,1] #combustion efficiency 98%
# about the engine design (gear ratio)
Gearf = ac.parg[igGearf] #gear ratio (1.0) for no gear
# about the engine design (design point pressure ratio, mass flow, and speed for map scaling purpose)
pifD  = ac.pare[iepifD,ipcruise1,1]
pilcD = ac.pare[iepilcD,ipcruise1,1]
pihcD = ac.pare[iepihcD,ipcruise1,1]
pihtD = ac.pare[iepihtD,ipcruise1,1]
piltD = ac.pare[iepiltD,ipcruise1,1]
#
mbfD  = ac.pare[iembfD,ipcruise1,1]
mblcD = ac.pare[iemblcD,ipcruise1,1]
mbhcD = ac.pare[iembhcD,ipcruise1,1]
mbhtD = ac.pare[iembhtD,ipcruise1,1]
mbltD = ac.pare[iembltD,ipcruise1,1]
#
NbfD  = ac.pare[ieNbfD,ipcruise1,1]
NblcD = ac.pare[ieNblcD,ipcruise1,1]
NbhcD = ac.pare[ieNbhcD,ipcruise1,1]
NbhtD = ac.pare[ieNbhtD,ipcruise1,1]
NbltD = ac.pare[ieNbltD,ipcruise1,1]
#
epolf  = ac.pare[ieepolf,ipcruise1,1] #89.48
epollc = ac.pare[ieepollc,ipcruise1,1] #88
epolhc = ac.pare[ieepolhc,ipcruise1,1] #87
epolht = ac.pare[ieepolht,ipcruise1,1] #88.9
epollt = ac.pare[ieepollt,ipcruise1,1] #89.9
# about the engine design (fixed geometry)
A2  = ac.pare[ieA2,ipcruise1,1]
A25 = ac.pare[ieA25,ipcruise1,1]
A5  = ac.pare[ieA5,ipcruise1,1]
A7  = ac.pare[ieA7,ipcruise1,1]
# for off-design mission operation, use thrust as mision requirement
opt_calc_call = TASOPT.engine.CalcMode.FixedFeOffDes
opt_cooling = TASOPT.engine.CoolingOpt.FixedCoolingFlowRatio
# about the engine design (fuel and combustor)
Tfuel = ac.pare[ieTfuel,ipcruise1,1] #280 K
ifuel = ac.options.ifuel
hvap = ac.pare[iehvapcombustor,ipcruise1,1] #358694.0 J/kg for jet fuel
# off-take
mofft = (ac.parg[igmofWpay] * ac.parg[igWpay] + ac.parg[igmofWMTO] * ac.parg[igWMTO]) / ac.parg[igneng]
Pofft = (ac.parg[igPofWpay] * ac.parg[igWpay] + ac.parg[igPofWMTO] * ac.parg[igWMTO]) / ac.parg[igneng] + ac.pare[ieHXrecircP,ipcruise1,1]
Tt9 = ac.pare[ieTt9,ipcruise1,1] #[K] offtake air discharge total temperature
pt9 = ac.pare[iept9,ipcruise1,1] #[Pa] offtake air discharge total pressure
# cooling
Mtexit = ac.pare[ieMtexit,ipcruise1,1] #turbine exit mach number for temperature calculation
dTstrk = ac.pare[iedTstrk,ipcruise1,1] #temperatrue gradient for heat transfer
StA = ac.pare[ieStA,ipcruise1,1] #Staton number for heat transfer
efilm = ac.pare[ieefilm,ipcruise1,1] #cooling efficiency
tfilm = ac.pare[ietfilm,ipcruise1,1] #film effectiveness
fc0   = ac.pare[iefc0,ipcruise1,1] #design point fraction cooling flow fraction for efficiency scaling
epht_fc = ac.pare[iedehtdfc,ipcruise1,1] #gradient of hpt efficiency with respect to the cooling flow fraction
M4a = ac.pare[ieM4a,ipcruise1,1] #cooling flow outlet effective mach number
ruc = ac.pare[ieruc,ipcruise1,1] #cooling flow outlet velocity ratio
ncrow = ncrowx #both are the number of blade rows to be coolled directly fixed by the tasopt indices (4 rows)
epsrow = ac.pare[ieepsc1:(ieepsc1+ncrowx-1),ipcruise1,1] #mass fraction of cooling air for the 4 rows (Driving parameters)
Tmrow = ac.pare[ieTmet1:(ieTmet1+ncrowx-1),ipcruise1,1] #metal temprature for the 4 rows              (Driven parameters to be updated)
# mision requirements
Fe = 130*1000.0 #ac.pare[ieFe,ipcruise1,1] #[N]thrust to be altered #single engine
# Initial guesses to be iterated
M2 = ac.pare[ieM2,ipcruise1,1]
M25 = ac.pare[ieM25,ipcruise1,1]
pif = max(ac.pare[iepif,ipcruise1,1], 1.1)
pilc = max(ac.pare[iepilc,ipcruise1,1], 1.1)
pihc = max(ac.pare[iepihc,ipcruise1,1], 1.1)
mbf = ac.pare[iembf,ipcruise1,1]
mblc = ac.pare[iemblc,ipcruise1,1]
mbhc = ac.pare[iembhc,ipcruise1,1]
Tt4 = ac.pare[ieTt4,ipcruise1,1]
pt5 = ac.pare[iept5,ipcruise1,1]
mcore = ac.pare[iemcore,ipcruise1,1]
#Heat exchanger variables (all 0 here dont worry)
Δh_PreC = ac.pare[iePreCDeltah,ipcruise1,1]
Δh_InterC = ac.pare[ieInterCDeltah,ipcruise1,1]
Δh_Regen = ac.pare[ieRegenDeltah,ipcruise1,1]
Δh_TurbC = ac.pare[ieTurbCDeltah,ipcruise1,1]
Δp_PreC = ac.pare[iePreCDeltap,ipcruise1,1]
Δp_InterC = ac.pare[ieInterCDeltap,ipcruise1,1]
Δp_Regen = ac.pare[ieRegenDeltap,ipcruise1,1]

println("ini Fe(kN): $(Fe/1000.0)")
println("ini OPR: $(pihc*pilc)")
println("ini Tt4: $(Tt4)")

####Run engine
TSFC, Fsp, hfuel, ff,
Fe, mcore,
pif, pilc, pihc,
mbf, mblc, mbhc,
Nbf, Nblc, Nbhc,
Tt0, ht0, pt0, cpt0, Rt0,
Tt18, ht18, pt18, cpt18, Rt18,
Tt19, ht19, pt19, cpt19, Rt19,
Tt19c, ht19c, pt19c, cpt19c, Rt19c,
Tt2, ht2, pt2, cpt2, Rt2,
Tt21, ht21, pt21, cpt21, Rt21,
Tt25, ht25, pt25, cpt25, Rt25,
Tt25c, ht25c, pt25c, cpt25c, Rt25c,
Tt3, ht3, pt3, cpt3, Rt3,
Tt4, ht4, pt4, cpt4, Rt4,
Tt41, ht41, pt41, cpt41, Rt41,
Tt45, ht45, pt45, cpt45, Rt45,
Tt49, ht49, pt49, cpt49, Rt49,
Tt5, ht5, pt5, cpt5, Rt5,
Tt7, ht7, pt7, cpt7, Rt7,
u0,
T2, u2, p2, cp2, R2, M2,
T25, u25, p25, cp25, R25, M25,
T5, u5, p5, cp5, R5, M5,
T6, u6, p6, cp6, R6, M6, A6,
T7, u7, p7, cp7, R7, M7,
T8, u8, p8, cp8, R8, M8, A8,
u9, A9,
epf, eplc, ephc, epht, eplt,
etaf, etalc, etahc, etaht, etalt,
Lconv = 
TASOPT.engine.tfoper!(
gee, M0, T0, p0, a0, Tref, pref,
Phiinl, Kinl, eng_has_BLI_cores,
pid, pib, pifn, pitn,
Gearf,
pifD, pilcD, pihcD, pihtD, piltD,
mbfD, mblcD, mbhcD, mbhtD, mbltD,
NbfD, NblcD, NbhcD, NbhtD, NbltD,
A2, A25, A5, A7,
opt_calc_call,
Tfuel, ifuel, hvap, etab,
epolf, epollc, epolhc, epolht, epollt,
mofft, Pofft,
Tt9, pt9,
epsl, epsh,
opt_cooling,
Mtexit, dTstrk, StA, efilm, tfilm,
fc0, epht_fc,
M4a, ruc,
ncrowx, ncrow,
epsrow, Tmrow,
Fe,
M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
Δp_PreC, Δp_InterC, Δp_Regen)

println()
println("Fe(kN): $(Fe/1000.0)")
println("OPR: $(pihc*pilc)")
println("Tt4: $(Tt4)")

