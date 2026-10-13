module RunEngine
using TASOPT
include(__TASOPTindices__)

"""
This module take the sized engine model from an aicraft model and run the off-design performance of the engine
"""

"""
Results for turbofan off-design operations
"""
mutable struct TFOperRes{T<:AbstractFloat}
    TSFC::T
    Fsp::T
    hfuel::T
    ff::T

    Fe::T
    mcore::T

    pif::T
    pilc::T
    pihc::T

    mbf::T
    mblc::T
    mbhc::T

    Nbf::T
    Nblc::T
    Nbhc::T

    Tt0::T
    ht0::T
    pt0::T
    cpt0::T
    Rt0::T

    Tt18::T
    ht18::T
    pt18::T
    cpt18::T
    Rt18::T

    Tt19::T
    ht19::T
    pt19::T
    cpt19::T
    Rt19::T

    Tt19c::T
    ht19c::T
    pt19c::T
    cpt19c::T
    Rt19c::T

    Tt2::T
    ht2::T
    pt2::T
    cpt2::T
    Rt2::T

    Tt21::T
    ht21::T
    pt21::T
    cpt21::T
    Rt21::T

    Tt25::T
    ht25::T
    pt25::T
    cpt25::T
    Rt25::T

    Tt25c::T
    ht25c::T
    pt25c::T
    cpt25c::T
    Rt25c::T

    Tt3::T
    ht3::T
    pt3::T
    cpt3::T
    Rt3::T

    Tt4::T
    ht4::T
    pt4::T
    cpt4::T
    Rt4::T

    Tt41::T
    ht41::T
    pt41::T
    cpt41::T
    Rt41::T

    Tt45::T
    ht45::T
    pt45::T
    cpt45::T
    Rt45::T

    Tt49::T
    ht49::T
    pt49::T
    cpt49::T
    Rt49::T

    Tt5::T
    ht5::T
    pt5::T
    cpt5::T
    Rt5::T

    Tt7::T
    ht7::T
    pt7::T
    cpt7::T
    Rt7::T

    u0::T

    T2::T
    u2::T
    p2::T
    cp2::T
    R2::T
    M2::T

    T25::T
    u25::T
    p25::T
    cp25::T
    R25::T
    M25::T

    T5::T
    u5::T
    p5::T
    cp5::T
    R5::T
    M5::T

    T6::T
    u6::T
    p6::T
    cp6::T
    R6::T
    M6::T
    A6::T

    T7::T
    u7::T
    p7::T
    cp7::T
    R7::T
    M7::T

    T8::T
    u8::T
    p8::T
    cp8::T
    R8::T
    M8::T
    A8::T

    u9::T
    A9::T

    epf::T
    eplc::T
    ephc::T
    epht::T
    eplt::T

    etaf::T
    etalc::T
    etahc::T
    etaht::T
    etalt::T

    BPR::T
    OPR::T
    mburner::T

    Lconv::Bool
end

"""
For BLI
assume zero Phiinl, Kinl, and no eng_has_BLI_cores
"""
function runOffDes(ac,M0_test,p0_test,T0_test,a0_test,Fe_test; zero_offtake::Bool=false)
    ####Extract parameters
    gee = TASOPT.gee
    # sea-level static inlet conditions
    M0 = M0_test#0.0
    p0 = p0_test#101320.0 #Pa
    T0 = T0_test#288.2 #K
    a0 = a0_test#340.2074661144284 #m/s
    # mision requirements
    Fe = Fe_test#130*1000.0 #ac.pare[ieFe,ipcruise1,1] #[N]thrust to be altered #single engine
    # for scaling
    pref = TASOPT.pref #101320.0 Pa
    Tref = TASOPT.Tref #288.2 K
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
    # off-take (zeroed when matching a no-offtake reference such as ICAO EEDB certification data)
    mofft = zero_offtake ? 0.0 : (ac.parg[igmofWpay] * ac.parg[igWpay] + ac.parg[igmofWMTO] * ac.parg[igWMTO]) / ac.parg[igneng]
    Pofft = zero_offtake ? 0.0 : (ac.parg[igPofWpay] * ac.parg[igWpay] + ac.parg[igPofWMTO] * ac.parg[igWMTO]) / ac.parg[igneng] + ac.pare[ieHXrecircP,ipcruise1,1]
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

    ####Run engine
    out = 
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

    return TFOperRes(out...)
end

end #RunEngine



# #Example to reproduce the inflight engine operations

# using TASOPT
# include(__TASOPTindices__)

# ####Load aircraft model to for engine test
# read_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_/Opti_Jet_NoACT_2900.jld2")
# ac2 = quickload_aircraft(read_dir)

# # print a reference case
# Fe_test1 = ac2.pare[ieFe,ipclimb1,1] #(N)
# M0_test1 = ac2.pare[ieM0,ipclimb1,1]
# p0_test1 = ac2.pare[iep0,ipclimb1,1] #(Pa)
# T0_test1 = ac2.pare[ieT0,ipclimb1,1] #(K)
# a0_test1 = ac2.pare[iea0,ipclimb1,1] #(m/s)
# OPR_test1 = ac2.pare[iepilc,ipclimb1,1]*ac2.pare[iepihc,ipclimb1,1]
# Tt4_test1 = ac2.pare[ieTt4,ipclimb1,1]
# TSFC_test1 = ac2.pare[ieTSFC,ipclimb1,1]
# println("engine thrust: $(Fe_test1/1000.0) N")
# println("Inlet mach number: $(M0_test1)")
# println("Inlet static pressure: $(p0_test1) Pa")
# println("Inlet static temperature: $(T0_test1) K")
# println("Inlet speed of sound: $(a0_test1) m/s")
# println("OPR: $(OPR_test1)")
# println("Tt4: $(Tt4_test1) K")
# println("TSFC: $(TSFC_test1)")

# result = RunEngine.runOffDes(ac2,M0_test1,p0_test1,T0_test1,a0_test1,Fe_test1)

# @assert result.Lconv

# println()
# println("Found engine thrust: $(result.Fe/1000.0) kN")
# println("Found OPR: $(result.pilc*result.pihc)")
# println("Found Tt4: $(result.Tt4) K")
# println("Found TSFC: $(result.TSFC)")