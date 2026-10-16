"""
This script extract the combustor SLS operating condtitions for both R1Size and OAG case.
Using the maximum recorded oag off-design missions's thrust or R1 sizing mission thrust.
Then Sweep from that maximum thrust to minimum thrust like in typical LTO cycle with no engine offtake.
Example for R1 Size:
    Data in: Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_/Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_300.jld2
    Case Dir: "Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_" and caseKey = "Opti_Eth_NoACT_V4_3_R1Sz_EtasEng" and missdir_or_rangelst = collect(300:100:3000)
Example for OAG Size:
    Data in: Opti_Eth_NoACT_OAG_MI_Tail_V3_/Opti_Eth_NoACT_OAG_Ml_4Se_TyC_24Bf_Tail_V3_50.jld2
    Case Dir: caseDir = "Opti_Eth_NoACT_OAG_MI_Tail_V3_" and caseKey = "Opti_Eth_NoACT_OAG_Ml_" and miss_dir = joinpath(model_dir,"OffDesignMissions_50_300_300_Tail.csv")
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(__TASOPTroot__,"../example/UtilityFunctions/objective_factory.jl"))
off_design_specified! = ObjectiveFactory.off_design_specified!
include(joinpath(__TASOPTroot__,"../example/UtilityFunctions/Extract.jl"))
using .Extract: read_oag,extract_combustion_inputs
include(joinpath(__TASOPTroot__,"../example/UtilityFunctions/run_engine.jl"))
using .RunEngine: runOffDes
using Glob

#### Setup IO
# Input case names - OAG seat-capacity sweep
model_dir = joinpath(__TASOPTroot__,"../example/ModelSaved")
# Sized Aircraft Model Directory
caseDir = "Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_"
caseKey = "Opti_Jet_NoACT_V4_4_R1Sz_EtasEng"
# Offdesign Mission Directory (oag seat mode) or range list (R1size range mode)
missdir_or_rangelst = collect(300:100:3000) #joinpath(model_dir,"OffDesignMissions_50_300_300_Tail.csv")
# Seat parameters only used for oag extraction
pass_bulk_frac = 0.825 
pass_tail_frac = 0.850 #For farthest distance
# Save Directory
saveDir = "Combustor_Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_SLS_/"
# Parameters
M0 = 0.0
P0 = 101320.0 #Pa
T0 = 288.2 #K
a0 = 340.2074661144284 #m/s
num_SLS_points = 50 # Number of sea-level-static, no-offtake points to sweep between the min and max off-design thrust
# Specified SLS fuel case. Not used for determining off-design thrust range, but specifically for SLS sweep run given an engine design.
# Overwrite fuel for copied ac model. Just enough to get the engine run using that fuel.
sls_fuel_name     = "Jet"    #"Eth"    , "Du31"   ,"Jet"    for output filename
sls_idx_fuel      = 24       # 32      ,  322431  ,24
sls_hvap_fuel_Jkg = 358694.0 # 918187.9,  586408.0,358694.0

#### Function defintion
"""
All SI units
Sweep through engine off-design performance
"""
function sweep_engine_offdesign(FnMin,FnMax,numPts,ac,M0,P0,T0,a0;zero_offtake=true)
    # Create list of thrust to test
    Fn_Lst_N = range(FnMin,FnMax,length=numPts)
    Fn_Lst_Act_kN = []
    Pt3_Lst_psi = []
    Tt3_Lst_R = []
    Pt4_Lst_psi = []
    Tt4_Lst_R = []
    mdotF_Lst_lbms = []
    mdotA_Lst_lbms = []
    for (idx,Fn_cur) in enumerate(Fn_Lst_N)
        res_cur = runOffDes(ac,M0,P0,T0,a0,Fn_cur;zero_offtake=zero_offtake)
        # Record the performance (Only the converged case)
        if !res_cur.Lconv
            @warn "Point Fn=$(Fn_cur/1000.0) kN did not converge"
            continue
        end
        push!(Fn_Lst_Act_kN, res_cur.Fe/1000.0) #[kN]
        push!(Pt3_Lst_psi, res_cur.pt3/6894.757) #[psi]
        push!(Pt4_Lst_psi, res_cur.pt4/6894.757) #[psi]
        push!(Tt3_Lst_R, res_cur.Tt3*1.8) #[R]
        push!(Tt4_Lst_R, res_cur.Tt4*1.8) #[R]
        push!(mdotF_Lst_lbms, (res_cur.mcore*res_cur.ff)*2.204622) #[lbm/s] single engine combustor fuel flow rate
        push!(mdotA_Lst_lbms, res_cur.mburner*2.204622) #[lbm/s] single engine combustor air flow rate
    end
    return (;Fn_Lst_Act_kN,Pt3_Lst_psi,Tt3_Lst_R,Pt4_Lst_psi,Tt4_Lst_R,mdotF_Lst_lbms,mdotA_Lst_lbms)
end

"""
Extract the maximum thrust from off-design and design mission
Assume the same fuel as design mission but five mission case
With the calibrated tail vs bulk passenger count
"""
function extract_max_min_des_offdes_thrust(ac,ran_Lst_off_nmi,pass_bulk_frac,pass_tail_frac)
    # First get the max thrust at the design mission
    Fn_max_N = maximum(ac.pare[ieFe,:,1])
    Fn_min_N = minimum(ac.pare[ieFe,ipclimb1:ipdescentn,1])
    # Get the payload weights
    wei_pay_off_N  = fill(ac.parg[igWpaymax]*pass_bulk_frac, length(ran_Lst_off_nmi))
    idx_max_range  = argmax(ran_Lst_off_nmi)
    wei_pay_off_N[idx_max_range] = ac.parg[igWpaymax]*pass_tail_frac
    # Test through each off-design (Will crash the code if any of the off-design mission were not converge)
    for (idx,ran_cur_nmi) in enumerate(ran_Lst_off_nmi)
        ac_used = deepcopy(ac)
        res_cur = off_design_specified!(ac_used, ac_used.options.ifuel, ac_used.parg[igrhofuel], ac_used.pare[iehvap,1,1], [ran_cur_nmi], [wei_pay_off_N[idx]]; mod_ac_inplace=true)
        length(res_cur.wei_pay_N)<=0 && error("Offdesign point for $(ran_cur_nmi) nmi did not converge")
        Fn_max_N = max(maximum(ac_used.pare[ieFe,:,2]),Fn_max_N)
        Fn_min_N = min(minimum(ac_used.pare[ieFe,ipclimb1:ipdescentn,2]),Fn_min_N)
    end
    return Fn_min_N,Fn_max_N
end

#### Main Operations
#### Create save directory
saveDirFull  = joinpath(model_dir,saveDir)
mkpath(saveDirFull)

#### Load OAG off-design mission data (ranges/weights per seat capacity)
if !(missdir_or_rangelst isa AbstractVector)
    miss_off_des = read_oag(missdir_or_rangelst)
    seat_cap_keys_all = sort(collect(keys(miss_off_des)))
else
    seat_cap_keys_all = missdir_or_rangelst
end

#### Process the seat groups based on the cycle data provided one at a time. (Skip the one with no aircraft model)
for sc in seat_cap_keys_all
    #### First see if the seat has a sized aircraft
    matches = glob(caseKey*"*_$(sc).jld2",joinpath(model_dir,caseDir))
    length(matches) == 1 || continue
    ac = quickload_aircraft(matches[1])
    println("Run Seat/Range: $(sc)")
    #### Indentify the extreme thrust
    if !(missdir_or_rangelst isa AbstractVector)
        # Load the corresponding missions
        ran_Lst_off_nmi = miss_off_des[sc].ranges_nmi
        # Get the thrust variation
        Fn_min_N, Fn_max_N = extract_max_min_des_offdes_thrust(ac,ran_Lst_off_nmi,pass_bulk_frac,pass_tail_frac)
    else
        Fn_max_N = maximum(ac.pare[ieFe,:,1])
        Fn_min_N = minimum(ac.pare[ieFe,ipclimb1:ipdescentn,1])
    end
    #### Override the fuel used for the SLS sweep itself (independent of the design/off-design primary fuel above)
    ac_fuel = deepcopy(ac)
    ac_fuel.options.ifuel = sls_idx_fuel
    ac_fuel.pare[iehvapcombustor,ipcruise1,1] = sls_hvap_fuel_Jkg
    #### Get the off-design engine performance sweep at SLS
    res_cur = sweep_engine_offdesign(Fn_min_N,Fn_max_N,num_SLS_points,ac_fuel,M0,P0,T0,a0;zero_offtake=true)
    WAR = fill(0.0, length(res_cur.Fn_Lst_Act_kN)) #Default water to air ratio
    #### Save the off-design engine performance
    savePath = joinpath(saveDirFull, splitext(basename(matches[1]))[1]*"_CombSLS_$(sls_fuel_name).csv")
    open(savePath, "w") do io
        println(io, "Thrust[kN],Pt3[psi],Pt4[psi],Tt3[R],Tt4[R],Wf[lbm/s],W3[lbm/s],WAR[m]")
        for idx in eachindex(res_cur.Fn_Lst_Act_kN)
            println(io, "$(res_cur.Fn_Lst_Act_kN[idx]),$(res_cur.Pt3_Lst_psi[idx]),$(res_cur.Pt4_Lst_psi[idx]),$(res_cur.Tt3_Lst_R[idx]),$(res_cur.Tt4_Lst_R[idx]),$(res_cur.mdotF_Lst_lbms[idx]),$(res_cur.mdotA_Lst_lbms[idx]),$(WAR[idx])")
        end
    end
end