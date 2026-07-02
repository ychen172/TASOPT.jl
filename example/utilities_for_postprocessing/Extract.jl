module Extract
export init_results_2Layers,extract_acModel,fill_results!,plot_cases,extract_combustion_inputs
using Plots
using DataFrames, CSV
using TASOPT
include(__TASOPTindices__)
using Statistics

function init_results_2Layers(num_cases::Int, fields::Tuple{Vararg{Symbol}})
    """
    This function create a two layer dictionary
    # num_cases: Int: number of cases for the current field
    # fields: tuple of symbols    
    """
    makevec() = fill!(Vector{Union{Missing,Float64}}(undef, num_cases), missing)
    return Dict{Symbol,Vector{Union{Missing,Float64}}}(f => makevec() for f in fields)
end

function fill_results!(results::Dict{Symbol,Vector{Union{Missing,Float64}}},
                       out::NamedTuple,
                       i::Int)
    for f in keys(results)
        if hasproperty(out, f)   # works for Dict; for NamedTuple, use `f in propertynames(out)` if needed
            results[f][i] = out[f]
        else
            @warn "Missing field $(f)"
        end
    end
end

function extract_acModel(ac_cur,idx_miss)
    """
    This function extract parameters from a aircraft model given a target mission
    Inputs:
        ac_cur: aircraft model
        idx_miss: index of the mission to extract data
    """
    #### overall performance data
    range_cur = ac_cur.parm[imRange, idx_miss] / 1852.0 #(nmi)
    PFEI_cur = ac_cur.parm[imPFEI, idx_miss] #(J/J)
    
    #### mass data
    massTO = ac_cur.parm[imWTO, idx_miss]/gee/1000.0 #Takeoff mass (Ton)
    massFuelTot = ac_cur.parm[imWfuel, idx_miss]/gee/1000.0 #Fuel mass (Ton)(Include reserved and burned)
    massPayload = ac_cur.parm[imWpay, idx_miss]/gee/1000.0 #(Ton)
    massEmpty = massTO - massFuelTot - massPayload #(Ton) empty weight
    massFuelMax = ac_cur.parg[igWfmax]/gee/1000.0 #(Ton) Maximum fuel mass stored in fuel tank
    
    #### fuel volume data
    volFuelTot = ac_cur.parm[imVfuel, idx_miss] #(m3) fuel volume integrated out using various density in flight
    volFuelMax = ac_cur.parg[igVfmax] #(m3) fuel volume get directly from geometry
    FuelVolumeFraction = volFuelTot/volFuelMax #fractional fuel volume used by that mission
    rhoFuelAve = massFuelTot*1000.0/volFuelTot #(kg/m3) averaged fuel density for that mission
    
    #### fuselage geometry
    lenFuseCyl = ac_cur.fuselage.layout.x_end_cylinder - ac_cur.fuselage.layout.x_start_cylinder #(m) might be different from direct reading due to ACT
    volACTFuelMax = ac_cur.fuse_tank.ACT_A*ac_cur.fuse_tank.ACT_l*ac_cur.fuse_tank.ACT_eta_vol
    volWingFuelMax = volFuelMax-volACTFuelMax
    FuelVolumeFractionACT = max(volFuelTot-volWingFuelMax,0)/volFuelTot #Fraction of fuel inside the ACT tank

    #### wing geometry
    span_wing = ac_cur.wing.span #Wing Span [m]
    sweep_wing = ac_cur.wing.layout.sweep #Wing Sweep Angle [deg]
    AR_wing = ac_cur.wing.layout.AR #Wing Aspect Ratio
    ThiCho_wing_in = ac_cur.wing.inboard.cross_section.thickness_to_chord
    ThiCho_wing_out = ac_cur.wing.outboard.cross_section.thickness_to_chord
    taper_wing_in = ac_cur.wing.inboard.λ
    taper_wing_out = ac_cur.wing.outboard.λ
    y_engine = ac_cur.parg[igyeng] #[m] spanwise engine location
    x_engine = ac_cur.parg[igxeng] #[m] axial engine location
    SpanBreak_wing = ac_cur.wing.layout.ηs #Wing span break
    cls_wing = ac_cur.para[iarcls, ipcruise1, 1]
    clt_wing = ac_cur.para[iarclt, ipcruise1, 1]
    PR_fan = ac_cur.pare[iepif, ipcruise1, 1]
    PR_LPC = ac_cur.pare[iepilc, ipcruise1, 1]
    PR_HPC = ac_cur.pare[iepihc,ipcruise1, 1]
    BPR = ac_cur.pare[ieBPR, ipcruise1, 1]
    Tt4 = ac_cur.pare[ieTt4, ipcruise1, 1] #[K]
    AR_vtail = ac_cur.vtail.layout.AR

    #### flight performance data
    CL_cruise = 0.5 * (ac_cur.para[iaCL, ipcruise1, idx_miss] + ac_cur.para[iaCL, ipcruise2, idx_miss]) #Averaged cruise CL
    Alt_cruise = 0.5 * (ac_cur.para[iaalt, ipcruise1, idx_miss] + ac_cur.para[iaalt, ipcruise2, idx_miss]) * 3.280839895 #Averaged cruise altitude [ft]
    LD_cruise = 0.5 * (ac_cur.para[iaCL, ipcruise1, idx_miss]/ac_cur.para[iaCD, ipcruise1, idx_miss] + 
                        ac_cur.para[iaCL, ipcruise2, idx_miss]/ac_cur.para[iaCD, ipcruise2, idx_miss]) #Averaged cruise lift-to-drag ratio
    LHV_cruise = 0.5 * (ac_cur.pare[iehfuel, ipcruise1, idx_miss] + ac_cur.pare[iehfuel, ipcruise2, idx_miss]) #Averaged cruise heating value (J/kg) (Include vaporization already)
    TSFC_cruise = 0.5 * (ac_cur.pare[ieTSFC, ipcruise1, idx_miss] + ac_cur.pare[ieTSFC, ipcruise2, idx_miss]) / gee #Averaged cruise thrust specfic heat consumption (kg/s/N)
    vel_cruise = 0.5 * (cos(ac_cur.para[iagamV, ipcruise1, idx_miss]) * ac_cur.pare[ieu0, ipcruise1, idx_miss] + 
                        cos(ac_cur.para[iagamV, ipcruise2, idx_miss]) * ac_cur.pare[ieu0, ipcruise2, idx_miss]) #Averaged cruise horizontal velocity (m/s)
    eta_total_cruise = (1.0/TSFC_cruise)*(vel_cruise/LHV_cruise) #total cruise engine efficiency

    #### derived data
    energy_flight = PFEI_cur * massPayload * (1000.0 * gee * 1852.0) * range_cur #(J) total flight energy 

    #### parameters for Breguet range
    frac_rese = ac_cur.parg[igfreserve] #W_reserveFuel / W_fuelburned
    PFEI_cru = (LHV_cruise*(ac_cur.para[iafracW, ipcruise1, idx_miss]-ac_cur.para[iafracW, ipcruise2, idx_miss])*ac_cur.parg[igWMTO]/gee)/
                (ac_cur.parm[imWpay, idx_miss] * ac_cur.parm[imRange, idx_miss]) #Cruise Only PFEI (J/J)
    range_cru = ac_cur.para[iaRange, ipcruise2, idx_miss] - ac_cur.para[iaRange, ipcruise1, idx_miss] #(m)
    @assert range_cru>0.0 "Find a negative cruise range, likely range to short"

    #### parameters for engine efficiency at cruise condition
    etaTherm_cru = Float64[]
    spePower_cru = Float64[]
    etaPropu_cru = Float64[]
    OPR_cru = Float64[]
    Tt_TurbIn_cru = Float64[]
    for phase in [ipcruise1,ipcruise2]
        ff_cru = ac_cur.pare[ieff, phase, idx_miss] #mdot_fuel / mdot_core
        BPR_cru = ac_cur.pare[ieBPR, phase, idx_miss] #mdot_BP / mdot_core
        mass_offtake_cru = ac_cur.pare[iemofft, phase, idx_miss] #kg/s single engine
        mass_core_cru = ac_cur.pare[iemcore, phase, idx_miss] #kg/s single engine
        u_coreExh_cru = ac_cur.pare[ieu6, phase, idx_miss] #m/s
        u_fanExh_cru = ac_cur.pare[ieu8, phase, idx_miss] #m/s
        u_inf_cru = ac_cur.pare[ieu0, phase, idx_miss] #m/s
        p_coreExh_cru = ac_cur.pare[iep6, phase, idx_miss] #Pa
        p_fanExh_cru = ac_cur.pare[iep8, phase, idx_miss] #Pa
        p_inf_cru = ac_cur.pare[iep0, phase, idx_miss] #Pa
        A_coreExh_cru = ac_cur.pare[ieA6, phase, idx_miss] #m2
        A_fanExh_cru = ac_cur.pare[ieA8, phase, idx_miss] #m2
        LHV_cru = ac_cur.pare[iehfuel, phase, idx_miss] #J/kg including vaporization heat
        Thrust_cru = ac_cur.pare[ieFsp, phase, idx_miss] * (u_inf_cru * mass_core_cru * (1.0 + BPR_cru)) #N
        Tt41_cru = ac_cur.pare[ieTt41, phase, idx_miss] #(K) turbine inlet temperature after cooling air
        P_Jet_cru = 0.5*(mass_core_cru*(1.0+ff_cru)-mass_offtake_cru)*u_coreExh_cru^2 +
                    0.5*mass_core_cru*BPR_cru*u_fanExh_cru^2 - 
                    0.5*mass_core_cru*(1.0 + BPR_cru)*u_inf_cru^2 + 
                    (p_coreExh_cru-p_inf_cru)*A_coreExh_cru*u_coreExh_cru +
                    (p_fanExh_cru-p_inf_cru)*A_fanExh_cru*u_fanExh_cru #Jet power (J/s)
        push!(OPR_cru, ac_cur.pare[iepid, phase, idx_miss]*ac_cur.pare[iepif, phase, idx_miss]*ac_cur.pare[iepilc, phase, idx_miss]*ac_cur.pare[iepihc, phase, idx_miss])
        push!(etaTherm_cru, P_Jet_cru/(mass_core_cru*ff_cru*LHV_cru)) #Thermal efficiency
        push!(spePower_cru, P_Jet_cru/(mass_core_cru*(1.0+ff_cru+BPR_cru)-mass_offtake_cru)) #J/kg
        push!(etaPropu_cru, (Thrust_cru*u_inf_cru)/P_Jet_cru)
        push!(Tt_TurbIn_cru, Tt41_cru)
    end
    spePower_cru = mean(spePower_cru) #J/kg
    etaTherm_cru = mean(etaTherm_cru)
    etaPropu_cru = mean(etaPropu_cru)            
    OPR_cru = mean(OPR_cru)
    Tt_TurbIn_cru = mean(Tt_TurbIn_cru) #K

    #### Return
    output = (
        range_nmi = range_cur,
        PFEI_JJ = PFEI_cur,
        massTO_Ton = massTO,
        massFuelTot_Ton = massFuelTot,
        massFuelMax_Ton = massFuelMax,
        massPayload_Ton = massPayload,
        massEmpty_Ton = massEmpty,
        rhoFuelAve_kgm3 = rhoFuelAve,
        volFuelTot_m3 = volFuelTot,
        volFuelMax_m3 = volFuelMax,
        LD_cru = LD_cruise,
        LHV_cru_Jkg = LHV_cruise,
        TSFC_cru_kgsN = TSFC_cruise,
        vel_cru_ms = vel_cruise,
        eta_total_cru = eta_total_cruise,
        ene_fli_J = energy_flight,
        frac_rese = frac_rese,
        PFEI_cru_JJ = PFEI_cru,
        range_cru_m = range_cru,
        eta_therm_cru = etaTherm_cru,
        spe_power_cru_Jkg = spePower_cru,
        eta_propu_cru = etaPropu_cru,
        OPR_cru = OPR_cru,
        Tt_turbin_cru_K = Tt_TurbIn_cru,
        FuelVolumeFraction = FuelVolumeFraction,
        lenFuseCyl_m = lenFuseCyl,
        volACTFuelMax_m3 = volACTFuelMax,
        volWingFuelMax_m3 = volWingFuelMax,
        FuelVolumeFractionACT = FuelVolumeFractionACT,
        lenACT_m = ac_cur.fuse_tank.ACT_l,
        span_wing_m = span_wing,
        sweep_wing_deg = sweep_wing,
        AR_wing = AR_wing,
        ThiCho_wing_in = ThiCho_wing_in,
        ThiCho_wing_out = ThiCho_wing_out,
        taper_wing_in = taper_wing_in,
        taper_wing_out = taper_wing_out,
        y_engine_m = y_engine,
        x_engine_m = x_engine,
        SpanBreak_wing = SpanBreak_wing,
        cls_wing = cls_wing,
        clt_wing = clt_wing,
        CL_cruise = CL_cruise,
        Alt_cruise_ft = Alt_cruise,
        PR_fan = PR_fan,
        PR_LPC = PR_LPC,
        PR_HPC = PR_HPC,
        BPR = BPR,
        Tt4 = Tt4,
        AR_vtail = AR_vtail,
    )
    return output
end

const LINESTYLES = [:solid, :dash, :dot, :dashdot, :dashdotdot]
const LINECOLORS = [:blue, :red, :green, :orange, :purple, :brown, :pink, :gray, :black, :cyan,
                     :magenta, :teal, :navy, :maroon, :olive, :gold, :coral, :turquoise, :lime, :indigo]
const MARKERS = [:rect, :circle, :diamond, :utriangle, :dtriangle]
function plot_cases(xlab::String,ylab::String,cases::AbstractVector{<:AbstractDict},xSym::Symbol,ySym::Symbol,datalab::AbstractVector{<:AbstractString},save_name::String,;dpi=800,lw=2,legend::Symbol=:best)
    @assert length(cases) == length(datalab) "cases and datalab length mismatch"
    p = plot(xlabel=xlab, ylabel=ylab, dpi=dpi, legend=legend)
    for (i,case) in enumerate(cases)
        c = LINECOLORS[mod1(i, length(LINECOLORS))]
        plot!(p, case[xSym], case[ySym], marker=MARKERS[mod1(i, length(MARKERS))], mc=c, msc=c, color=c, lw=lw, linestyle=LINESTYLES[mod1(i, length(LINESTYLES))], label=datalab[i])
    end
    savefig(p,save_name)
    return p
end

"""
    extract_combustion_inputs(ac,idx_miss::Int,save_name::AbstractString,save_dir::AbstractString)

`extract_combustion_inputs` extract combustor input parameters and save to csv for Pycaso

    **Inputs**
        - ac: TASOPT model
        - idx_miss: Mission index to extract combustor operating conditions
        - save_name: Name for the csv file without the .csv extension
        - save_dir: Folder to save the csv. Will be created if it did not exist.
    **Ouputs**
        - Saved combustor input conditions
          All parameters are for single engine
"""
function extract_combustion_inputs(ac,idx_miss::Int,save_name::AbstractString,save_dir::AbstractString)
    #### Extract combustor input parameters
    phase_miss = ["C1","C2","C3","C4","C5","R1","R2","D1","D2","D3","D4","D5"]
    time_miss = ac.para[iatime, ipclimb1:ipdescentn, idx_miss]  #[s] Time of the phases
    thrust_miss = ac.pare[ieFe, ipclimb1:ipdescentn, idx_miss] ./ 1000.0 #[kN] Mission thrust
    Pt3_miss = ac.pare[iept3, ipclimb1:ipdescentn, idx_miss] ./ 6894.757 #[psi] combustor inlet pressure
    Pt4_miss = ac.pare[iept4, ipclimb1:ipdescentn, idx_miss] ./ 6894.757 #[psi] combustor outlet pressure
    Tt3_miss = ac.pare[ieTt3, ipclimb1:ipdescentn, idx_miss] .* 1.8 #[R] combustor inlet temperature
    Tt4_miss = ac.pare[ieTt4, ipclimb1:ipdescentn, idx_miss] .* 1.8 #[R] combustor outlet temperature
    _mdot_core_miss_single_ = ac.pare[iemcore, ipclimb1:ipdescentn, idx_miss] #[kg/s] single engine core mass flow range
    mdot_fuel_miss_single = (_mdot_core_miss_single_ .* ac.pare[ieff, ipclimb1:ipdescentn, idx_miss]) .* 2.204622 #[lbm/s] single engine combustor fuel flow rate
    mdot_air_miss_single  = (_mdot_core_miss_single_ .* (1.0 .- ac.pare[iefc, ipclimb1:ipdescentn, idx_miss]) .- ac.pare[iemofft, ipclimb1:ipdescentn, idx_miss]) .* 2.204622 #[lbm/s] single engine combustor air flow rate
    water_air_ratio = fill(0.0, length(time_miss))
    #### Store data into csv
    mission_param = DataFrame()
    mission_param[!, Symbol("Phase")] = phase_miss
    mission_param[!, Symbol("Time[s]")] = time_miss
    mission_param[!, Symbol("Thrust[kN]")] = thrust_miss
    mission_param[!, Symbol("Pt3[psi]")] = Pt3_miss
    mission_param[!, Symbol("Pt4[psi]")] = Pt4_miss
    mission_param[!, Symbol("Tt3[R]")] = Tt3_miss
    mission_param[!, Symbol("Tt4[R]")] = Tt4_miss
    mission_param[!, Symbol("Wf[lbm/s]")] = mdot_fuel_miss_single
    mission_param[!, Symbol("W3[lbm/s]")] =  mdot_air_miss_single
    mission_param[!, Symbol("WAR[m]")] = water_air_ratio
    #### Save csv
    mkpath(save_dir)
    CSV.write(joinpath(save_dir,save_name*".csv"), mission_param)
end

end #module