"""
This script compare the off-design performance between cases in terms of PFEI
and flight range
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
using Plots

#### Helpers
function interp_range_at_threshold(range, y, i_peak, eps; side::Symbol)
    yth = y[i_peak] * (1 - eps)

    if side == :right
        # first point to the RIGHT of peak below threshold
        k2 = findnext(v -> v < yth, y, i_peak + 1)
        isnothing(k2) && return nothing
        k1 = k2 - 1
    elseif side == :left
        # first point to the LEFT of peak below threshold
        k1 = findprev(v -> v < yth, y, i_peak - 1)
        isnothing(k1) && return nothing
        k2 = k1 + 1
    else
        error("side must be :right or :left")
    end

    x1, x2 = range[k1], range[k2]
    y1, y2 = y[k1], y[k2]

    # linear interpolation: x at y=yth
    if y2 == y1
        return x1  # degenerate case
    end
    xth = x1 + (yth - y1) * (x2 - x1) / (y2 - y1)
    return xth
end

#### Setup IO
# I
case_keywords = ["off_designJet", "jetfuel_match_payload", "jetfuel_to_ethanolJet"]
case_names    = ["Base", "Match", "Retrofit"]
model_dir     = "ModelSaved"
des_ranges    = [3000] #float.(collect(300:100:3000)) #design range to compare (nmi)
# O
save_dir      = "ModelProcessed"
save_name     = "compare" #sub_folder will be created
# Test conditions
offdes_ranges = float.(collect(0:100:8000)) # (nmi)
idx_R1R2IdxCorrection = [[2,3]] #3->2
range_comparison_index = [3,1] #(case 1 - case 3) / case1
epsR1R2 = 0.005

#### Save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
range_lst = [] #(nmi)
PFEI_lst = [] #(J/J)
massEmp_lst = [] #(Ton)
voluFuel_lst = [] #(m3)
voluFuelMax_lst = [] #(m3)
massTO_lst = [] #(Ton)
massTOMax_lst = [] #(Ton)
massPay_lst = [] #(Ton)
LD_cru_lst = [] #Cruise lift to drag ratio
eta_tot_cru_lst = [] #Total engine efficiency at cruise
LHV_lst = [] #Cruise heating value including evaporation (J/kg)
idx_R1_lst = []
idx_R2_lst = []
R1_lst = []
R2_lst = []

#### Extract data for each scenerio
for (i, keyword_cur) in enumerate(case_keywords)
    model_dir_sub = joinpath(model_dir, keyword_cur)
    range_lst_sub = [] #(nmi)
    PFEI_lst_sub = [] #(J/J)
    massEmp_lst_sub = [] #(Ton)
    voluFuel_lst_sub = [] #(m3)
    voluFuelMax_lst_sub = [] #(m3)
    massTO_lst_sub = [] #(Ton)
    massTOMax_lst_sub = [] #(Ton)
    massPay_lst_sub = [] #(Ton)
    LD_cru_lst_sub = [] #Cruise lift to drag ratio
    eta_tot_cru_lst_sub = [] #Total engine efficiency at cruise
    LHV_lst_sub = [] #Cruise heating value including evaporation (J/kg)
    idx_R1_lst_sub = []
    idx_R2_lst_sub = []
    R1_lst_sub = []
    R2_lst_sub = []
    # Extract data for each design case
    for (j, des_range_cur) in enumerate(des_ranges)
        model_dir_sub_sub = joinpath(model_dir_sub, keyword_cur*"$(round(Int,des_range_cur))")
        range_lst_sub_sub = [] #(nmi)
        PFEI_lst_sub_sub = [] #(J/J)
        massEmp_lst_sub_sub = [] #(Ton)
        voluFuel_lst_sub_sub = [] #(m3)
        voluFuelMax_lst_sub_sub = [] #(m3)
        massTO_lst_sub_sub = [] #(Ton)
        massTOMax_lst_sub_sub = [] #(Ton)
        massPay_lst_sub_sub = [] #(Ton)
        LD_cru_lst_sub_sub = [] #Cruise lift to drag ratio
        eta_tot_cru_lst_sub_sub = [] #Total engine efficiency at cruise
        LHV_lst_sub_sub = [] #Cruise heating value including evaporation (J/kg)
        ## extract for each off-design case
        for (k, offdes_range_cur) in enumerate(offdes_ranges)
            model_dir_sub_sub_sub = joinpath(model_dir_sub_sub, keyword_cur*"$(round(Int,des_range_cur))_$(round(Int,offdes_range_cur)).jld2")
            ## reading
            println("Attempt to read from $(model_dir_sub_sub_sub)")
            try
                ac_cur = quickload_aircraft(model_dir_sub_sub_sub)
                ## extract data
                range_cur = ac_cur.parm[imRange,2] / 1852.0 #(nmi)
                PFEI_cur = ac_cur.parm[imPFEI, 2] #(J/J)
                
                ## Check on masses
                massTO = ac_cur.parm[imWTO,2]/gee/1000.0 #Takeoff mass (Ton)
                massFuelTot = ac_cur.parm[imWfuel,2]/gee/1000.0 #Fuel mass (Ton)(Include reserved and burned)
                massPayload = ac_cur.parm[imWpay, 2]/gee/1000.0 #(Ton)
                massEmpty = massTO - massFuelTot - massPayload #(Ton) empty weight
                massTOMax = ac_cur.parg[igWMTO] / gee / 1000.0 #Maximum takeoff mass (Ton)
                ## Check on volume
                rhoFuel = ac_cur.parg[igrhofuel] #kg/m3
                volFuel = massFuelTot * 1000.0 / rhoFuel #m3
                volFuelMax = ac_cur.parg[igWfmax] / gee / rhoFuel #m3 (The design mission fuel mass might be different from the maximum fuel mass with off-design fuel density)
                ## Check on flight performance
                LD_cruise = 0.5 * (ac_cur.para[iaCL, ipcruise1, 2]/ac_cur.para[iaCD, ipcruise1, 2] + 
                                   ac_cur.para[iaCL, ipcruise2, 2]/ac_cur.para[iaCD, ipcruise2, 2]) #Averaged cruise lift-to-drag ratio
                LHV_cruise = 0.5 * (ac_cur.pare[iehfuel, ipcruise1, 2] + ac_cur.pare[iehfuel, ipcruise2, 2]) #Averaged cruise heating value (J/kg) (Include vaporization already)
                TSFC_cruise = 0.5 * (ac_cur.pare[ieTSFC, ipcruise1, 2] + ac_cur.pare[ieTSFC, ipcruise2, 2]) / gee #Averaged cruise thrust specfic heat consumption (kg/s/N)
                vel_cruise = 0.5 * (cos(ac_cur.para[iagamV, ipcruise1, 2]) * ac_cur.pare[ieu0, ipcruise1, 2] + 
                                    cos(ac_cur.para[iagamV, ipcruise2, 2]) * ac_cur.pare[ieu0, ipcruise2, 2]) #Averaged cruise horizontal velocity (m/s)
                eta_total_cruise = (1.0/TSFC_cruise)*(vel_cruise/LHV_cruise) #total cruise engine efficiency

                ## store
                push!(range_lst_sub_sub, range_cur)
                push!(PFEI_lst_sub_sub, PFEI_cur)
                push!(massEmp_lst_sub_sub, massEmpty)
                push!(voluFuel_lst_sub_sub, volFuel)
                push!(voluFuelMax_lst_sub_sub, volFuelMax)
                push!(massTO_lst_sub_sub, massTO)
                push!(massTOMax_lst_sub_sub, massTOMax)
                push!(massPay_lst_sub_sub, massPayload)
                push!(LD_cru_lst_sub_sub, LD_cruise)
                push!(eta_tot_cru_lst_sub_sub, eta_total_cruise)
                push!(LHV_lst_sub_sub, LHV_cruise)
            catch
                nothing
            end
        end
        # Range threshold identification
        m,l = findmax(massPay_lst_sub_sub)
        idx_R1 = findnext(x -> x < m*(1.0-epsR1R2), massPay_lst_sub_sub, l) - 1 #index for R1 range (Not applicable to matching case)
        R1 = interp_range_at_threshold(range_lst_sub_sub, massPay_lst_sub_sub, l, epsR1R2; side=:right)
        m,l = findmax(voluFuel_lst_sub_sub)
        idx_R2 = findprev(x -> x < m*(1.0-epsR1R2), voluFuel_lst_sub_sub, l) + 1 #index for R2 range (Matching case should follow retrofit)
        R2 = interp_range_at_threshold(range_lst_sub_sub, voluFuel_lst_sub_sub, l, epsR1R2; side=:left)
        # Store
        push!(range_lst_sub,range_lst_sub_sub)
        push!(PFEI_lst_sub,PFEI_lst_sub_sub)
        push!(massEmp_lst_sub,massEmp_lst_sub_sub)
        push!(voluFuel_lst_sub,voluFuel_lst_sub_sub)
        push!(voluFuelMax_lst_sub,voluFuelMax_lst_sub_sub)
        push!(massTO_lst_sub,massTO_lst_sub_sub)
        push!(massTOMax_lst_sub,massTOMax_lst_sub_sub)
        push!(massPay_lst_sub,massPay_lst_sub_sub)
        push!(LD_cru_lst_sub,LD_cru_lst_sub_sub)
        push!(eta_tot_cru_lst_sub,eta_tot_cru_lst_sub_sub)
        push!(LHV_lst_sub,LHV_lst_sub_sub)
        push!(idx_R1_lst_sub,idx_R1)
        push!(idx_R2_lst_sub,idx_R2)
        push!(R1_lst_sub,R1)
        push!(R2_lst_sub,R2)
    end
    # Store
    push!(range_lst,range_lst_sub)
    push!(PFEI_lst,PFEI_lst_sub)
    push!(massEmp_lst,massEmp_lst_sub)
    push!(voluFuel_lst,voluFuel_lst_sub)
    push!(voluFuelMax_lst,voluFuelMax_lst_sub)
    push!(massTO_lst,massTO_lst_sub)
    push!(massTOMax_lst,massTOMax_lst_sub)
    push!(massPay_lst,massPay_lst_sub)
    push!(LD_cru_lst,LD_cru_lst_sub)
    push!(eta_tot_cru_lst,eta_tot_cru_lst_sub)
    push!(LHV_lst,LHV_lst_sub)
    push!(idx_R1_lst,idx_R1_lst_sub)
    push!(idx_R2_lst,idx_R2_lst_sub)
    push!(R1_lst,R1_lst_sub)
    push!(R2_lst,R2_lst_sub)
end

#### Correct the R1 R2 Index for the matching case
for idxSwap in idx_R1R2IdxCorrection
    idx_R1_lst[idxSwap[1]] = idx_R1_lst[idxSwap[2]]
    idx_R2_lst[idxSwap[1]] = idx_R2_lst[idxSwap[2]]
end

#### Calculate range reduction
idx_Rcomp_base = range_comparison_index[2]
idx_Rcomp_targ = range_comparison_index[1]
frac_change_R1 = []
frac_change_R2 = []
des_range_R1R2 = [] #extracted design range just to compare R1 and R2
for (i, des_range_cur) in enumerate(des_ranges)
    if (idx_R1_lst[idx_Rcomp_base][i] != 1) && (idx_R1_lst[idx_Rcomp_targ][i] != 1) #R1 touch the left bound
        R1_base = R1_lst[idx_Rcomp_base][i]
        R1_targ = R1_lst[idx_Rcomp_targ][i]
        push!(frac_change_R1,(R1_targ-R1_base)/(R1_base))
        R2_base = R2_lst[idx_Rcomp_base][i]
        R2_targ = R2_lst[idx_Rcomp_targ][i]
        push!(frac_change_R2,(R2_targ-R2_base)/(R2_base))
        push!(des_range_R1R2, des_range_cur)
    end
end

#### Calculate PFEI reduction
PFEI_change = []
PFEI_ranges = []
PFEI_R1R2_idx = []
energy_flight_change = [] 
for (i, des_range_cur) in enumerate(des_ranges)
    # Extract and PFEI and range
    range_base = range_lst[idx_Rcomp_base][i]
    range_targ = range_lst[idx_Rcomp_targ][i]
    PFEI_base  = PFEI_lst[idx_Rcomp_base][i] #vector of off-design
    PFEI_targ  = PFEI_lst[idx_Rcomp_targ][i] #vector of off-design
    mass_payload_base = massPay_lst[idx_Rcomp_base][i] #Ton
    mass_payload_targ = massPay_lst[idx_Rcomp_targ][i]
    # Filter the PFEI and range for common subset
    min_range  = max(minimum(range_base),minimum(range_targ)) #common range bound (assume same spacing)
    max_range  = min(maximum(range_base),maximum(range_targ))
    msk_base   = (range_base .>= min_range) .& (range_base .<= max_range)
    msk_targ   = (range_targ .>= min_range) .& (range_targ .<= max_range)
    range_base = range_base[msk_base] 
    PFEI_base  = PFEI_base[msk_base] 
    range_targ = range_targ[msk_targ] 
    PFEI_targ  = PFEI_targ[msk_targ]
    mass_payload_base = mass_payload_base[msk_base] #Ton
    mass_payload_targ = mass_payload_targ[msk_targ]
    # find the PFEI change
    push!(PFEI_change, (PFEI_targ .- PFEI_base) ./ PFEI_base)
    push!(PFEI_ranges, range_base)
    # Identify the point closes to R1 and closest to R2
    R1_targ = R1_lst[idx_Rcomp_targ][i]
    R2_targ = R2_lst[idx_Rcomp_targ][i]
    push!(PFEI_R1R2_idx, [argmin(abs.(range_base .- R1_targ)) , argmin(abs.(range_base .- R2_targ))])
    # find the flight energy from payload mass
    energy_flight_base = PFEI_base .* mass_payload_base .* (1000.0*gee*1852.0) .* range_base #(J)
    energy_flight_target = PFEI_targ .* mass_payload_targ .* (1000.0*gee*1852.0) .* range_targ #(J)
    push!(energy_flight_change, (energy_flight_target .- energy_flight_base) ./ energy_flight_base)
end

#### Plotting
# PFEI
plot_PFEI = plot(xlabel="Range (nmi)", ylabel="PFEI (J/J)", dpi=800, yscale=:log10)
# plot!(plot_PFEI, xlims=(100,1600),ylims=(0.6, 1.0))
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_PFEI, range_lst[i][j], PFEI_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_PFEI, [range_lst[i][j][idx_R1_lst[i][j]]], [PFEI_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_PFEI, [range_lst[i][j][idx_R2_lst[i][j]]], [PFEI_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
savefig(plot_PFEI, joinpath(save_dir_sub, "PFEI_comparison.png"))

# PFEI_Zoom
plot_PFEI2 = plot(xlabel="Range (nmi)", ylabel="PFEI (J/J)", dpi=800)
plot!(plot_PFEI2, xlims=(0.0,2000),ylims=(0.65, 1.1))
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_PFEI2, range_lst[i][j], PFEI_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_PFEI2, [range_lst[i][j][idx_R1_lst[i][j]]], [PFEI_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_PFEI2, [range_lst[i][j][idx_R2_lst[i][j]]], [PFEI_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
savefig(plot_PFEI2, joinpath(save_dir_sub, "PFEI_comparison_zoom.png"))

# Empty weight
plot_mEmpty = plot(xlabel="Range (nmi)", ylabel="Empty Mass (Ton)", dpi=800)
plot!(plot_mEmpty, ylims=(0.0, massEmp_lst[1][1][1]+30.0))
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_mEmpty, range_lst[i][j], massEmp_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_mEmpty, [range_lst[i][j][idx_R1_lst[i][j]]], [massEmp_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_mEmpty, [range_lst[i][j][idx_R2_lst[i][j]]], [massEmp_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
savefig(plot_mEmpty, joinpath(save_dir_sub, "mass_empty_comparison.png"))

# Fuel volume
plot_volfuel = plot(xlabel="Range (nmi)", ylabel="Fuel Volume (m3)", dpi=800)
# plot!(plot_volfuel, ylims=(0.0, massEmp_lst[1][1][1]+30.0))
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_volfuel, range_lst[i][j], voluFuel_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_volfuel, [range_lst[i][j][idx_R1_lst[i][j]]], [voluFuel_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_volfuel, [range_lst[i][j][idx_R2_lst[i][j]]], [voluFuel_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
# for (i, keyword_cur) in enumerate(case_keywords)
#     for (j, des_range_cur) in enumerate(des_ranges)
#         plot!(plot_volfuel, range_lst[i][j], voluFuelMax_lst[i][j], marker=:none, lw=0.5, label=label=case_names[i]*"_$(round(Int,des_ranges[j]))")
#     end
# end
savefig(plot_volfuel, joinpath(save_dir_sub, "volume_fuel_comparison.png"))

# Takeeoff weight
plot_MTO = plot(xlabel="Range (nmi)", ylabel="Takeoff Mass (Ton)", dpi=800)
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_MTO, range_lst[i][j], massTO_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_MTO, [range_lst[i][j][idx_R1_lst[i][j]]], [massTO_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_MTO, [range_lst[i][j][idx_R2_lst[i][j]]], [massTO_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_MTO, range_lst[i][j], massTOMax_lst[i][j], marker=:none, lw=0.5, label=false)
    end
end
savefig(plot_MTO, joinpath(save_dir_sub, "takeoff_mass_comparison.png"))

# Empty weight
plot_mPay = plot(xlabel="Range (nmi)", ylabel="Payload Mass (Ton)", dpi=800)
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_mPay, range_lst[i][j], massPay_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_mPay, [range_lst[i][j][idx_R1_lst[i][j]]], [massPay_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_mPay, [range_lst[i][j][idx_R2_lst[i][j]]], [massPay_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
savefig(plot_mPay, joinpath(save_dir_sub, "mass_payload_comparison.png"))

# Lift-to-drag ratio
plot_LD = plot(xlabel="Range (nmi)", ylabel="Lift-to-drag at Cruise", dpi=800)
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_LD, range_lst[i][j], LD_cru_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_LD, [range_lst[i][j][idx_R1_lst[i][j]]], [LD_cru_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_LD, [range_lst[i][j][idx_R2_lst[i][j]]], [LD_cru_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
savefig(plot_LD, joinpath(save_dir_sub, "Lift_to_drag_comparison.png"))

# Engine total efficiency
plot_etaEng = plot(xlabel="Range (nmi)", ylabel="Engine Total Efficiency at Cruise", dpi=800)
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_etaEng, range_lst[i][j], eta_tot_cru_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_etaEng, [range_lst[i][j][idx_R1_lst[i][j]]], [eta_tot_cru_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_etaEng, [range_lst[i][j][idx_R2_lst[i][j]]], [eta_tot_cru_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
savefig(plot_etaEng, joinpath(save_dir_sub, "Engine_total_efficiency_comparison.png"))

# Heating value
plot_LHV = plot(xlabel="Range (nmi)", ylabel="Fuel Heating Value (J/kg)", dpi=800)
for (i, keyword_cur) in enumerate(case_keywords)
    for (j, des_range_cur) in enumerate(des_ranges)
        plot!(plot_LHV, range_lst[i][j], LHV_lst[i][j], marker=:cross, lw=2, label=case_names[i]*"_$(round(Int,des_ranges[j]))")
        # Mark down R1 R2
        scatter!(plot_LHV, [range_lst[i][j][idx_R1_lst[i][j]]], [LHV_lst[i][j][idx_R1_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(plot_LHV, [range_lst[i][j][idx_R2_lst[i][j]]], [LHV_lst[i][j][idx_R2_lst[i][j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    end
end
savefig(plot_LHV, joinpath(save_dir_sub, "Fuel_heating_value_comparison.png"))

# Percentage range reduction
plot_R1R2Red = plot(xlabel="Design Range (nmi)", ylabel="Change of Ranges WRT Baseline (%)", dpi=800)
plot!(plot_R1R2Red, des_range_R1R2, frac_change_R1 .* 100.0, marker=:cross, lw=2, label="R1")
plot!(plot_R1R2Red, des_range_R1R2, frac_change_R2 .* 100.0, marker=:cross, lw=2, label="R2")
savefig(plot_R1R2Red, joinpath(save_dir_sub, "R1R2Change.png"))

# PFEI Change
plot_PFEI_Change = plot(xlabel="Range (nmi)", ylabel="Change of PFEI WRT Baseline (%)", dpi=800)
plot!(plot_PFEI_Change,ylims=(-1, 25))
for (i, des_range_cur) in enumerate(des_ranges)
    plot!(plot_PFEI_Change, PFEI_ranges[i], PFEI_change[i] .* 100.0, marker=:cross, lw=2, label="$(round(Int,des_ranges[i]))")
    # Mark down R1 R2
    scatter!(plot_PFEI_Change, [PFEI_ranges[i][PFEI_R1R2_idx[i][1]]], [PFEI_change[i][PFEI_R1R2_idx[i][1]] * 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    scatter!(plot_PFEI_Change, [PFEI_ranges[i][PFEI_R1R2_idx[i][2]]], [PFEI_change[i][PFEI_R1R2_idx[i][2]] * 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
end
savefig(plot_PFEI_Change, joinpath(save_dir_sub, "PFEI_change.png"))

# Flight Energy Change
plot_Ene_Change = plot(xlabel="Range (nmi)", ylabel="Change of Flight Energy WRT Baseline (%)", dpi=800)
plot!(plot_Ene_Change,ylims=(-1, 25))
for (i, des_range_cur) in enumerate(des_ranges)
    plot!(plot_Ene_Change, PFEI_ranges[i], energy_flight_change[i] .* 100.0, marker=:cross, lw=2, label="$(round(Int,des_ranges[i]))")
    # Mark down R1 R2
    scatter!(plot_Ene_Change, [PFEI_ranges[i][PFEI_R1R2_idx[i][1]]], [energy_flight_change[i][PFEI_R1R2_idx[i][1]] * 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
    scatter!(plot_Ene_Change, [PFEI_ranges[i][PFEI_R1R2_idx[i][2]]], [energy_flight_change[i][PFEI_R1R2_idx[i][2]] * 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
end
savefig(plot_Ene_Change, joinpath(save_dir_sub, "Ene_change.png"))