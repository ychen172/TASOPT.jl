"""
This script compare the off-design performance between cases in terms of PFEI
and flight range
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"Breguet_range_solve_offdes.jl"))
using .Breguet: Bre_off_des

# Other constants
epsR1R2 = 0.005

#### Helpers Functions
# Interpolation function for threshold range determination
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

# Design parameters output
function init_results(case_keywords, des_ranges, offdes_ranges, fields_to_read)
    # case_keywords, des_ranges, offdes_ranges: Vector
    # fields_to_read: tuple of symbols
    noff = length(offdes_ranges)
    makevec() = fill!(Vector{Union{Missing,Float64}}(undef, noff), missing)
    Dict(
        kw => Dict(
            d => Dict(
                f => makevec() for f in fields_to_read
            ) for d in eachindex(des_ranges)
        ) for kw in eachindex(case_keywords)
    )
end

#### Setup IO
# Input case names
case_keywords = ["off_designJet", "jetfuel_match_payload", "jetfuel_to_ethanolJet"]
case_names    = ["Baseline", "Matched Baseline", "Retrofit"]
model_dir     = "ModelSaved"
des_ranges    = [3000] #float.(collect(300:100:3000)) #design range to compare (Has to be integer (No 0.1 nmi)) (nmi) Make sure all cases have these design ranges
offdes_ranges = float.(collect(300:100:8000)) #Off-design ranges to search through (Has to be integer (No 0.1 nmi)) (can be wider than what are available)
# For R1 and R2 calculation
idx_R1R2Skip  = [2] #Case to skip R1 R2 determination
idx_base_R1R2 = 1 #Index of the base case
idx_targ_R1R2 = 3 #Index of the target case
# For PFEI comparison
idx_base_PFEI = 1
idx_targ_PFEI = 3
# For Breguet range plot
flg_plot_Breguet = true

# Output directory
save_dir      = "ModelProcessed"
save_name     = "Compare_Retrofit" #sub_folder will be created
# Fields to read out
# const fields_to_read = (:range_nmi,   :PFEI_JJ, :massEmp_Ton, :voluFuel_m3, :voluFuelMax_m3, :massTO_Ton, :massTOMax_Ton,
#                         :massPay_Ton, :LD_cru,  :eta_tot_cru, :LHV_Jkg)
const fields_to_read = (:range_nmi,   :PFEI_JJ, :massEmp_Ton, :voluFuel_m3, :massTO_Ton,
                        :massPay_Ton, :LD_cru,  :eta_tot_cru, :LHV_Jkg, :EneFli_J, :frac_rese, :rhoFuel_kgm3, :PFEI_cru_JJ, :range_cru_m)
#### Save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

#### Initialization
results = init_results(case_keywords, des_ranges, offdes_ranges, fields_to_read)

#### Extract data for each scenerio
for (i, keyword_cur) in enumerate(case_keywords)
    model_dir_sub = joinpath(model_dir, keyword_cur)
    # Extract data for each design case
    for (j, des_ranges_cur) in enumerate(des_ranges)
        dkey = round(Int, des_ranges_cur)
        model_dir_sub_sub = joinpath(model_dir_sub, keyword_cur*"$(dkey)")
        results_cur = results[i][j] #field * off-design
        # extract for each off-design case
        k = 0
        for offdes_range_cur in offdes_ranges
            odkey = round(Int, offdes_range_cur) #off-design range key
            model_name_to_read = joinpath(model_dir_sub_sub, keyword_cur*"$(dkey)_$(odkey).jld2")
            #### reading in the model
            if !isfile(model_name_to_read) #check if the file exist
                continue
            end
            println("Read: $(model_name_to_read)")
            # Get the aircraft model
            ac_cur = quickload_aircraft(model_name_to_read)
            
            #### overall performance data
            range_cur = ac_cur.parm[imRange,2] / 1852.0 #(nmi)
            PFEI_cur = ac_cur.parm[imPFEI, 2] #(J/J)
            
            #### mass data
            massTO = ac_cur.parm[imWTO,2]/gee/1000.0 #Takeoff mass (Ton)
            massFuelTot = ac_cur.parm[imWfuel,2]/gee/1000.0 #Fuel mass (Ton)(Include reserved and burned)
            massPayload = ac_cur.parm[imWpay, 2]/gee/1000.0 #(Ton)
            massEmpty = massTO - massFuelTot - massPayload #(Ton) empty weight
            # massTOMax = ac_cur.parg[igWMTO] / gee / 1000.0 #Maximum takeoff mass (Ton)
            
            #### fuel volume data
            rhoFuel = ac_cur.parg[igrhofuel] #kg/m3
            volFuel = massFuelTot * 1000.0 / rhoFuel #m3
            # volFuelMax = ac_cur.parg[igWfmax] / gee / rhoFuel #m3 (The design mission fuel mass might be different from the maximum fuel mass with off-design fuel density)
            
            #### flight performance data
            LD_cruise = 0.5 * (ac_cur.para[iaCL, ipcruise1, 2]/ac_cur.para[iaCD, ipcruise1, 2] + 
                                ac_cur.para[iaCL, ipcruise2, 2]/ac_cur.para[iaCD, ipcruise2, 2]) #Averaged cruise lift-to-drag ratio
            LHV_cruise = 0.5 * (ac_cur.pare[iehfuel, ipcruise1, 2] + ac_cur.pare[iehfuel, ipcruise2, 2]) #Averaged cruise heating value (J/kg) (Include vaporization already)
            TSFC_cruise = 0.5 * (ac_cur.pare[ieTSFC, ipcruise1, 2] + ac_cur.pare[ieTSFC, ipcruise2, 2]) / gee #Averaged cruise thrust specfic heat consumption (kg/s/N)
            vel_cruise = 0.5 * (cos(ac_cur.para[iagamV, ipcruise1, 2]) * ac_cur.pare[ieu0, ipcruise1, 2] + 
                                cos(ac_cur.para[iagamV, ipcruise2, 2]) * ac_cur.pare[ieu0, ipcruise2, 2]) #Averaged cruise horizontal velocity (m/s)
            eta_total_cruise = (1.0/TSFC_cruise)*(vel_cruise/LHV_cruise) #total cruise engine efficiency

            #### derived data
            energy_flight = PFEI_cur * massPayload * (1000.0 * gee * 1852.0) * range_cur #(J) total flight energy 

            #### parameters for Breguet range
            frac_rese = ac_cur.parg[igfreserve] #W_reserveFuel / W_fuelburned
            PFEI_cru = (LHV_cruise*(ac_cur.para[iafracW, ipcruise1, 2]-ac_cur.para[iafracW, ipcruise2, 2])*ac_cur.parg[igWMTO]/gee)/
                       (ac_cur.parm[imWpay, 2] * ac_cur.parm[imRange, 2]) #Cruise Only PFEI (J/J)
            range_cru = ac_cur.para[iaRange, ipcruise2, 2] - ac_cur.para[iaRange, ipcruise1, 2] #(m)
            @assert range_cru>0.0 "Find a negative cruise range, likely range to short"

            ## store
            k += 1
            results_cur[:range_nmi][k] = range_cur
            results_cur[:PFEI_JJ][k] = PFEI_cur
            results_cur[:massEmp_Ton][k] = massEmpty
            results_cur[:voluFuel_m3][k] = volFuel
            results_cur[:massTO_Ton][k] = massTO
            results_cur[:massPay_Ton][k] = massPayload
            results_cur[:LD_cru][k] = LD_cruise
            results_cur[:eta_tot_cru][k] = eta_total_cruise
            results_cur[:LHV_Jkg][k] = LHV_cruise
            results_cur[:EneFli_J][k] = energy_flight
            results_cur[:frac_rese][k] = frac_rese
            results_cur[:rhoFuel_kgm3][k] = rhoFuel
            results_cur[:PFEI_cru_JJ][k] = PFEI_cru
            results_cur[:range_cru_m][k] = range_cru
        end
        #### Trim the trailing missing
        for f in fields_to_read
            resize!(results_cur[f], k)
        end
    end
end

#### Identify the R1 and R2 for each case and each design range
R1 = fill!(Matrix{Union{Missing,Float64}}(undef,length(case_keywords),length(des_ranges)), missing)
R2 = fill!(Matrix{Union{Missing,Float64}}(undef,length(case_keywords),length(des_ranges)), missing)
R1Idx = fill!(Matrix{Union{Missing,Int}}(undef,length(case_keywords),length(des_ranges)), missing)
R2Idx = fill!(Matrix{Union{Missing,Int}}(undef,length(case_keywords),length(des_ranges)), missing)
for i in eachindex(case_keywords)
    if i in idx_R1R2Skip
        continue
    end
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        # R1
        massPay_Max_Val, massPay_Max_Idx = findmax(results_cur[:massPay_Ton]) #assume your off-design mission reach maximum payload (if not, the first point would be picked)
        idx_R1_approx = findnext(massPay -> massPay < massPay_Max_Val*(1.0-epsR1R2), results_cur[:massPay_Ton], massPay_Max_Idx) - 1 #approximated R1 index for plot
        R1_interp = interp_range_at_threshold(results_cur[:range_nmi], results_cur[:massPay_Ton], massPay_Max_Idx, epsR1R2; side=:right) #inteprolated R1
        # R2
        volFuel_Max_Val, volFuel_Max_Idx = findmax(results_cur[:voluFuel_m3]) #assume your off-design reach maximum fuel volume (threortically always for payload range diagram)
        idx_R2_approx = findprev(volFuel -> volFuel < volFuel_Max_Val*(1.0-epsR1R2), results_cur[:voluFuel_m3], volFuel_Max_Idx) + 1
        R2_interp = interp_range_at_threshold(results_cur[:range_nmi], results_cur[:voluFuel_m3], volFuel_Max_Idx, epsR1R2; side=:left) #inteprolated R2
        # Store R1 and R2
        R1[i,j] = R1_interp
        R2[i,j] = R2_interp
        R1Idx[i,j] = idx_R1_approx
        R2Idx[i,j] = idx_R2_approx
    end
end
# Calculate fractional change of R1 and R2
frac_change_R1 = (R1[idx_targ_R1R2,:] .- R1[idx_base_R1R2,:]) ./ R1[idx_base_R1R2,:]
frac_change_R2 = (R2[idx_targ_R1R2,:] .- R2[idx_base_R1R2,:]) ./ R2[idx_base_R1R2,:]

#### Calculate fractional change of PFEI
ranges_PFEI = []
PFEI_change = []
EneFli_change = []
R1Idx_PFEI = []
R2Idx_PFEI = []
for i in eachindex(des_ranges)
    #### Extract the base and target case
    results_base = results[idx_base_PFEI][i]
    results_targ = results[idx_targ_PFEI][i]
    # Extract and PFEI and range
    range_base = results_base[:range_nmi]
    range_targ = results_targ[:range_nmi]
    PFEI_base  = results_base[:PFEI_JJ]
    PFEI_targ  = results_targ[:PFEI_JJ]
    EneFli_base = results_base[:EneFli_J]
    EneFli_targ = results_targ[:EneFli_J]
    # Filter the PFEI and range for common subset
    min_range  = max(minimum(range_base),minimum(range_targ)) #common range bound (assume same spacing)
    max_range  = min(maximum(range_base),maximum(range_targ))
    msk_base   = (range_base .>= min_range) .& (range_base .<= max_range)
    msk_targ   = (range_targ .>= min_range) .& (range_targ .<= max_range)
    range_base = range_base[msk_base] 
    PFEI_base  = PFEI_base[msk_base] 
    EneFli_base = EneFli_base[msk_base]
    range_targ = range_targ[msk_targ] 
    PFEI_targ  = PFEI_targ[msk_targ]
    EneFli_targ = EneFli_targ[msk_targ]
    @assert round.(Int, range_base) == round.(Int, range_targ) "Base/target ranges do not match after filtering"

    #### Calculate the changes
    push!(ranges_PFEI, range_base)
    push!(PFEI_change, (PFEI_targ .- PFEI_base) ./ PFEI_base)
    push!(EneFli_change, (EneFli_targ .- EneFli_base) ./ EneFli_base)
    
    #### Use the target R1 and R2 ranges to remap the current R1 R2 indices
    @assert !ismissing(R1[idx_targ_PFEI,i]) "Target case for PFEI change calculation has to have a R1 R2 indicator computed"
    push!(R1Idx_PFEI, argmin(abs.(range_base .- R1[idx_targ_PFEI,i])))
    push!(R2Idx_PFEI, argmin(abs.(range_base .- R2[idx_targ_PFEI,i])))
end

#### Use Breguet Range to calculate the PRD for checking
const fields_breguet = (:range_Bre_nmi,:PFEI_Bre_JJ,:range_cru_Bre_nmi,:PFEI_cru_Bre_JJ)
results_Bre = init_results(case_keywords, des_ranges, offdes_ranges, fields_breguet)
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        #### Use the parameter for each case compute a Breguet range counter part
        results_cur = results[i][j]
        results_Bre_cur = results_Bre[i][j]
        # Full mission PFEI
        for k in eachindex(results_cur[:massPay_Ton])
            outDict   = Bre_off_des(results_cur[:range_nmi][k] * 1852.0, results_cur[:massPay_Ton][k] * 1000.0 * gee;
                                    LD=results_cur[:LD_cru][k], eta=results_cur[:eta_tot_cru][k], LHV_Jkg=results_cur[:LHV_Jkg][k],
                                    wEmp_N=results_cur[:massEmp_Ton][k] * 1000.0 * gee, rhoFuel_kgm3=results_cur[:rhoFuel_kgm3][k], frac_rese=results_cur[:frac_rese][k],
                                    wTO_Max_N=results_cur[:massTO_Ton][k] * 1000.0 * gee, wPay_Max_N=results_cur[:massPay_Ton][k] * 1000.0 * gee, volFuel_Max_m3=results_cur[:voluFuel_m3][k],
                                    gee=gee) #Due to cruise flight effciency use, it is believe the predicted missions requires less fuel and less TO weight.
            range_Bre_nmi = outDict["range_m_out"] / 1852.0
            massPay_Bre_Ton = outDict["wPay_N_out"] / gee / 1000.0
            @assert abs(range_Bre_nmi-results_cur[:range_nmi][k])<1 "Breguet computed range is different from requested"
            @assert abs(massPay_Bre_Ton-results_cur[:massPay_Ton][k])<0.001 "Breguet computed payload mass is more than 1 kg different"
            results_Bre_cur[:range_Bre_nmi][k] = range_Bre_nmi
            results_Bre_cur[:PFEI_Bre_JJ][k] = outDict["PFEI_JJ_out"]
        end
        # Cruise only PFEI (Normalize still by full mission flight range)
        for k in eachindex(results_cur[:massPay_Ton])
            outDict   = Bre_off_des(results_cur[:range_cru_m][k], results_cur[:massPay_Ton][k] * 1000.0 * gee;
                        LD=results_cur[:LD_cru][k], eta=results_cur[:eta_tot_cru][k], LHV_Jkg=results_cur[:LHV_Jkg][k],
                        wEmp_N=results_cur[:massEmp_Ton][k] * 1000.0 * gee, rhoFuel_kgm3=results_cur[:rhoFuel_kgm3][k], frac_rese=results_cur[:frac_rese][k],
                        wTO_Max_N=results_cur[:massTO_Ton][k] * 1000.0 * gee, wPay_Max_N=results_cur[:massPay_Ton][k] * 1000.0 * gee, volFuel_Max_m3=results_cur[:voluFuel_m3][k],
                        gee=gee) #Due to cruise flight effciency use, it is believe the predicted missions requires less fuel and less TO weight.
            massPay_Bre_Ton = outDict["wPay_N_out"] / gee / 1000.0
            @assert abs(outDict["range_m_out"]-results_cur[:range_cru_m][k])<1 "Breguet computed range is different from requested"
            @assert abs(massPay_Bre_Ton-results_cur[:massPay_Ton][k])<0.001 "Breguet computed payload mass is more than 1 kg different"
            PFEI_cru_Bre_JJ = outDict["PFEI_JJ_out"]*results_cur[:range_cru_m][k]/(results_cur[:range_nmi][k] * 1852.0) #Renormalized by total flight range
            results_Bre_cur[:range_cru_Bre_nmi][k] = outDict["range_m_out"] / 1852.0
            results_Bre_cur[:PFEI_cru_Bre_JJ][k] = PFEI_cru_Bre_JJ
        end
        for f in fields_breguet #Cut off the missing part
            resize!(results_Bre_cur[f], length(results_cur[:massPay_Ton]))
        end
    end
end

#### Plotting - PFEI and energy
linestyles = repeat([:solid, :dash, :dot, :dashdot, :dashdotdot],1000)
linecolors = repeat([:blue, :red, :green, :orange, :purple],1000)
# Large PFEI
p1 = plot(xlabel="Off-design Range (nmi)", ylabel="PFEI (J/J)", dpi=800, yscale=:log10)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p1, results_cur[:range_nmi], results_cur[:PFEI_JJ], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p1, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:PFEI_JJ][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p1, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:PFEI_JJ][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
        if flg_plot_Breguet
            results_Bre_cur = results_Bre[i][j]
            plot!(p1, results_Bre_cur[:range_Bre_nmi], results_Bre_cur[:PFEI_Bre_JJ], marker=:none, color=linecolors[il], lw=0.75, linestyle=linestyles[il], label="$(case_names[i]) (Breguet), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        end
    end
end
savefig(p1, joinpath(save_dir_sub, "PFEI_Large.png"))
# plot!(pExample, xlims=(100,1600),ylims=(0.6, 1.0))

# Zoom in plot PFEI
p2 = plot(xlabel="Off-design Range (nmi)", ylabel="PFEI (J/J)", dpi=800)
global il = 0
plot!(p2, xlims=(0.0,2000),ylims=(0.65, 1.1))
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p2, results_cur[:range_nmi], results_cur[:PFEI_JJ], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p2, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:PFEI_JJ][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p2, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:PFEI_JJ][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p2, joinpath(save_dir_sub, "PFEI_Small.png"))
if flg_plot_Breguet
    # Plot cruise PFEI comparison with Breguet range prediction
    p2_4 = plot(xlabel="Off-design Cruise Range (nmi)", ylabel="Cruise only PFEI (J/J)", dpi=800)
    global il = 0
    plot!(p2_4, xlims=(0.0,2000),ylims=(0.0, 1.5))
    for i in eachindex(case_keywords)
        for j in eachindex(des_ranges)
            results_cur = results[i][j]
            results_Bre_cur = results_Bre[i][j]
            global il += 1
            plot!(p2_4, results_cur[:range_cru_m] / 1852.0, results_cur[:PFEI_cru_JJ], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
            plot!(p2_4, results_Bre_cur[:range_cru_Bre_nmi], results_Bre_cur[:PFEI_cru_Bre_JJ], marker=:none, color=:black, lw=0.75, linestyle=linestyles[il], label="$(case_names[i]) (Breguet), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        end
    end
    savefig(p2_4, joinpath(save_dir_sub, "PFEI_Cruise.png"))
end

# Large Flight Energy
p2_5 = plot(xlabel="Off-design Range (nmi)", ylabel="Flight Energy (J)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p2_5, results_cur[:range_nmi], results_cur[:EneFli_J], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p2_5, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:EneFli_J][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p2_5, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:EneFli_J][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p2_5, joinpath(save_dir_sub, "Flight_Energy.png"))

#### Plotting - Payload Range
p3 = plot(xlabel="Off-design Range (nmi)", ylabel="Payload Mass (Ton)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p3, results_cur[:range_nmi], results_cur[:massPay_Ton], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p3, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:massPay_Ton][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p3, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:massPay_Ton][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p3, joinpath(save_dir_sub, "Payload.png"))

#### Plotting - Other Weight and Volumes (sanity check)
# Emtpy weight
p4 = plot(xlabel="Off-design Range (nmi)", ylabel="Empty Mass (Ton)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p4, results_cur[:range_nmi], results_cur[:massEmp_Ton], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p4, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:massEmp_Ton][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p4, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:massEmp_Ton][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p4, joinpath(save_dir_sub, "Empty.png"))

# Takeoff weight
p5 = plot(xlabel="Off-design Range (nmi)", ylabel="Takeoff Mass (Ton)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p5, results_cur[:range_nmi], results_cur[:massTO_Ton], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p5, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:massTO_Ton][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p5, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:massTO_Ton][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p5, joinpath(save_dir_sub, "TakeoffWeight.png"))

#Fuel volume
p6 = plot(xlabel="Off-design Range (nmi)", ylabel="Fuel Volume (m3)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p6, results_cur[:range_nmi], results_cur[:voluFuel_m3], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p6, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:voluFuel_m3][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p6, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:voluFuel_m3][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p6, joinpath(save_dir_sub, "FuelVolume.png"))

#### Plotting - Breguet Range Related
# Flight Efficiency
p7 = plot(xlabel="Off-design Range (nmi)", ylabel="Cruise Lift-to-drag Ratio", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p7, results_cur[:range_nmi], results_cur[:LD_cru], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p7, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:LD_cru][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p7, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:LD_cru][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p7, joinpath(save_dir_sub, "LD_Ratio.png"))

# Engine Total Efficiency
p8 = plot(xlabel="Off-design Range (nmi)", ylabel="Cruise Engine Total Efficiency", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p8, results_cur[:range_nmi], results_cur[:eta_tot_cru], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p8, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:eta_tot_cru][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p8, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:eta_tot_cru][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p8, joinpath(save_dir_sub, "eta_tot_Engine.png"))

# Fuel Heating Value (sanity check)
p9 = plot(xlabel="Off-design Range (nmi)", ylabel="Fuel Heating Value (kg/J)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p9, results_cur[:range_nmi], results_cur[:LHV_Jkg], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p9, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:LHV_Jkg][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p9, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:LHV_Jkg][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p9, joinpath(save_dir_sub, "Fuel_heating_value.png"))

#### Plotting - Change between cases
# Change of R1 and R2
p10 = plot(xlabel="Baseline R2 Range (nmi)", ylabel="Range Reduction from Retrofitting (%)", dpi=800)
scatter!(p10, R2[idx_base_R1R2,:], -frac_change_R1 .* 100.0, marker=:cross, ms=4, msw=2.5, label="R1 Reduction")
scatter!(p10, R2[idx_base_R1R2,:], -frac_change_R2 .* 100.0, marker=:cross, ms=4, msw=2.5, label="R2 Reduction")
savefig(p10, joinpath(save_dir_sub, "R1_R2_Change.png"))

# PFEI Increase
p11 = plot(xlabel="Off-design Range (nmi)", ylabel="PFEI Increase from Retrofitting(%)", dpi=800)
global il = 0
for i in eachindex(des_ranges)
    global il += 1
    plot!(p11, ranges_PFEI[i], PFEI_change[i] .* 100.0, marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="R₂: $(round(Int,R2[idx_base_R1R2,i])) nmi")
    if !ismissing(R1Idx_PFEI[i])
        # Mark down R1 R2 location
        scatter!(p11, [ranges_PFEI[i][R1Idx_PFEI[i]]], [PFEI_change[i][R1Idx_PFEI[i]] .* 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(p11, [ranges_PFEI[i][R2Idx_PFEI[i]]], [PFEI_change[i][R2Idx_PFEI[i]] .* 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
    end
end
savefig(p11, joinpath(save_dir_sub, "PFEI_Change.png"))

# Flight Energy Increase
p12 = plot(xlabel="Off-design Range (nmi)", ylabel="Flight Energy Increase from Retrofitting(%)", dpi=800)
global il = 0
for i in eachindex(des_ranges)
    global il += 1
    plot!(p12, ranges_PFEI[i], EneFli_change[i] .* 100.0, marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="R₂: $(round(Int,R2[idx_base_R1R2,i])) nmi")
    if !ismissing(R1Idx_PFEI[i])
        # Mark down R1 R2 location
        scatter!(p12, [ranges_PFEI[i][R1Idx_PFEI[i]]], [EneFli_change[i][R1Idx_PFEI[i]] .* 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(p12, [ranges_PFEI[i][R2Idx_PFEI[i]]], [EneFli_change[i][R2Idx_PFEI[i]] .* 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
    end
end
savefig(p12, joinpath(save_dir_sub, "Flight_Energy_Change.png"))