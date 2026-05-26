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

function find_frac_change(results,idx_base,idx_targ,numDesRan,xSymb,ySymb)
    #results: case*desRan*symb*offdes
    #idx_base,idx_targ: cases
    #numDesRan: number of design range
    #xSymb,ySymb: symbols to find difference, Ex. x: range, y: PFEI
    x_common = []
    y_change = []
    for i in 1:numDesRan
        #### Extract the base and target case
        results_base = results[idx_base][i]
        results_targ = results[idx_targ][i]
        x_base = results_base[xSymb]
        x_targ = results_targ[xSymb]
        y_base = results_base[ySymb]
        y_targ = results_targ[ySymb]
        #### Filter out the common cases using x
        minx = max(minimum(x_base),minimum(x_targ))
        maxx = min(maximum(x_base),maximum(x_targ))
        msk_base = (x_base .>= minx) .& (x_base .<= maxx)
        msk_targ = (x_targ .>= minx) .& (x_targ .<= maxx)
        x_base = x_base[msk_base]
        x_targ = x_targ[msk_targ]
        y_base = y_base[msk_base]
        y_targ = y_targ[msk_targ]
        @assert (round.(Int, x_base) == round.(Int, x_targ)) "Base$(x_base)/target$(x_targ) x do not match after filtering"
        #### Calculate the change
        push!(x_common, x_base)
        push!(y_change, (y_targ .- y_base) ./ y_base)
    end
    return x_common,y_change
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
idx_base_PFEI = 2
idx_targ_PFEI = 3
# For Breguet range plot
flg_plot_Breguet = true
idx_constraints_Breguet = 1 #Index of case to obtain limiting parameters including fuel volume, maximum takeoff weight, and maximum payload weight

# Output directory
save_dir      = "ModelProcessed"
save_name     = "Compare_Retrofit" #sub_folder will be created
# Fields to read out
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

#### Use Breguet Range to calculate the PRD for checking
if flg_plot_Breguet
    const fields_breguet = (:range_Bre_nmi,:PFEI_Bre_JJ,:range_cru_Bre_nmi,:PFEI_cru_Bre_JJ,:range_fix_Bre_nmi,:PFEI_fix_Bre_JJ)
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
            #Cut off the missing part
            for f in fields_breguet
                resize!(results_Bre_cur[f], length(results_cur[:massPay_Ton]))
            end
        end
    end
    #### Use fixed parameters to compute Breguet range performance
    for i in eachindex(case_keywords)
        for j in eachindex(des_ranges)
            results_cur = results[i][j]
            results_fix = results[idx_base_PFEI][j]
            results_lim = results[idx_constraints_Breguet][j] #Case to extract aircraft limitations
            results_Bre_cur = results_Bre[i][j]
            idxOD_fix = argmin(results_fix[:PFEI_JJ]) #Index of the off-design range for the fixed case where parameter will be extracted
            global kend = length(results_cur[:massPay_Ton])
            for k in eachindex(results_cur[:massPay_Ton])
                outDict   = Bre_off_des(results_cur[:range_nmi][k] * 1852.0, results_cur[:massPay_Ton][k] * 1000.0 * gee;
                                        LD=results_fix[:LD_cru][idxOD_fix], eta=results_fix[:eta_tot_cru][idxOD_fix], LHV_Jkg=results_cur[:LHV_Jkg][k],
                                        wEmp_N=results_fix[:massEmp_Ton][idxOD_fix] * 1000.0 * gee, rhoFuel_kgm3=results_cur[:rhoFuel_kgm3][k], frac_rese=results_fix[:frac_rese][idxOD_fix],
                                        wTO_Max_N=maximum(results_lim[:massTO_Ton]) * 1000.0 * gee, wPay_Max_N=maximum(results_lim[:massPay_Ton]) * 1000.0 * gee, volFuel_Max_m3=maximum(results_lim[:voluFuel_m3]),
                                        gee=gee) #Due to cruise flight effciency use, it is believe the predicted missions requires less fuel and less TO weight.
                range_fix_Bre_nmi = outDict["range_m_out"] / 1852.0
                massPay_Bre_Ton = outDict["wPay_N_out"] / gee / 1000.0
                (abs(massPay_Bre_Ton-results_cur[:massPay_Ton][k]) < 0.001) || println("Warning: Breguet range leads to change in payload range, probably limited by flight efficiency")
                if range_fix_Bre_nmi > 0
                    results_Bre_cur[:range_fix_Bre_nmi][k] = range_fix_Bre_nmi
                    results_Bre_cur[:PFEI_fix_Bre_JJ][k] = outDict["PFEI_JJ_out"]
                else
                    global kend = k-1
                    println("Warning: Breguet range infeasible solution at $(results_cur[:range_nmi][k]) nmi")
                    break
                end
            end
            #Cut off the missing part just for this fixed parameter calculation
            for f in (:range_fix_Bre_nmi,:PFEI_fix_Bre_JJ)
                resize!(results_Bre_cur[f], kend)
            end
        end
    end
end

#### Calculate fractional change of PFEI
R1Idx_PFEI = []
R2Idx_PFEI = []
ranges_PFEI, PFEI_change = find_frac_change(results,idx_base_PFEI,idx_targ_PFEI,length(des_ranges),:range_nmi,:PFEI_JJ)
_, EneFli_change = find_frac_change(results,idx_base_PFEI,idx_targ_PFEI,length(des_ranges),:range_nmi,:EneFli_J)
for i in eachindex(des_ranges)
    #### Use the target R1 and R2 ranges to remap the current R1 R2 indices
    @assert !ismissing(R1[idx_targ_PFEI,i]) "Target case for PFEI change calculation has to have a R1 R2 indicator computed"
    push!(R1Idx_PFEI, argmin(abs.(ranges_PFEI[i] .- R1[idx_targ_PFEI,i])))
    push!(R2Idx_PFEI, argmin(abs.(ranges_PFEI[i] .- R2[idx_targ_PFEI,i])))
end

#### Calculate fractional change of PFEI from Breguet Range Analysis
if flg_plot_Breguet
    ranges_PFEI_Bre_Fix, PFEI_change_Bre_Fix = find_frac_change(results_Bre,idx_base_PFEI,idx_targ_PFEI,length(des_ranges),:range_fix_Bre_nmi,:PFEI_fix_Bre_JJ)
    ranges_PFEI_Bre_Var, PFEI_change_Bre_Var = find_frac_change(results_Bre,idx_base_PFEI,idx_targ_PFEI,length(des_ranges),:range_Bre_nmi,:PFEI_Bre_JJ)
end

#### Plotting - PFEI and energy
linestyles = repeat([:solid, :dash, :dot, :dashdot, :dashdotdot],1000)
linecolors = repeat([:blue, :red, :green, :orange, :purple],1000)
# Large PFEI
p1_1 = plot(xlabel="Off-design Range (nmi)", ylabel="PFEI (J/J)", dpi=800, yscale=:log10)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p1_1, results_cur[:range_nmi], results_cur[:PFEI_JJ], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p1_1, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:PFEI_JJ][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p1_1, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:PFEI_JJ][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
        if flg_plot_Breguet
            results_Bre_cur = results_Bre[i][j]
            plot!(p1_1, results_Bre_cur[:range_Bre_nmi], results_Bre_cur[:PFEI_Bre_JJ], marker=:none, color=linecolors[il], lw=0.75, linestyle=linestyles[il], label="$(case_names[i]) (Breguet), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        end
    end
end
savefig(p1_1, joinpath(save_dir_sub, "PFEI_Large.png"))
# plot!(pExample, xlims=(100,1600),ylims=(0.6, 1.0))

# Zoom in plot PFEI
p1_2 = plot(xlabel="Off-design Range (nmi)", ylabel="PFEI (J/J)", dpi=800)
global il = 0
plot!(p1_2, xlims=(0.0,2000),ylims=(0.65, 1.1))
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p1_2, results_cur[:range_nmi], results_cur[:PFEI_JJ], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p1_2, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:PFEI_JJ][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p1_2, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:PFEI_JJ][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p1_2, joinpath(save_dir_sub, "PFEI_Small.png"))

# Plot cruise PFEI comparison with Breguet range prediction
p1_3 = plot(xlabel="Off-design Cruise Range (nmi)", ylabel="Cruise only PFEI (J/J)", dpi=800)
global il = 0
plot!(p1_3, xlims=(0.0,2000),ylims=(0.0, 1.5))
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        global il += 1
        results_cur = results[i][j]
        plot!(p1_3, results_cur[:range_cru_m] / 1852.0, results_cur[:PFEI_cru_JJ], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if flg_plot_Breguet
            results_Bre_cur = results_Bre[i][j]
            plot!(p1_3, results_Bre_cur[:range_cru_Bre_nmi], results_Bre_cur[:PFEI_cru_Bre_JJ], marker=:none, color=:black, lw=0.75, linestyle=linestyles[il], label="$(case_names[i]) (Breguet), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        end
    end
end
savefig(p1_3, joinpath(save_dir_sub, "PFEI_Cruise.png"))

# Large Flight Energy
p1_4 = plot(xlabel="Off-design Range (nmi)", ylabel="Flight Energy (J)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p1_4, results_cur[:range_nmi], results_cur[:EneFli_J], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p1_4, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:EneFli_J][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p1_4, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:EneFli_J][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p1_4, joinpath(save_dir_sub, "Flight_Energy.png"))

#### Plotting - Payload Range
p2_1 = plot(xlabel="Off-design Range (nmi)", ylabel="Payload Mass (Ton)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p2_1, results_cur[:range_nmi], results_cur[:massPay_Ton], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p2_1, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:massPay_Ton][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p2_1, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:massPay_Ton][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p2_1, joinpath(save_dir_sub, "Payload.png"))

#### Plotting - Other Weight and Volumes (sanity check)
# Emtpy weight
p3_1 = plot(xlabel="Off-design Range (nmi)", ylabel="Empty Mass (Ton)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p3_1, results_cur[:range_nmi], results_cur[:massEmp_Ton], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p3_1, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:massEmp_Ton][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p3_1, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:massEmp_Ton][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p3_1, joinpath(save_dir_sub, "Empty.png"))

# Takeoff weight
p3_2 = plot(xlabel="Off-design Range (nmi)", ylabel="Takeoff Mass (Ton)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p3_2, results_cur[:range_nmi], results_cur[:massTO_Ton], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p3_2, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:massTO_Ton][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p3_2, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:massTO_Ton][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p3_2, joinpath(save_dir_sub, "TakeoffWeight.png"))

#Fuel volume
p3_3 = plot(xlabel="Off-design Range (nmi)", ylabel="Fuel Volume (m3)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p3_3, results_cur[:range_nmi], results_cur[:voluFuel_m3], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p3_3, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:voluFuel_m3][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p3_3, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:voluFuel_m3][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p3_3, joinpath(save_dir_sub, "FuelVolume.png"))

#### Plotting - Fuel Efficiency Related
# Flight Efficiency
p4_1 = plot(xlabel="Off-design Range (nmi)", ylabel="Cruise Lift-to-drag Ratio", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p4_1, results_cur[:range_nmi], results_cur[:LD_cru], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p4_1, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:LD_cru][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p4_1, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:LD_cru][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p4_1, joinpath(save_dir_sub, "LD_Ratio.png"))

# Engine Total Efficiency
p4_2 = plot(xlabel="Off-design Range (nmi)", ylabel="Cruise Engine Total Efficiency", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p4_2, results_cur[:range_nmi], results_cur[:eta_tot_cru], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p4_2, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:eta_tot_cru][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p4_2, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:eta_tot_cru][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p4_2, joinpath(save_dir_sub, "eta_tot_Engine.png"))

# Fuel Heating Value (sanity check)
p4_3 = plot(xlabel="Off-design Range (nmi)", ylabel="Fuel Heating Value (J/kg)", dpi=800)
global il = 0
for i in eachindex(case_keywords)
    for j in eachindex(des_ranges)
        results_cur = results[i][j]
        global il += 1
        plot!(p4_3, results_cur[:range_nmi], results_cur[:LHV_Jkg], marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="$(case_names[i]), R₂: $(round(Int,R2[idx_base_R1R2,j])) nmi")
        if !ismissing(R1Idx[i,j])
            # Mark down R1 R2 location
            scatter!(p4_3, [results_cur[:range_nmi][R1Idx[i,j]]], [results_cur[:LHV_Jkg][R1Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
            scatter!(p4_3, [results_cur[:range_nmi][R2Idx[i,j]]], [results_cur[:LHV_Jkg][R2Idx[i,j]]], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
        end
    end
end
savefig(p4_3, joinpath(save_dir_sub, "Fuel_heating_value.png"))

#### Plotting - Change between cases
# Change of R1 and R2
p5_1 = plot(xlabel="Baseline R2 Range (nmi)", ylabel="Range Reduction from Retrofitting (%)", dpi=800)
scatter!(p5_1, R2[idx_base_R1R2,:], -frac_change_R1 .* 100.0, marker=:cross, ms=4, msw=2.5, label="R1 Reduction")
scatter!(p5_1, R2[idx_base_R1R2,:], -frac_change_R2 .* 100.0, marker=:cross, ms=4, msw=2.5, label="R2 Reduction")
savefig(p5_1, joinpath(save_dir_sub, "R1_R2_Change.png"))

# PFEI Increase
p5_2 = plot(xlabel="Off-design Range (nmi)", ylabel="PFEI Increase from Retrofitting(%)", dpi=800)
# plot!(p5_2,ylims=(-1, 25))
global il = 0
for i in eachindex(des_ranges)
    global il += 1
    plot!(p5_2, ranges_PFEI[i], PFEI_change[i] .* 100.0, marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="Variable Parameters R₂: $(round(Int,R2[idx_base_R1R2,i])) nmi")
    if !ismissing(R1Idx_PFEI[i])
        # Mark down R1 R2 location
        scatter!(p5_2, [ranges_PFEI[i][R1Idx_PFEI[i]]], [PFEI_change[i][R1Idx_PFEI[i]] .* 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(p5_2, [ranges_PFEI[i][R2Idx_PFEI[i]]], [PFEI_change[i][R2Idx_PFEI[i]] .* 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
    end
    if flg_plot_Breguet
        plot!(p5_2, ranges_PFEI_Bre_Fix[i], PFEI_change_Bre_Fix[i] .* 100.0, marker=:diamond, ms=3, msw=0, color=linecolors[il], lw=0.75, linestyle=linestyles[il], label="Fixed Parameters (Breguet), R₂: $(round(Int,R2[idx_base_R1R2,i])) nmi")
        plot!(p5_2, ranges_PFEI_Bre_Var[i], PFEI_change_Bre_Var[i] .* 100.0, marker=:utriangle, ms=3, msw=0, color=linecolors[il], lw=0.75, linestyle=linestyles[il], label="Variable Parameters (Breguet), R₂: $(round(Int,R2[idx_base_R1R2,i])) nmi")
    end
end
savefig(p5_2, joinpath(save_dir_sub, "PFEI_Change.png"))

# Flight Energy Increase
p5_3 = plot(xlabel="Off-design Range (nmi)", ylabel="Flight Energy Increase from Retrofitting(%)", dpi=800)
global il = 0
for i in eachindex(des_ranges)
    global il += 1
    plot!(p5_3, ranges_PFEI[i], EneFli_change[i] .* 100.0, marker=:cross, color=linecolors[il], lw=2, linestyle=linestyles[il], label="R₂: $(round(Int,R2[idx_base_R1R2,i])) nmi")
    if !ismissing(R1Idx_PFEI[i])
        # Mark down R1 R2 location
        scatter!(p5_3, [ranges_PFEI[i][R1Idx_PFEI[i]]], [EneFli_change[i][R1Idx_PFEI[i]] .* 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R1
        scatter!(p5_3, [ranges_PFEI[i][R2Idx_PFEI[i]]], [EneFli_change[i][R2Idx_PFEI[i]] .* 100.0], marker=:cross, ms=4, msw=2.5, mc=:black, msc=:black, label=false) #Mark R2
    end
end
savefig(p5_3, joinpath(save_dir_sub, "Flight_Energy_Change.png"))