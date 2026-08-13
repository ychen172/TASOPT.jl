"""
This script compares the OAG seat-capacity-swept optimization campaigns
produced by opt_from_multi_warm_starts_para_oag.jl.
Off-design performance is compared via an OAG-weighted PFEI, computed over the same
route-frequency-weighted mission set (OAG_Data_2024.csv) used to optimize each case,
re-run here on the saved (2nd-mission) off-design missions.
Picker the better (Low off-design averaged PFEI) or feasible case from the provided flight models across the seat capacities
Assume the models for filtering all have the same fuel type
Assume all model folders contain the same set of files regardless whether some cases are converged(Have ac jld2) or not (Only have optpara and other history)
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../utilities_for_optimization/objective_factory.jl"))
off_design_specified! = ObjectiveFactory.off_design_specified!
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel_compact!, init_results_2Layers, plot_cases_specified, read_oag

#### Setup IO
# Input case names - OAG seat-capacity sweep
model_dir  = "../ModelSaved"
caseKeys   = ["Opti_Jet_NoACT_OAG_6Seats_TypeC_","Opti_Jet_NoACT_OAG_6Seats_TypeC_"]
# Off-design fuel properties, aligned with caseKeys (must match what each campaign was optimized/run with)
idx_fuel_case      = 24        # Eth: 32, Jet: 24
rho_fuel_case_kgm3 = 817.0     # Eth: 789.0, Jet: 817.0 kg/m3
hvap_fuel_case_Jkg = 358694.0  # Eth: 918187.9, Jet: 358694.0 J/kg
pass_load_frac_off = 0.825 # Off-design payload load factor, matches opt_from_multi_warm_starts_para_oag.jl
# OAG route-frequency mission data (off-design ranges/weights, keyed by seat_capacity)
miss_dir = joinpath(@__DIR__,"../ModelSaved/OAG_Data_2024/OAG_Data_2024.csv")
# Output directory
save_name     = "Opti_Jet_NoACT_OAG_6Seats_TypeC_V2" #sub_folder will be created

#### Create save directory
save_dir  = joinpath(model_dir,save_name)
mkpath(save_dir)

#### Load OAG off-design mission data (ranges/weights per seat capacity)
miss_off_des = read_oag(miss_dir)
seat_cap_keys_all = sort(collect(keys(miss_off_des)))

#### Determine which seat capacities have a saved optimized model, per case
# (Jet and Ethanol need not match exactly - e.g. Jet succeeded at 80 seats while Ethanol did not)
seat_caps_avail = Vector{Vector{Int}}()
for caseKey in caseKeys
    avail = Int[]
    for sc in seat_cap_keys_all
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(sc).jld2")
        isfile(ac_dir) && push!(avail,sc)
    end
    push!(seat_caps_avail,avail)
end
println("Available seat capacities per case:")
for (caseName,avail) in zip(caseNames,seat_caps_avail)
    println("  $(caseName): $(avail)")
end

#### Initialization
oag_weighted_PFEI = [fill(NaN,length(seat_caps_avail[j])) for j in eachindex(caseKeys)] #[J/J]

#### Extract design-mission data, and compute OAG route-frequency-weighted off-design PFEI
for (j, caseKey) in enumerate(caseKeys)
    for (i, sc) in enumerate(seat_caps_avail[j])
        # Read in the case
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(sc).jld2")
        ac = quickload_aircraft(ac_dir)
        println("File, $(ac_dir), read successfully")

        # Run the OAG-weighted off-design missions (2nd mission) for this seat capacity
        ranges_off_nmi = miss_off_des[sc].ranges_nmi
        weights_off    = miss_off_des[sc].weights
        wei_pay_off_N  = fill(ac.parg[igWpaymax]*pass_load_frac_off, length(ranges_off_nmi))
        out = off_design_specified!(ac, idx_fuel_case, rho_fuel_case_kgm3, hvap_fuel_case_Jkg,
                                    ranges_off_nmi, wei_pay_off_N;
                                    mod_ac_inplace=false, itermax=150, constraints=[:WPay,:MWTO,:VolFuel], save_model=false)

        # Average PFEI = weighted sum(PFEI*payload*range) / weighted sum(payload*range), using the OAG weight
        # for each mission and the actual off-design (2nd mission) payload/range/PFEI returned above.
        # off_design_specified! silently drops missions that fail to size/converge, preserving relative order,
        # so match each returned mission back to its OAG weight via an order-preserving two-pointer walk.
        k_ref = 1
        n_matched = 0
        num_energy   = 0.0 # sum(w_i * PFEI_i * payload_i * range_i)
        den_payrange = 0.0 # sum(w_i * payload_i * range_i)
        for idx_out in eachindex(out.ran_nmi)
            while k_ref <= length(ranges_off_nmi) && !isapprox(out.ran_nmi[idx_out], ranges_off_nmi[k_ref]; atol=1e-3, rtol=1e-6)
                k_ref += 1
            end
            k_ref <= length(ranges_off_nmi) || error("Could not match off-design output range $(out.ran_nmi[idx_out]) back to the OAG mission list for $(caseKey) at seat capacity $(sc)")
            w_i = weights_off[k_ref]
            payrange_i = out.wei_pay_N[idx_out]*out.ran_nmi[idx_out] # off-design (2nd) mission's actual payload*range
            num_energy   += w_i*out.PFEI_JJ[idx_out]*payrange_i
            den_payrange += w_i*payrange_i
            k_ref += 1
            n_matched += 1
        end
        n_fail = length(ranges_off_nmi) - n_matched

        if den_payrange > 0.0
            oag_weighted_PFEI[j][i] = num_energy/den_payrange
        else
            @warn "All OAG off-design missions failed for $(caseKey) at seat capacity $(sc); OAG-weighted PFEI left as missing"
        end
        n_fail>0 && println("  $(n_fail)/$(length(ranges_off_nmi)) OAG off-design missions failed for $(caseKey) at seat capacity $(sc)")
        println("Case $(caseNames[j]) at seat capacity $(sc): OAG-weighted PFEI = $(oag_weighted_PFEI[j][i]) J/J")
    end
end

#### Plotting: Design (R1) mission comparisons vs seat capacity
xdata_design = [Int64.(avail) for avail in seat_caps_avail]

#### Screen out the best performance at each range
lookup = [Dict(zip(xdata_design[i],oag_weighted_PFEI[i])) for i in eachindex(xdata_design)]
idx_col_best = fill(missing,length(seat_cap_keys_all))
for i in eachindex(idx_col_best)
    PFEI_Collect = fill(Inf,length(lookup))
    for j in eachindex(lookup)
        try
            PFEI_Collect[j] = lookup[j][seat_cap_keys_all[i]]
        catch
        end
    end
    if all(PFEI_Collect .== Inf)
        idx_col_best[i] = 1
    else
        idx_col_best[i] = argmin(PFEI_Collect)
    end
end

#### Save the copy the corresponding best model from each case into the new folder
for (idx,idx_best) in enumerate(idx_col_best)
    sc = round(Int, seat_cap_keys_all[idx])
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(sc)_design_constraints.csv")
    dst = joinpath(save_dir,
                   save_name * "$(sc)_design_constraints.csv")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(sc)_global_bounds.csv")
    dst = joinpath(save_dir,
                   save_name * "$(sc)_global_bounds.csv")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(sc)_mission_requirements.csv")
    dst = joinpath(save_dir,
                   save_name * "$(sc)_mission_requirements.csv")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(sc)_optimization_history.jld2")
    dst = joinpath(save_dir,
                   save_name * "$(sc)_optimization_history.jld2")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(sc)_optimized_parameters.csv")
    dst = joinpath(save_dir,
                   save_name * "$(sc)_optimized_parameters.csv")
    cp(src, dst; force=true)
    #
    src = joinpath(model_dir,
                   caseKeys[idx_best],
                   caseKeys[idx_best] * "$(sc)_OptLog.txt")
    dst = joinpath(save_dir,
                   save_name * "$(sc)_OptLog.txt")
    cp(src, dst; force=true)
    #
    try
        src = joinpath(model_dir,
                    caseKeys[idx_best],
                    caseKeys[idx_best] * "$(sc).jld2")
        dst = joinpath(save_dir,
                    save_name * "$(sc).jld2")
        cp(src, dst; force=true)
    catch
    end
end