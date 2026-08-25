"""
This script compares the OAG seat-capacity-swept optimization campaigns (Jet vs Ethanol)
produced by opt_from_multi_warm_starts_para_oag.jl.

Design/R1-mission metrics are compared directly across seat capacity.
Off-design performance is compared via an OAG-weighted PFEI, computed over the same
route-frequency-weighted mission set (OAG_Data_2024.csv) used to optimize each case,
re-run here on the saved (2nd-mission) off-design missions.
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
caseKeys   = ["Opti_Jet_NoACT_OAG_6Seats_TypeC_V2_","Opti_Eth_NoACT_OAG_6Seats_TypeC_V3P_","Opti_Jet_NoACT_OAG_6Seats_TypeE_V1_","Opti_Eth_NoACT_OAG_6Seats_TypeE_V1_"]
caseNames  = ["Jet Fuel"                           ,"Ethanol"                             ,"Jet Fuel with Type E Span"          ,"Ethanol with Type E Span"           ]
# Off-design fuel properties, aligned with caseKeys (must match what each campaign was optimized/run with)
idx_fuel_case      = [24       ,32       ,24       ,32       ] # Jet, Eth ,32       ,24
rho_fuel_case_kgm3 = [817.0    ,789.0    ,817.0    ,789.0    ] # kg/m3.   ,789.0    ,817.0
hvap_fuel_case_Jkg = [358694.0 ,918187.9 ,358694.0 ,918187.9 ] # J/kg     ,918187.9 ,358694.0
pass_load_frac_off = 0.825 # Off-design payload load factor, matches opt_from_multi_warm_starts_para_oag.jl
constraints        = [[:WPay,:MWTO,:VolFuel],[:WPay,:MWTO,:VolFuel],[:WPay,:MWTO,:VolFuel],[:WPay,:MWTO,:VolFuel]] #Constraints for off-design
# OAG route-frequency mission data (off-design ranges/weights, keyed by seat_capacity)
miss_dir = joinpath(@__DIR__,"../ModelSaved/OAG_Data_2024/OAG_Data_2024.csv")
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "OAG_Jet_Eth_WingSpanE" #sub_folder will be created
iter_max      = 150 #max iteration for off-design calculation
# Fields to read out for the design (R1) mission
const fields = [:(parm[imRange,1]),:(parm[imPFEI,1]),:(parm[imVfuel,1]),:(parg[igVfmax]),
                :(wing.layout.span),:(wing.outboard.λ),:(wing.layout.ηs),:(wing.layout.AR),
                :(para[iaCL,ipcruise1,1]),:(para[iaCL,ipcruise2,1]),
                :(para[iaCD,ipcruise1,1]),:(para[iaCD,ipcruise2,1]),
                :(pare[iehfuel,ipcruise1,1]),:(pare[iehfuel,ipcruise2,1]),
                :(pare[ieTSFC,ipcruise1,1]),:(pare[ieTSFC,ipcruise2,1]),
                :(para[iagamV,ipcruise1,1]),:(para[iagamV,ipcruise2,1]),
                :(pare[ieu0,ipcruise1,1]),:(pare[ieu0,ipcruise2,1]),
                :(parm[imWTO,1]),:(parm[imWfuel,1]),:(parm[imWpay,1]),:(parg[igrhofuel]),:(parg[igfreserve]),
                :(wing.layout.S)]



#### Create save directory
save_dir_sub  = joinpath(save_dir,save_name)
mkpath(save_dir_sub)

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
dataset = [init_results_2Layers(length(seat_caps_avail[j]), fields) for j in eachindex(caseKeys)] #[dataset[Expr][:],...]
oag_weighted_PFEI = [fill(NaN,length(seat_caps_avail[j])) for j in eachindex(caseKeys)] #[J/J]

#### Extract design-mission data, and compute OAG route-frequency-weighted off-design PFEI
for (j, caseKey) in enumerate(caseKeys)
    for (i, sc) in enumerate(seat_caps_avail[j])
        # Read in the case
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(sc).jld2")
        ac = quickload_aircraft(ac_dir)
        println("File, $(ac_dir), read successfully")

        # Extract design (R1, 1st mission) data
        extract_acModel_compact!(ac, dataset[j], i)

        # Run the OAG-weighted off-design missions (2nd mission) for this seat capacity
        ranges_off_nmi = miss_off_des[sc].ranges_nmi
        weights_off    = miss_off_des[sc].weights
        wei_pay_off_N  = fill(ac.parg[igWpaymax]*pass_load_frac_off, length(ranges_off_nmi))
        out = off_design_specified!(ac, idx_fuel_case[j], rho_fuel_case_kgm3[j], hvap_fuel_case_Jkg[j],
                                    ranges_off_nmi, wei_pay_off_N;
                                    mod_ac_inplace=false, itermax=iter_max, constraints=constraints[j], save_model=false)

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
xdata_design = [Float64.(avail) for avail in seat_caps_avail]

plot_cases_specified("Seat Capacity", "Design/R1 Flight PFEI [J/J]",
                     xdata_design,
                     [d[:(parm[imPFEI,1])] for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"Design_PFEI.png"))
#
plot_cases_specified("Seat Capacity", "Fuel Level [%]",
                     xdata_design,
                     [d[:(parm[imVfuel,1])] .* 100.0 ./ d[:(parg[igVfmax])] for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"Fuel_Level.png"))
#
plot_cases_specified("Seat Capacity", "Wingspan [m]",
                     xdata_design,
                     [d[:(wing.layout.span)] for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"Wingspan.png"))
#
plot_cases_specified("Seat Capacity", "Taper Ratio of Outer Wing Panel",
                     xdata_design,
                     [d[:(wing.outboard.λ)] for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"OuterTaperWing.png"))
#
plot_cases_specified("Seat Capacity", "Wing Span Break Location [%]",
                     xdata_design,
                     [d[:(wing.layout.ηs)] .* 100.0 for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"SpanBreakWing.png"))
#
plot_cases_specified("Seat Capacity", "Wing Aspect Ratio",
                     xdata_design,
                     [d[:(wing.layout.AR)]  for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"AspectRatioWing.png"))
#
plot_cases_specified("Seat Capacity", "Lift-to-drag Ratio at Start of Cruise",
                     xdata_design,
                     [d[:(para[iaCL,ipcruise1,1])] ./ d[:(para[iaCD,ipcruise1,1])]  for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"LD_StartOfCruise.png"))
#
plot_cases_specified("Seat Capacity", "Lift-to-drag Ratio at End of Cruise",
                     xdata_design,
                     [d[:(para[iaCL,ipcruise2,1])] ./ d[:(para[iaCD,ipcruise2,1])]  for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"LD_EndOfCruise.png"))
#
plot_cases_specified("Seat Capacity", "Engine Total Efficiency at Start of Cruise [%]",
                     xdata_design,
                     [100.0 .* (1.0 ./ (d[:(pare[ieTSFC,ipcruise1,1])] ./ gee)) .* ((cos.(d[:(para[iagamV,ipcruise1,1])]) .* d[:(pare[ieu0,ipcruise1,1])]) ./ d[:(pare[iehfuel,ipcruise1,1])])  for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"TotEffEngine_StartOfCruise.png"))
#
plot_cases_specified("Seat Capacity", "Engine Total Efficiency at End of Cruise [%]",
                     xdata_design,
                     [100.0 .* (1.0 ./ (d[:(pare[ieTSFC,ipcruise2,1])] ./ gee)) .* ((cos.(d[:(para[iagamV,ipcruise2,1])]) .* d[:(pare[ieu0,ipcruise2,1])]) ./ d[:(pare[iehfuel,ipcruise2,1])])  for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"TotEffEngine_EndOfCruise.png"))
#
plot_cases_specified("Seat Capacity", "Aircraft Empty Weight [Ton]",
                     xdata_design,
                     [(d[:(parm[imWTO,1])] .- d[:(parm[imWfuel,1])] .- d[:(parm[imWpay,1])]) ./ gee ./ 1000.0  for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"EmptyWeight.png"))
#
plot_cases_specified("Seat Capacity", "Fuel Weight Fraction of Total Takeoff Weight [%]",
                     xdata_design,
                     [100.0 .* (d[:(parm[imWfuel,1])] ./ d[:(parm[imWTO,1])]) for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"FuelWeightFraction.png"))
#
plot_cases_specified("Seat Capacity", "Wing Surface Area [m²]",
                     xdata_design,
                     [d[:(wing.layout.S)] for d in dataset],
                     caseNames,
                     joinpath(save_dir_sub,"WingSurfaceArea.png"))
#

#### Plotting: OAG route-frequency-weighted off-design PFEI comparison
plot_cases_specified("Seat Capacity", "OAG-Weighted Off-Design PFEI [J/J]",
                     xdata_design,
                     oag_weighted_PFEI,
                     caseNames,
                     joinpath(save_dir_sub,"OAG_Weighted_PFEI.png"))
#

# #### Ethanol OAG-weighted PFEI penalty vs Jet, restricted to seat capacities both cases have
# iRJ = 1 #Index of the jet case
# iRE = 2 #Index of the ethanol case
# common_caps = sort(collect(intersect(Set(seat_caps_avail[iRJ]), Set(seat_caps_avail[iRE]))))
# PFEIEthPenalty_OAG = Float64[]
# for sc in common_caps
#     iJ = findfirst(==(sc), seat_caps_avail[iRJ])
#     iE = findfirst(==(sc), seat_caps_avail[iRE])
#     push!(PFEIEthPenalty_OAG, 100.0*(oag_weighted_PFEI[iRE][iE]-oag_weighted_PFEI[iRJ][iJ])/oag_weighted_PFEI[iRJ][iJ])
# end

# plot_cases_specified("Seat Capacity", "Ethanol OAG-Weighted PFEI Penalty [%]",
#                      [Float64.(common_caps)],
#                      [PFEIEthPenalty_OAG],
#                      [nothing],
#                      joinpath(save_dir_sub,"OAG_Weighted_PFEI_Penalty.png"))
# #
