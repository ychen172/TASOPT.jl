"""
This script compares the PFEI between ethanol and jet fuel aircraft from R1 sizing process
    - Compare PFEI between ethanol and jet fuel aircraft the sizing mission
    - Compare PFEI difference from TASOPT vs Breguet range (pure fuel energy density effect)
"""

using TASOPT
include(__TASOPTindices__)
using Plots
using Glob
include(joinpath(@__DIR__,"../UtilityFunctions/Breguet_range_solve_offdes.jl"))
using .Breguet: Bre_off_des

#### Setup IO
# Input case names - Retrofit
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
save_dir   = joinpath(__TASOPTroot__,"../example/OutputAndPlots")
jet_key    = "Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_"
eth_key    = "Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_"
jet_Nam    = "Jet Fuel"
eth_Nam    = "Ethanol"
save_prefi = "R1_Size_PFEI"
ref_fuel   = :jet #Which fuel's airframe (L/D, engine efficiency, empty weight, reserve fraction) is
                   #held fixed as the reference for the Breguet range sweep below (:jet or :eth)
range_ref  = 2000 #nmi, the reference range where reference airframe performance and each fuel's
                   #own LHV/density are extracted for the Breguet calculation

#### Create save directory if not exist
mkpath(save_dir)

#### Helper: discover available design ranges for a given case key (glob + regex on saved .jld2 files,
#### matching "<key><range>.jld2" only, excluding "..._optimization_history.jld2")
function available_ranges(model_dir, key)
    case_dir = joinpath(model_dir, key)
    pattern = Regex("^" * key * "(\\d+)\\.jld2\$")
    rans_avail = Int[]
    for f in glob(key*"*.jld2", case_dir)
        m = match(pattern, basename(f))
        if m !== nothing
            push!(rans_avail, parse(Int, m.captures[1]))
        end
    end
    sort!(rans_avail)
    return rans_avail
end

rans_jet = available_ranges(model_dir, jet_key)
rans_eth = available_ranges(model_dir, eth_key)
println("$(jet_key): found $(length(rans_jet)) available design ranges")
println("$(eth_key): found $(length(rans_eth)) available design ranges")

#### Only compare ranges that actually converged for BOTH fuels
ranges_common = sort(intersect(rans_jet, rans_eth))
println("Common ranges available: $(length(ranges_common))")

#### Loop through the common ranges once, loading both aircraft at each range and extracting PFEI
PFEI_jet = Float64[]
PFEI_eth = Float64[]
for ran_i in ranges_common
    ac_jet = quickload_aircraft(joinpath(model_dir, jet_key, jet_key*"$(ran_i).jld2"))
    ac_eth = quickload_aircraft(joinpath(model_dir, eth_key, eth_key*"$(ran_i).jld2"))
    println("Range $(ran_i) nmi: Jet and Ethanol loaded")

    push!(PFEI_jet, ac_jet.parm[imPFEI,1])
    push!(PFEI_eth, ac_eth.parm[imPFEI,1])
end

#### Actual (TASOPT-optimized) PFEI penalty of Ethanol vs Jet, at matching ranges
PFEIEthPenalty_Opt = 100.0 .* (PFEI_eth .- PFEI_jet) ./ PFEI_jet

#### Collect reference values for the simplified Breguet range calculation: pull L/D, engine total
#### efficiency, empty weight, and reserve fuel fraction from a SINGLE reference aircraft (ref_fuel
#### at range_ref), to be held FIXED across the whole Breguet sweep below. This isolates the pure
#### fuel-chemistry (LHV + density) effect from the optimizer's compensating airframe redesign
#### (bigger wing, heavier structure, etc.) that the TASOPT-optimized PFEI above already includes.
ac_jet_ref = quickload_aircraft(joinpath(model_dir, jet_key, jet_key*"$(range_ref).jld2"))
ac_eth_ref = quickload_aircraft(joinpath(model_dir, eth_key, eth_key*"$(range_ref).jld2"))
ac_ref = ref_fuel == :jet ? ac_jet_ref : ac_eth_ref

LD_ref        = ac_ref.para[iaCL, ipcruise1, 1] / ac_ref.para[iaCD, ipcruise1, 1]
eta_ref       = (1.0/(ac_ref.pare[ieTSFC, ipcruise1, 1]/gee)) *
                ((cos(ac_ref.para[iagamV, ipcruise1, 1])*ac_ref.pare[ieu0, ipcruise1, 1])/ac_ref.pare[iehfuel, ipcruise1, 1])
wEmp_ref      = ac_ref.parm[imWTO,1] - ac_ref.parm[imWfuel,1] - ac_ref.parm[imWpay,1]
frac_rese_ref = ac_ref.parg[igfreserve]

#### Collect fuel information for Breguet range prediction (LHV and density are the only things
#### that differ between the two Bre_off_des calls below - everything else above is held constant)
LHV_jet     = ac_jet_ref.pare[iehfuel, ipcruise1, 1]
LHV_eth     = ac_eth_ref.pare[iehfuel, ipcruise1, 1]
rhofuel_jet = ac_jet_ref.parg[igrhofuel]
rhofuel_eth = ac_eth_ref.parg[igrhofuel]

#### Loop through the reference fuel's own ranges/payloads for Breguet range estimation: for each
#### point, ask "if the SAME reference airframe (fixed L/D, eta, empty weight) simply burned Jet
#### fuel vs Ethanol via the Breguet equation, how would PFEI compare?" - no optimizer redesign involved.
ref_key = ref_fuel == :jet ? jet_key : eth_key
PFEIEthPenalty_Breguet = zeros(length(ranges_common))
for (i, ran_i) in enumerate(ranges_common)
    ac_i = quickload_aircraft(joinpath(model_dir, ref_key, ref_key*"$(ran_i).jld2"))
    range_m = ac_i.parm[imRange,1]
    wPay_N  = ac_i.parm[imWpay,1]

    outJet = Bre_off_des(range_m, wPay_N;
                        LD = LD_ref, eta = eta_ref, LHV_Jkg = LHV_jet, wEmp_N = wEmp_ref,
                        rhoFuel_kgm3 = rhofuel_jet, frac_rese = frac_rese_ref,
                        wTO_Max_N=4e10, wPay_Max_N=1e10, volFuel_Max_m3=1e10/1e5) #Not considering physical limitations
    outEth = Bre_off_des(range_m, wPay_N;
                        LD = LD_ref, eta = eta_ref, LHV_Jkg = LHV_eth, wEmp_N = wEmp_ref,
                        rhoFuel_kgm3 = rhofuel_eth, frac_rese = frac_rese_ref,
                        wTO_Max_N=4e10, wPay_Max_N=1e10, volFuel_Max_m3=1e10/1e5) #Not considering physical limitations

    PFEIEthPenalty_Breguet[i] = 100.0*(outEth["PFEI_JJ_out"]-outJet["PFEI_JJ_out"])/outJet["PFEI_JJ_out"]
end

#### Plots
#### Single shared color sequence, same as empty_weight_R1S.jl/aero_per_R1S.jl's BASE_COLORS
const BASE_COLORS = [:teal, :royalblue]

#### PFEI comparison: actual TASOPT-optimized PFEI, Jet vs Ethanol
p1 = plot(ranges_common, PFEI_jet, label=jet_Nam, color=BASE_COLORS[1], lw=2.5, marker=:circle, markersize=3,
          xlabel="Design/R1 Range [nmi]", ylabel="Flight PFEI [J/J]", dpi=800,
          legend=:topleft, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18,
          left_margin=16Plots.mm, bottom_margin=12Plots.mm)
plot!(p1, ranges_common, PFEI_eth, label=eth_Nam, color=BASE_COLORS[2], lw=2.5, marker=:circle, markersize=3)
savefig(p1, joinpath(save_dir, "$(save_prefi)_PFEI.png"))

#### PFEI difference comparison: actual (optimized) penalty vs the pure-fuel-chemistry-only
#### Breguet-based penalty, to show how much of Ethanol's PFEI penalty is pure chemistry vs the
#### optimizer's compensating airframe redesign (the rest of this conversation's main finding)
p2 = plot(ranges_common, PFEIEthPenalty_Opt, label="Optimized Parameters", color=BASE_COLORS[1], lw=2.5, marker=:circle, markersize=3,
          xlabel="Design/R1 Range [nmi]", ylabel="Ethanol Flight PFEI Penalty [%]", dpi=800,
          legend=:topleft, size=(1400, 900), titlefontsize=21, guidefontsize=24, tickfontsize=18, legendfontsize=18,
          left_margin=16Plots.mm, bottom_margin=12Plots.mm)
plot!(p2, ranges_common, PFEIEthPenalty_Breguet, label="Constant Efficiencies (Breguet)", color=BASE_COLORS[2], lw=2.5, marker=:circle, markersize=3)
savefig(p2, joinpath(save_dir, "$(save_prefi)_EthPFEIPenalty.png"))
