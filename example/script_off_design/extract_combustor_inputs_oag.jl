"""
This script extract the combustor operating condtitions for oag off-design missions.
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(__TASOPTroot__,"../example/utilities_for_optimization/objective_factory.jl"))
off_design_specified! = ObjectiveFactory.off_design_specified!
include(joinpath(__TASOPTroot__,"../example/utilities_for_postprocessing/Extract.jl"))
using .Extract: read_oag,extract_combustion_inputs

#### Setup IO
# Input case names - OAG seat-capacity sweep
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
caseKey   = "Opti_Jet_NoACT_OAG_6Seats_TypeC_V4_"
idx_fuel_case      = 24
rho_fuel_case_kgm3 = 817.0
hvap_fuel_case_Jkg = 358694.0
pass_load_frac_off = 0.825 # Off-design payload load factor, matches opt_from_multi_warm_starts_para_oag.jl
miss_dir = joinpath(@__DIR__,"../ModelSaved/OAG_Data_2024/OAG_Data_2024.csv")
# Output folder name (within model_dir)
save_name     = "Combustor_$(caseKey)" #sub_folder will be created

#### Create save directory
save_dir_sub  = joinpath(model_dir,save_name)
mkpath(save_dir_sub)

#### Load OAG off-design mission data (ranges/weights per seat capacity)
miss_off_des = read_oag(miss_dir)
seat_cap_keys_all = sort(collect(keys(miss_off_des)))

#### Determine which seat capacities have a saved optimized model, per case
seat_caps_avail = Vector{Int}()
for sc in seat_cap_keys_all
    ac_dir = joinpath(model_dir,caseKey,caseKey*"$(sc).jld2")
    isfile(ac_dir) && push!(seat_caps_avail,sc)
end
println("  Available Seat Capacities: $(seat_caps_avail)")

#### Compute the oag missions engine and combustor performance and save
for (i, sc) in enumerate(seat_caps_avail)
    # Read in the case
    ac_dir = joinpath(model_dir,caseKey,caseKey*"$(sc).jld2")
    ac = quickload_aircraft(ac_dir)
    println("File, $(ac_dir), read successfully")

    # Run the OAG-weighted off-design missions (2nd mission) for this seat capacity
    ranges_off_nmi = miss_off_des[sc].ranges_nmi
    weights_off    = miss_off_des[sc].weights # This is the statistic weighting between missions instead of any physical weight
    wei_pay_off_N  = fill(ac.parg[igWpaymax]*pass_load_frac_off, length(ranges_off_nmi))
    
    # Test the off-design range one by one
    for idx_off in eachindex(ranges_off_nmi)
        ac_cur = deepcopy(ac)
        out = off_design_specified!(ac_cur, idx_fuel_case, rho_fuel_case_kgm3, hvap_fuel_case_Jkg, [ranges_off_nmi[idx_off]], [wei_pay_off_N[idx_off]]; mod_ac_inplace=true)
        length(out.wei_pay_N)<=0 && error("Offdesign point for $(ranges_off_nmi[idx_off]) nmi and $(sc) seats did not converge")
        # Extract the combustor operating condition
        extract_combustion_inputs(ac_cur,2,save_name*"$(sc)Seats_Miss$(idx_off)",save_dir_sub)
    end
end