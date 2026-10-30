"""
This script extract the combustor operating condtitions for oag off-design missions.
"""

using TASOPT
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"../example/UtilityFunctions/objective_factory.jl"))
off_design_specified! = ObjectiveFactory.off_design_specified!
include(joinpath(__TASOPTroot__,"../example/UtilityFunctions/Extract.jl"))
using .Extract: read_oag,extract_combustion_inputs
using Glob

#### Setup IO
# Input case names - OAG seat-capacity sweep
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
caseFolder = "Opti_Du31_NoACT_OAG_MI_TailCorr_V3_"
caseKey = r"^Opti_Du31_NoACT_OAG_Ml_.*_TailCorr_V3_(\d+)\.jld2$"
idx_fuel_case      = 322431
rho_fuel_case_kgm3 = 805.368
hvap_fuel_case_Jkg = 586408.0
pass_load_frac_off = 0.825 # Off-design payload load factor, matches opt_from_multi_warm_starts_para_oag.jl
pass_load_frac_tail = 0.850 # The farthest (tail) mission payload fraction overwrite, matches opt_from_multi_warm_starts_para_oag_Rerun.jl
miss_dir = joinpath(model_dir,"OffDesignMissions_50_300_300_Tail.csv")
# Output folder name (within model_dir)
save_name     = "CombFLI_$(caseFolder)" #sub_folder will be created

#### Create save directory
save_dir_sub  = joinpath(model_dir,save_name)
mkpath(save_dir_sub)

#### Load OAG off-design mission data (ranges/weights per seat capacity)
miss_off_des = read_oag(miss_dir)

#### Determine which seat capacities have a saved optimized model in caseFolder
sc = Int[]
ac_name = String[]
model_dir_sub = joinpath(model_dir,caseFolder)
matches = glob([caseKey],model_dir_sub)
for match_cur in matches
    match_cur = basename(match_cur)
    pattern = match(caseKey,match_cur)
    sc_cur = parse(Int, pattern[1])
    push!(sc, sc_cur)
    push!(ac_name, match_cur)
end
idx_sort = sortperm(sc)
permute!(sc,idx_sort)
permute!(ac_name,idx_sort)
println("Available Seat Capacities: $(sc)")

#### Compute the oag missions engine and combustor performance and save
for (sc_cur,ac_name_cur) in zip(sc,ac_name)
    # Read in the case
    ac_dir = joinpath(model_dir_sub,ac_name_cur)
    ac = quickload_aircraft(ac_dir)
    println("File, $(ac_dir), read successfully")

    # Run the OAG-weighted off-design missions (2nd mission) for this seat capacity
    ranges_off_nmi = miss_off_des[sc_cur].ranges_nmi
    wei_pay_off_N  = fill(ac.parg[igWpaymax]*pass_load_frac_off, length(ranges_off_nmi))
    idx_max_range  = argmax(ranges_off_nmi)
    wei_pay_off_N[idx_max_range] = ac.parg[igWpaymax]*pass_load_frac_tail

    # Test the off-design range one by one
    for idx_off in eachindex(ranges_off_nmi)
        ac_cur = deepcopy(ac)
        out = off_design_specified!(ac_cur, idx_fuel_case, rho_fuel_case_kgm3, hvap_fuel_case_Jkg, [ranges_off_nmi[idx_off]], [wei_pay_off_N[idx_off]]; mod_ac_inplace=true)
        length(out.wei_pay_N)<=0 && error("Offdesign point for $(ranges_off_nmi[idx_off]) nmi and $(sc_cur) seats did not converge")
        # Extract the combustor operating condition
        extract_combustion_inputs(ac_cur,2,save_name*"$(sc_cur)Seats_Miss$(idx_off)",save_dir_sub)
    end
end