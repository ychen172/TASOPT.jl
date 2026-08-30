"""
This script performs a single full-aircraft R1Size re-optimization (geometry + engine cycle),
warm-started from the original Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_2400.jld2 (same base airframe as
V1 -- NOT V1's own reopt result), with the 8 technology parameters fixed at the
test_engine_opt_LEAP1B_hard_ratio_M78_210S.jl double-loop calibration's result (TechRun#190,
Mach 0.78, 210-seat exit limit, hard-fixed 3-stage-LPC/10-stage-HPC ratio in that calibration --
this script itself keeps the engine cycle free/split, see below).

The 5 engine-cycle Parameters (pif/pilc/pihc/BPR/Tt4) use that same calibration's own converged
engine-cycle values (bestSol_eng) as their initial guess for the local search, instead of
extracting a starting value from the base airframe (which still carries its pre-reopt,
pre-tech-level cycle state -- a poor guess once the technology level has been swapped out from
under it). Since the calibration itself hard-fixed the LPC/HPC split, its single `per_stage`
value was expanded back into separate pilc/pihc initial guesses via `per_stage^3`/`per_stage^10`
(equal pressure ratio per stage, matching the calibration's numStageLC=3/numStageHC=10) --
this script's own search still treats pilc/pihc as independent free variables from there.
All other Parameters (geometry) still warm-start from the base airframe as in V1. Single aircraft
model, no job array.
"""

using TASOPT
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils","sensitivity.jl"))
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
using .ObjectiveFactory: Constraint, Parameter, OptHistory, Requirement, optimize_singlePt_PFEI!, optimizer_wrapper_global_local, size_aircraft_w_param!,
                         save_vec_struct_csv, load_csv_constraints, load_csv_parameters, load_csv_requirements, save_jld2, load_jld2, OffDesMission
const min_fuse_radius = TASOPT.structures.find_minimum_radius_for_seats_per_row
const upd_fuse_pax! = TASOPT.structures.update_fuse_for_pax!
const success_statuses = ObjectiveFactory.success_statuses

#### Optimization parameters
# Warm-start model (also the prefix for the aircraft model) -- the original base airframe, same as V1
par_path_base_prefix = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_/Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_")
# Path to save the model from optimization
save_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/")
save_key = "Opti_Jet_NoACT_V4_FreSpl_Max8_210S_M078_Glo"
# Optimization configuration parameters
flag_skip_global = false #Switch on to skip the global search if confident that warm start model can converge
max_iter_sizing = 150 # Maximum iterations for TASOPT sizing
optimizers = [:GN_CRS2_LM, :LN_NELDERMEAD] # Optimizer choice. [Global,Local]
max_num_iter_opt = [30000, 800, 10000] # Maximum number of optimization steps [Global, Local Coarse, Local Fine]
max_num_round_loc = [120, 30] #Maximum number of adaptive bounds refinements rounds for local search [Local Coarse, Local Fine]
span_glo_to_loc = 0.25 # span_local_search/span_global_search
# Single design range for this run
range_opti_nmi = 2621.0 #[nmi]
# Setup an entry (value does not matter) for test ratio to be sweep later in optimization loop
mis_opt = Requirement[]
push!(mis_opt, Requirement(:(parm[imRange,1]), range_opti_nmi*1852.0)) #[m]
push!(mis_opt, Requirement(:(options.ifuel), Int(24))) #Eth: 32       , Jet: 24
push!(mis_opt, Requirement(:(parg[igrhofuel]), 817.0)) #Eth: 789.0    , Jet: 817.0 #kg/m3
push!(mis_opt, Requirement(:(pare[iehvap, :, 1]), 358694.0)) #Eth: 918187.9 , Jet: 358694.0 #J/kg
push!(mis_opt, Requirement(:(pare[iehvapcombustor, :, 1]), 358694.0)) #Eth: 918187.9 , Jet: 358694.0 #J/kg
push!(mis_opt, Requirement(:(options.has_ACT_fuel), false)) # Whether to allow additional center fuel tank (ACT)
push!(mis_opt, Requirement(:(options.compensate_ACT), false)) # Whether to increase the aircraft length to accmondate for the cargo space taken by ACT
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_vol), 1.00)) # Volumetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_wei), 1.00)) # Gravimetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(para[iaMach, ipclimbn:ipdescent1, 1]), 0.78)) # Cruise Mach number
push!(mis_opt, Requirement(:(vtail.opt_sizing), TailSizing.OEI)) # Vertical tail sizing mode (Engine out or fixed volume)
push!(mis_opt, Requirement(:(parg[igCLveout]), 0.47)) # Vertical tail CL at engine out condition
push!(mis_opt, Requirement(:(htail.opt_sizing), TailSizing.CLmaxFwdCG)) # Horizontal tail sizing mode (Forward CG or fixed volume)
push!(mis_opt, Requirement(:(htail.CL_max_fwd_CG), -0.7)) # Horizontal tail CL at forward CG condition
push!(mis_opt, Requirement(:(parg[igxeng]), 1e10)) # [m] Engine axial location (Though fixed value for each optimization, to be)
push!(mis_opt, Requirement(:(parg[igdxeng2wbox]), 1e10)) # [m] Same as above   (set by the default model before each optimization)
# Fixed technology level -- test_engine_opt_LEAP1B_hard_ratio_M78_210S.jl calibration result
# (TechRun#190, penalty 7.021458306154994, Mach 0.78 / 210-seat exit limit), broadcast across
# every mission point (`:`), matching UpdAcTecLvl!'s own broadcast convention.
push!(mis_opt, Requirement(:(pare[iepib, :, 1]), 0.9571536034218671)) #Combustor pressure ratio
push!(mis_opt, Requirement(:(pare[ieepolf, :, 1]), 0.9331577273491901)) #Fan Poly Eff
push!(mis_opt, Requirement(:(pare[ieepollc, :, 1]), 0.913697684362504)) #LPC Poly Eff
push!(mis_opt, Requirement(:(pare[ieepolhc, :, 1]), 0.9184150776226279)) #HPC Poly Eff
push!(mis_opt, Requirement(:(pare[ieepolht, :, 1]), 0.9430533930681784)) #HPT Poly Eff
push!(mis_opt, Requirement(:(pare[ieepollt, :, 1]), 0.9366599150846535)) #LPT Poly Eff
push!(mis_opt, Requirement(:(pare[ieetab, :, 1]), 0.986180493884576)) #Combustion Eff
push!(mis_opt, Requirement(:(parg[igTmetal]), 1318.619764671319)) #Maximum metal temperature

# Constraints for this optimization
con_opt = Constraint[]
push!(con_opt, Constraint(:(wing.layout.span); pen_sca=1e4, lim_up=35.814)) # [m] Type C wing span constraint
push!(con_opt, Constraint(:(parm[imlBF, 1]); pen_sca=1e4, lim_up=2400.0)) # [m] Maximum balanced field length for takeoff
push!(con_opt, Constraint(:(para[iagamV, ipclimbn, 1]); pen_sca=1e4, lim_lo=0.015)) # [rad] Minimum cruise climb angle at TOC
push!(con_opt, Constraint(:(pare[ieTt3, :, 1]); pen_sca=1e4, lim_up=900.0)) # [K] Maximum compressor outlet temperature
push!(con_opt, Constraint(:(pare[ieTmet1, :, 1]); pen_sca=1e4, lim_up=:(parg[igTmetal]), eps_buff=1e-4)) # [K] Maximum turbine metal temperature
push!(con_opt, Constraint(:(parg[igdfan]); pen_sca=1e4, lim_up=2.0)) # [m] Maximum fan diameter
push!(con_opt, Constraint(:(parm[imWTO, 1]); pen_sca=1e4, lim_up=:(parg[igWMTO]), eps_buff=1e-4)) # [N] Maximum takeoff weight
push!(con_opt, Constraint(:(parm[imVfuel, 1]); pen_sca=1e4, lim_up=:(parg[igVfmax]), eps_buff=1e-4)) # [m3] Maximum fuel volume
# Parameters to be optimized (Initial guess(Only for local search), Upper bound, Lower bound, Initial Step Size(Only for local search))
# Engine-cycle parameters (pif/pilc/pihc/BPR/Tt4) stay free here -- only the technology level above is fixed.
bound_glob = Parameter[] #Step size = 1/5 of local search span = 1/5 * 1/4 * global search span specified below
push!(bound_glob, Parameter(:(wing.layout.sweep), 30.0, 60.0, 0.0, 5.0)) # [deg] Wing sweep angle                            #1
push!(bound_glob, Parameter(:(wing.layout.AR), 10.0, 20.0, 5.0, 1.5)) # Wing aspect ratio                                    #2
push!(bound_glob, Parameter(:(wing.inboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05)) # Wing inner t/c        #3
push!(bound_glob, Parameter(:(wing.outboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05))                       #4
push!(bound_glob, Parameter(:(wing.inboard.λ), 0.7, 1.0, 0.1, 0.1)) # Wing inner taper ratio                                 #5
push!(bound_glob, Parameter(:(wing.outboard.λ), 0.3, 1.0, 0.1, 0.1))                                                        #6
push!(bound_glob, Parameter(:(parg[igyeng]), -1e10, -1e10, -1e10, 1e10)) # [m] Span wise engine location (Placeholder)       #7
push!(bound_glob, Parameter(:(para[iarcls, ipclimb2:ipdescent4, 1]), 1.0, 2.0, 0.4, 0.15)) # Wing cl ratio break to root     #8
push!(bound_glob, Parameter(:(para[iarclt, ipclimb2:ipdescent4, 1]), 1.0, 2.0, 0.4, 0.15)) # Wing cl ratio tip to root       #9
push!(bound_glob, Parameter(:(para[iaCL, ipclimb2:ipdescent4, 1]), 0.6, 1.00, 0.3, 0.07)) # Wing CL at cruise                #10
push!(bound_glob, Parameter(:(para[iaalt, ipcruise1, 1]), 10000.0, 20000.0, 4000.0, 1600.0)) # [m] Cruise altitude           #11
push!(bound_glob, Parameter(:(pare[iepif, ipclimb2:ipdescent4, 1]), 2.0, 4.0, 1.25, 0.2)) # Fan PR at cruise                 #12
push!(bound_glob, Parameter(:(pare[iepilc, ipclimb2:ipdescent4, 1]), 3.0, 10.0, 1.25, 0.4)) # LPC PR at cruise               #13
push!(bound_glob, Parameter(:(pare[iepihc,ipclimb2:ipdescent4, 1]), 10.0, 60.0, 1.25, 1.0)) # HPC PR at cruise               #14
push!(bound_glob, Parameter(:(pare[ieBPR, ipclimb2:ipdescent4, 1]), 6.0, 12.0, 1.0, 0.5)) # Fan BPR at cruise                #15
push!(bound_glob, Parameter(:(pare[ieTt4, ipclimb2:ipdescent4, 1]), 1500.0, 2000.0, 1000.0, 100.0)) #[K] Tt4 at cruise       #16
push!(bound_glob, Parameter(:(vtail.layout.AR), 2.0, 5.0, 1.0, 0.2)) #Vertical tail aspec ratio                             #17

# test_engine_opt_LEAP1B_hard_ratio_M78_210S.jl calibration's own converged engine-cycle+geometry
# solution (bestSol_eng, TechRun#190), keyed by bound_glob index (1=Sweep,2=AR,10=CL,12=pif,
# 13=pilc,14=pihc,15=BPR,16=Tt4). pilc/pihc were split from that calibration's single per_stage
# value (1.3266058396370546) via per_stage^3/per_stage^10 -- see docstring. Used below as the
# local-search initial guess for these Parameters instead of extracting a value from the
# (pre-tech-level-swap) warm-start aircraft.
freesplit_eng_x0 = Dict(
    1  => 24.797218150225078,#Sweep
    2  => 10.008717264427581,#AR
    10 => 0.6296627190292609,#CLcruise
    12 => 1.7216518740797253,# pif
    13 => 2.33467113619516,# pilc
    14 => 16.881811513717103,# pihc
    15 => 7.693550401127311,# BPR
    16 => 1521.6812558715578,# Tt4
)

#### Setup the save directory
save_dir_actual = joinpath(save_dir,save_key*"_") #Actual directory to save data
mkpath(save_dir_actual)

#### Create log file
status_log = joinpath(save_dir_actual, "$(save_key)_$(round(Int,range_opti_nmi))_OptLog.txt")
open(status_log, "w") do io
    println(io, "Range_nmi,Status")
end

#### Load the warm-start model
ac = quickload_aircraft(par_path_base_prefix*"$(round(Int,2400)).jld2")
@assert ac.is_sized[1] #Make sure mission 1 is sized

#### Material
ac.wing.inboard.caps.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.outboard.caps.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.inboard.webs.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.wing.outboard.webs.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.fuselage.skin.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=3.0)
ac.fuselage.cone.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=2.0)
ac.fuselage.bendingmaterial_h.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)
ac.fuselage.bendingmaterial_v.material = ac.fuselage.bendingmaterial_h.material
ac.fuselage.floor.material = TASOPT.materials.StructuralAlloy("Al-2024-T4"; max_avg_stress=1.1, safety_factor=1.5)

ac.para[iaMach, ipclimbn:ipdescent1, 1] .= 0.78

#### Update the fuselage radius and layout
num_pass_row = Int(6) # Number of passengers per row
num_pax = 210 # 737-8's FAA exit limit (Boeing's own "Maximum Seats" figure), matching M78_210S
wei_per_pass_N = 956.36773
ac.parm[imWperpax,:] .= wei_per_pass_N 
ac.parg[igWpaymax] = wei_per_pass_N*num_pax 
ac.parm[imWpay,:] .= ac.parg[igWpaymax]*(1-1e-10) 
ac.fuselage.cabin.exit_limit = num_pax 
ac.fuselage.APU.W = 0.035*ac.parg[igWpaymax]
ac.fuselage.seat.W = 0.10*ac.parg[igWpaymax]
ac.fuselage.added_payload.W = 0.35*ac.parg[igWpaymax]
fuse_radius = min_fuse_radius(num_pass_row, ac) #Minimum fuselage radius given the number of seats per row
ac.fuselage.layout.cross_section.radius = fuse_radius
ac.fuselage.cabin.front_seat_offset = 0.0  # [m] Front seat offset from the front of cylinder
ac.fuselage.cabin.rear_seat_offset = 0.0 # [m] Back seat offset from the back of cylinder
upd_fuse_pax!(ac) #Update cabin geometry

#### Use the fuselage radius to update the mission and bound requirements
mis_opt[15].val = ac.wing.layout.box_x - fuse_radius
mis_opt[16].val = fuse_radius
bound_glob[7].val = fuse_radius * 3.0
bound_glob[7].bon_up = fuse_radius * 10.0
bound_glob[7].bon_lo = fuse_radius * 1.75 #unrealistic but should not have happened
bound_glob[7].d_val = fuse_radius * 0.2

#### Setup special initial optimization parameters if choose to not run global search first but by using a warm start from a default model to directly do local search
par_opt = deepcopy(bound_glob) #Use directly global bound if starting from global search
if flag_skip_global
    for (j,para_cur) in enumerate(par_opt)
        if haskey(freesplit_eng_x0, j)
            # Engine-cycle parameter -- use the Free Split candidate's own value as the initial
            # guess, not the (pre-tech-level-swap) warm-start aircraft's stale cycle state.
            para_cur.val = freesplit_eng_x0[j]
        else
            # Get starting value from the warm start aircraft model
            try # Many of the variable (Need to check why cruise value were set to many phases)
                para_cur.val = getNestedProp_fromExpr(ac; parm_sym=para_cur.name)[5] #Hopefully get the cruise phase value if possible
            catch
                para_cur.val = getNestedProp_fromExpr(ac; parm_sym=para_cur.name)[1]
            end
        end
        # Correct the local bounds be at least contained within global bounds
        bon_up = para_cur.bon_up
        bon_lo = para_cur.bon_lo
        para_cur.val = clamp(para_cur.val, bon_lo, bon_up)
        (span_glo_to_loc<1.0) || error("Local bound span needs to be smaller than global")
        dx = (bon_up - bon_lo)*span_glo_to_loc
        para_cur.bon_up = para_cur.val+0.5*dx
        para_cur.bon_lo = para_cur.val-0.5*dx
        if para_cur.bon_up>bon_up
            para_cur.bon_up = bon_up
            para_cur.bon_lo = bon_up - dx
        elseif para_cur.bon_lo<bon_lo
            para_cur.bon_lo = bon_lo
            para_cur.bon_up = bon_lo + dx
        end
    end
end

#### Save the mission requirement
save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(round(Int,range_opti_nmi))_mission_requirements.csv"), mis_opt)
save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(round(Int,range_opti_nmi))_design_constraints.csv"), con_opt)
save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(round(Int,range_opti_nmi))_global_bounds.csv"), bound_glob)

#### Run the global local optimization process
par_opt_found, status_found, hist_found =
optimizer_wrapper_global_local(ac, par_opt, bound_glob; miss_req=mis_opt, constraints=con_opt, max_iter_sizing=max_iter_sizing,
                               optimizer_global=optimizers[1], max_iter_glo=max_num_iter_opt[1],   span_glo_to_loc=span_glo_to_loc, run_global=!flag_skip_global,
                               optimizer_local =optimizers[2], max_iter_loc_C=max_num_iter_opt[2], max_round_loc_C=max_num_round_loc[1],
                                                               max_iter_loc_F=max_num_iter_opt[3], max_round_loc_F=max_num_round_loc[2])

#### Post-process the optimization results
if status_found in success_statuses
    # attempt to update an aircraft model to save
    ac_copy = deepcopy(ac)
    sized_succeeded = size_aircraft_w_param!(ac_copy; mission_req=mis_opt, parameters=par_opt_found, max_iter_sizing=max_iter_sizing)
    # correct the status
    status_found = sized_succeeded ? status_found : (:FAILED_TO_REPRODUCE)
    # Only save the new model if successfuly sized
    if sized_succeeded
        quicksave_aircraft(ac_copy, joinpath(save_dir_actual, "$(save_key)_$(round(Int,range_opti_nmi)).jld2"))
    end
end

#### save the results
open(status_log, "a") do io
    println(io, "$(range_opti_nmi),$(string(status_found))")
end
save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(round(Int,range_opti_nmi))_optimized_parameters.csv"), par_opt_found)
save_jld2(joinpath(save_dir_actual, "$(save_key)_$(round(Int,range_opti_nmi))_optimization_history.jld2"), hist_found)

println("Status: ", status_found)
println("DONE")
