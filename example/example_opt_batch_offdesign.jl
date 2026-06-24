"""
This script is a modified version of example_opt_batch to specially design an aircraft to match the R1 R2 R3 mission requirements for an jet fuel retrofitted back to ethanoal aircraft
"""

using CSV
using DataFrames
using TASOPT
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__, "../example/utilities_for_optimization/objective_factory.jl"))
using .ObjectiveFactory: Constraint, Parameter, OptHistory, Requirement, optimize_singlePt_PFEI!, optimizer_wrapper_global_local, size_aircraft_w_param!,
                         save_vec_struct_csv, load_csv_constraints, load_csv_parameters, load_csv_requirements, save_jld2, load_jld2, OffDesMission
const success_statuses = ObjectiveFactory.success_statuses

#### Optimization parameters
# Prefix for baseline aircraft parameters as a starting point
par_path_base_prefix = joinpath(__TASOPTroot__,"../example/ModelSaved/Opti_Jet_NoACT_/Opti_Jet_NoACT_") #Also the prefixed for aircraft model
# CSV for R1 R2 R3 Missions requirements
R1_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R1.csv")
R2_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R2.csv")
R3_csv_path = joinpath(__TASOPTroot__,"../example/ModelSaved/R1R2R3_Jet_NoACT_to_Eth_/R1R2R3_Jet_NoACT_to_Eth_R3.csv")
# Path to save the models from optimization
save_dir = joinpath(__TASOPTroot__,"../example/ModelSaved/")
save_key = "Opti_Eth_NoACT_Mat_Jet_NoACT_to_Eth"
# 1.1) Optimization configuration parameters
max_iter_sizing = 150 # Maximum iterations for TASOPT sizing
optimizers = [:GN_CRS2_LM, :LN_NELDERMEAD] # Optimizer choice. [Global,Local]
max_num_iter_opt = [30000, 800, 10000] # Maximum number of optimization steps [Global, Local Coarse, Local Fine]
max_num_round_loc = [120, 30] #Maximum number of adaptive bounds refinements rounds for local search [Local Coarse, Local Fine]
span_glo_to_loc = 0.25 # span_local_search/span_global_search
# 1.2) Setup an entry (value does not matter) for test ratio to be sweep later in optimization loop
mis_opt = Requirement[]
push!(mis_opt, Requirement(:(options.has_ACT_fuel), false))  # Whether to allow additional center fuel tank (ACT)
push!(mis_opt, Requirement(:(options.compensate_ACT), false)) # Whether to increase the aircraft length to accmondate for the cargo space taken by ACT
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_vol), 1.00))  # Volumetric efficiency of ACT fuel tank
push!(mis_opt, Requirement(:(fuse_tank.ACT_eta_wei), 1.00))  # Gravimetric efficiency of ACT fuel tank
# 3) Constraints for this optimization
con_opt = Constraint[]
push!(con_opt, Constraint(:(wing.layout.span); pen_sca=1e4, lim_up=35.814)) # [m] Type C wing span constraint
push!(con_opt, Constraint(:(parm[imlBF, 1]); pen_sca=1e4, lim_up=2400.0)) # [m] Maximum balanced field length for takeoff
push!(con_opt, Constraint(:(para[iagamV, ipclimbn, 1]); pen_sca=1e4, lim_lo=0.015)) # [rad] Minimum cruise climb angle at TOC
push!(con_opt, Constraint(:(pare[ieTt3, :, 1]); pen_sca=1e4, lim_up=900.0)) # [K] Maximum compressor outlet temperature
push!(con_opt, Constraint(:(pare[ieTmet1, :, 1]); pen_sca=1e4, lim_up=1333.33)) # [K] Maximum turbine metal temperature
push!(con_opt, Constraint(:(parg[igdfan]); pen_sca=1e4, lim_up=2.0)) # [m] Maximum fan diameter
push!(con_opt, Constraint(:(parm[imWTO, 1]); pen_sca=1e4, lim_up=:(parg[igWMTO]), eps_buff=1e-4)) # [N] Maximum takeoff weight
push!(con_opt, Constraint(:(parm[imVfuel, 1]); pen_sca=1e4, lim_up=:(parg[igVfmax]), eps_buff=1e-4)) # [m3] Maximum fuel volume
# 4) Parameters to be optimized (Initial guess(Only for local search), Upper bound, Lower bound, Initial Step Size(Only for local search))
bound_glob = Parameter[] #Step size = 1/5 of local search span = 1/5 * 1/4 * global search span specified below
push!(bound_glob, Parameter(:(wing.layout.sweep), 30.0, 60.0, 0.0, 5.0)) # [deg] Wing sweep angle 
push!(bound_glob, Parameter(:(wing.layout.AR), 10.0, 20.0, 5.0, 1.5)) # Wing aspect ratio
push!(bound_glob, Parameter(:(wing.inboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05)) # Wing inner thickness to chord
push!(bound_glob, Parameter(:(wing.outboard.cross_section.thickness_to_chord), 0.2, 0.6, 0.04, 0.05))
push!(bound_glob, Parameter(:(wing.inboard.λ), 0.7, 1.0, 0.1, 0.1)) # Wing inner taper ratio
push!(bound_glob, Parameter(:(wing.outboard.λ), 0.3, 1.0, 0.1, 0.1))
push!(bound_glob, Parameter(:(parg[igyeng]), -1e10, -1e10, -1e10, 1e10)) # [m] Span wise engine location (Placeholder value here)
push!(bound_glob, Parameter(:(wing.layout.ηs), 0.3, 0.4, 0.2, 0.02)) # Wing span break location over half span
push!(bound_glob, Parameter(:(para[iarclt, ipclimb2:ipdescent4, 1]), 1.0, 2.0, 0.4, 0.15)) # Wing cl ratio of tip to root at cruise
push!(bound_glob, Parameter(:(para[iaCL, ipclimb2:ipdescent4, 1]), 0.6, 1.00, 0.3, 0.07)) # Wing CL at cruise
push!(bound_glob, Parameter(:(para[iaalt, ipcruise1, 1]), 10000.0, 20000.0, 4000.0, 1600.0)) # [m] Cruise altitude
push!(bound_glob, Parameter(:(pare[iepif, ipclimb2:ipdescent4, 1]), 2.0, 4.0, 1.25, 0.2)) # Fan PR at cruise 
push!(bound_glob, Parameter(:(pare[iepihc,ipclimb2:ipdescent4, 1]), 10.0, 50.0, 1.25, 1.0)) # HPC PR at cruise
push!(bound_glob, Parameter(:(pare[ieBPR, ipclimb2:ipdescent4, 1]), 8.0, 30.0, 1.0, 2.0)) # Fan BPR at cruise
push!(bound_glob, Parameter(:(pare[ieTt4, ipclimb2:ipdescent4, 1]), 1500.0, 2000.0, 1000.0, 100.0)) #[K] Turbine inlet temperature  at cruise
push!(bound_glob, Parameter(:(vtail.layout.AR), 2.0, 5.0, 1.0, 0.2)) #Vertical tail aspec ratio 

#### Save the directory to save the optimized solution
save_dir_actual = joinpath(save_dir,save_key*"_") #Actual directory to save data
mkpath(save_dir_actual)

#### Status log initialization
status_log = joinpath(save_dir_actual, "$(save_key)_optimization_log.txt")
open(status_log, "w") do io
    println(io, "Range_nmi,Converge_status")
end

# Initial mission requirement (Only flight range is assumed mutated later on)
save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_mission_requirements.csv"), mis_opt)
# Constraints (Fixed throughout the optimization)
save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_design_constraints.csv"), con_opt)

#### Read in the R1 R2 R3 requirements
df_R1 = CSV.read(R1_csv_path, DataFrame)
df_R2 = CSV.read(R2_csv_path, DataFrame)
df_R3 = CSV.read(R3_csv_path, DataFrame)
ran_des_nmi  = Vector{Float64}(df_R1.design_range_nmi)
ran_des_nmi  = round.(Int, ran_des_nmi) #Round the original jet aircraft design range just for name in saving
wei_pay_N_R1 = Vector{Float64}(df_R1.payload_weight_N)
ran_nmi_R1   = Vector{Float64}(df_R1.range_nmi)
wei_pay_N_R2 = Vector{Float64}(df_R2.payload_weight_N)
ran_nmi_R2   = Vector{Float64}(df_R2.range_nmi)
wei_pay_N_R3 = Vector{Float64}(df_R3.payload_weight_N)
ran_nmi_R3   = Vector{Float64}(df_R3.range_nmi)

#### Optimization
for (i,name_cur_range) in enumerate(ran_des_nmi)
    # Load the baseline aircraft model as a starting guess
    ac = quickload_aircraft(par_path_base_prefix*"$(name_cur_range).jld2")
    
    # Uset the baseline fuselage radius to update the y engine bounds as optimization parameters
    global bound_glob[7].val = ac.fuselage.layout.cross_section.radius * 3.0
    global bound_glob[7].bon_up = ac.fuselage.layout.cross_section.radius * 4.0
    global bound_glob[7].bon_lo = ac.fuselage.layout.cross_section.radius * 2.0
    global bound_glob[7].d_val = ac.fuselage.layout.cross_section.radius * 0.2
    
    # Finalize global bounds
    save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(name_cur_range)_global_bounds.csv"), bound_glob)
    
    # Load basline design parameters
    global par_opt = load_csv_parameters(par_path_base_prefix*"$(name_cur_range)_optimized_parameters.csv")
    for (idx_cur, par_cur) in enumerate(par_opt) #Update the initial step size
        par_cur.d_val = bound_glob[idx_cur].d_val
    end

    # Design mission follows the baseline setup. Specify off-design missions requirements here
    range_off_des_nmi = [ran_nmi_R1[i], ran_nmi_R2[i], ran_nmi_R3[i]]
    wei_pay_off_des_N = [wei_pay_N_R1[i], wei_pay_N_R2[i], wei_pay_N_R3[i]]
    idx_fuel_off_des = fill(Int(32), 3)
    rho_fuel_off_des_kgm3 = fill(789.0, 3)
    hvap_fuel_off_des_Jkg = fill(918187.9, 3)
    off_des_miss = OffDesMission{Float64}(range_off_des_nmi, wei_pay_off_des_N, idx_fuel_off_des, rho_fuel_off_des_kgm3, hvap_fuel_off_des_Jkg)
    off_des_constraints = [:WPay,:MWTO,:VolFuel]
    PFEI_Weighting = [0.0,0.0,1.0,0.0] # Here we only use R2 PFEI for optimization

    # Run the global local optimization process
    par_opt_found, status_found, hist_found = 
    optimizer_wrapper_global_local(ac, par_opt, bound_glob; miss_req=mis_opt, constraints=con_opt, max_iter_sizing=max_iter_sizing,
                                   off_des_miss=off_des_miss, off_des_constraints=off_des_constraints, PFEI_Weighting=PFEI_Weighting,
                                   optimizer_global=optimizers[1], max_iter_glo=max_num_iter_opt[1],   span_glo_to_loc=span_glo_to_loc, run_global=false,
                                   optimizer_local =optimizers[2], max_iter_loc_C=max_num_iter_opt[2], max_round_loc_C=max_num_round_loc[1],
                                                                   max_iter_loc_F=max_num_iter_opt[3], max_round_loc_F=max_num_round_loc[2])
    
    # Post-process the optimization results
    if status_found in success_statuses
        # attempt to update an aircraft model to save
        ac_copy = deepcopy(ac)
        sized_succeeded = size_aircraft_w_param!(ac_copy; mission_req=mis_opt, parameters=par_opt_found, max_iter_sizing=max_iter_sizing)
        # correct the status
        status_found = sized_succeeded ? status_found : (:FAILED_TO_REPRODUCE)
        # Only save the new model if successfuly sized
        if sized_succeeded
            quicksave_aircraft(ac_copy, joinpath(save_dir_actual, "$(save_key)_$(name_cur_range).jld2"))
        end
    end
    
    # save the results
    open(status_log, "a") do io
        println(io, "$(name_cur_range),$(string(status_found))")
    end
    save_vec_struct_csv(joinpath(save_dir_actual, "$(save_key)_$(name_cur_range)_optimized_parameters.csv"), par_opt_found)
    save_jld2(joinpath(save_dir_actual, "$(save_key)_$(name_cur_range)_optimization_history.jld2"), hist_found)
end