"""
This script read in a baseline model and a series of off-design cases from maybe another model.
Then, the off-design cases are run on the basedline model.
New off-design models are then saved.
"""

using TASOPT
# using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(@__DIR__, "offdesign.jl"))
using .PRD: off_design_specified

#### Setup IO
model_dir    = "ModelSaved" # Setup Directory
des_model_prefix = "acOptimized_BatOptJet" # Setup baseline model for offdesign to read
offdes_model_prefix = "acOptimized_BatOptJet" # Setup off-design model to read
save_model_prefix = "acOptimized_BatOptJet_rerun" # Setup off-design result to save to 
range_des = [3000.0] #Float64.(collect(300:100:3000)) # Setup design mission
range_offdes = [3000.0]#Float64.(collect(300:100:3000)) # Setup offdesign mission
# Setup fuel to run the off-design missions
fuel_idx = 24        #Eth: 32 ,       Jet: 24
rho_fuel = 817.0     #Eth: 789.0 ,    Jet: 817.0 (kg/m3)
hvap_fuel = 358694.0 #Eth: 918187.9 , Jet: 358694.0 (J/kg)

#### Run the off-design missions
for (i, range_des_cur) in enumerate(range_des)
    #### Collect off-design missions requirements
    weight_payload_to_test = Vector{Float64}() #(N)
    range_to_test = Vector{Float64}() #(nmi)
    for (j, range_offdes_cur) in enumerate(range_offdes)
        offdes_file_name = joinpath(model_dir,offdes_model_prefix,offdes_model_prefix*"$(round(Int,range_offdes_cur)).jld2")
        println("Attempt to read: $(offdes_file_name)")
        try
            ac_offdes = quickload_aircraft(offdes_file_name)
            push!(weight_payload_to_test, ac_offdes.parm[imWpay,1]) #(N) off-design mission payload
            push!(range_to_test, ac_offdes.parm[imRange,1]/1852.0 ) #(nmi) off-design range
        catch
            ac_offdes = nothing
        end
    end
    println("For design range $(range_des_cur) nmi, read off-design ranges $(range_to_test) nmi")

    #### Collect the design mission model
    des_file_name = joinpath(model_dir,des_model_prefix ,des_model_prefix*"$(round(Int,range_des_cur)).jld2")
    println("Attempt to get the design model: $(des_file_name)")
    ac_des = nothing
    try
        ac_des = quickload_aircraft(des_file_name)
    catch
        ac_des = nothing
        continue
    end

    #### Use the design model to run off-design missions
    # Create the sub-directory for this saving
    save_dir_sub = joinpath(model_dir,save_model_prefix)
    # Run offdesign
    println("Attempt to run off-design for the design range model of R: $(round(Int,range_des_cur))")
    try
        out_off = off_design_specified(ac_des, fuel_idx, rho_fuel, hvap_fuel, range_to_test, weight_payload_to_test;
                                       save_dir = save_dir_sub, save_name = save_model_prefix*"$(round(Int,range_des_cur))_")
        println("Get off-design for ranges: $(out_off["range_nmi"])")
        println("Get off-design with PFEI: $(out_off["PFEI_JJ"])")
    catch err
        showerror(stderr, err, catch_backtrace())
    end
end