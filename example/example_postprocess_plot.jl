"""
This script compare the design point performance between ethanol and jet
Adaptive design parameter bounds across different cases
Fixed and default constraints across different cases
"""

using TASOPT
using DataFrames, CSV
include(__TASOPTindices__)
include(joinpath(@__DIR__, "postprocess.jl"))
using .PostProcess: ExtractDes
using .PostProcess.OptimizeRangeFuel: BoundsOpt, ConstraintsOpt, extract_opt_para
using Plots

#### File Name Creation
function modelFileName(model_prefix::String, curFuel::String, curRange::String)
    #return "$(model_prefix)$(curFuel)$(curRange)"
    return "$(curFuel)$(model_prefix)$(curRange)"
end

#### Setup IO
model_dir    = "ModelSaved"
model_prefix = "_300_3000_SpanLar_FanLar" #Frontal key name for the models(FuelRange)
save_dir     = "ModelProcessed" #Outer Directory for saving models
save_prefix  = "Eth_Jet_Comparison" #The total save name for jet and ethanol comparison
mkpath(model_dir)
mkpath(save_dir)
save_dirSub = joinpath(save_dir,save_prefix) #Sub-directory to save comparison data
mkpath(save_dirSub)

#### Missions Parameters to Compare
Fuels = ["JetFuel","Ethanol"] #These corresponding to the model file name
Ranges = collect(300:100:3000) #Prefixex+Fuel+Range.jld2
flgSaveIndividual = false #Whether to save an output for individual case
flgPlotBounds = false #Whether to plots the optimization bounds
flgPlotConstraints = true #Whether to plots the constraints

#### Load the default constraints which is believed at this point to be common to call cases
constraints_opt = ConstraintsOpt()
constraints_opt.DiaFan_max = 3.0 #Overwrite
constraints_opt.span_max = 65.0 #Overwrite
constraints_names = ["diaFan", "TMetalMax", "Tt3Max", "gamTOC", "lenFieldBalanced", "spanWing"]
constraints = [constraints_opt.DiaFan_max,
               constraints_opt.TMetal_max,
               constraints_opt.Tt3_max,
               constraints_opt.TOCGamma_min*180.0/pi, #deg for consistency
               constraints_opt.lenField_max,
               constraints_opt.span_max]

#### Extract data for each case
paraSet = [] # Each element is one fuel
boundSet = [] # Each element is one fuel
optParSet = []
for (i,curFuel) in enumerate(Fuels)
    paraSetSub = [] # Each element is one range
    boundSetSub = [] # Each element is one range
    optParSetSub = []
    for (j,curRange) in enumerate(Ranges)
        caseCur = "$(curFuel)$(round(Int,curRange))"
        println("Attempt to read case: Fuel type: $(curFuel), Range: $(curRange)")
        try
            # Load the aircraft model
            ac = quickload_aircraft(joinpath(model_dir,"$(modelFileName(model_prefix, curFuel, string(round(Int,curRange)))).jld2"))
            # Load the parameter bounds
            bdcsv = CSV.read(joinpath(model_dir,"$(modelFileName(model_prefix, curFuel, string(round(Int,curRange))))_BoundLocal.csv"), DataFrame)
            bounds_opt = BoundsOpt()
            for fName in fieldnames(typeof(bounds_opt))
                fNameSym = Symbol(fName)
                setfield!(bounds_opt, fName, (Float64(bdcsv[1, fNameSym]), Float64(bdcsv[2, fNameSym]), Float64(bdcsv[3, fNameSym])))
            end
            # Extract the design parameters and optionally save individual case
            design_para = ExtractDes(ac, save_dir, "$(model_prefix)$(caseCur)"; flg_save=flgSaveIndividual,
                          bounds_opt=bounds_opt, constraints_opt=constraints_opt)
            # Extract the subset of design parameters that have bounds corresponds to
            optimi_para = extract_opt_para(ac) #vector{Float64} Those parameters having bounds. Same order as the bound
            # Append the data
            push!(paraSetSub, (; range_id = string(round(Int,curRange)), fuel_id = curFuel, design_para...))
            push!(boundSetSub, (bounds_opt))
            push!(optParSetSub, (optimi_para))
        catch err
            println("Reading case failed for Fuel type: $(curFuel), Range: $(curRange)")
        end
    end
    push!(paraSet, (paraSetSub))
    push!(boundSet, (boundSetSub))
    push!(optParSet, (optParSetSub))
end

####Plot out the operation parameter comparison between the fuels over the ranges
fNamesOptPara = fieldnames(BoundsOpt)
numOptParas = length(fNamesOptPara)
# Get the range out
range_matrix = []
for idxFuel in eachindex(Fuels) #Through each fuel
    # Read out the ranges first
    range_lst = []
    for idxRange in eachindex(paraSet[idxFuel])
        push!(range_lst, paraSet[idxFuel][idxRange].range) #[nmi]
    end
    push!(range_matrix, range_lst)
end
# Get and compare each design parameters
for idxPara = 1:numOptParas #Loop through all the parameters
    p = plot(xlabel="Range [nmi]", ylabel="$(fNamesOptPara[idxPara])", dpi=800)
    for idxFuel in eachindex(Fuels) #Through each fuel    
        # Get parameters
        value_lst = []
        for idxRange in eachindex(optParSet[idxFuel]) # Assume optparaset have the same number of valid fuel and ranges as the paraset
            push!(value_lst, optParSet[idxFuel][idxRange][idxPara])
        end
        # Plot parameters
        plot!(p, range_matrix[idxFuel], value_lst, marker=:cross, lw=2, label=Fuels[idxFuel])
        # Get bounds
        if flgPlotBounds
            LB_lst = []
            UB_lst = []
            for idxRange in eachindex(boundSet[idxFuel])
                push!(LB_lst, getfield(boundSet[idxFuel][idxRange], fNamesOptPara[idxPara])[1])
                push!(UB_lst, getfield(boundSet[idxFuel][idxRange], fNamesOptPara[idxPara])[2])
            end
            # Plot bounds
            plot!(p, range_matrix[idxFuel], LB_lst, lw=0.5, label="$(Fuels[idxFuel]) Low Bound")
            plot!(p, range_matrix[idxFuel], UB_lst, lw=0.5, label="$(Fuels[idxFuel]) High Bound")
        end
    end
    savefig(p, joinpath(save_dirSub, "CompBound_$(fNamesOptPara[idxPara]).png"))
end

#### Plot out the constraints parameter comparison between the fuels
for idxPara = 1:length(constraints_names) #Loop through all the parameters
    p = plot(xlabel="Range [nmi]", ylabel="$(constraints_names[idxPara])", dpi=800)
    for idxFuel in eachindex(Fuels) #Through each fuel    
        # Get constrained parameters
        value_lst = []
        for idxRange in eachindex(paraSet[idxFuel]) # Assume optparaset have the same number of valid fuel and ranges as the paraset
            push!(value_lst, getfield(paraSet[idxFuel][idxRange], Symbol(constraints_names[idxPara])))
        end
        # Plot parameters
        plot!(p, range_matrix[idxFuel], value_lst, marker=:cross, lw=2, label=Fuels[idxFuel])
        # Plot constraints
        if flgPlotConstraints
            plot!(p, range_matrix[idxFuel], fill(constraints[idxPara], length(range_matrix[idxFuel])), lw=0.5, label="$(Fuels[idxFuel]) Constraints")
        end
    end
    savefig(p, joinpath(save_dirSub, "CompConstr_$(constraints_names[idxPara]).png"))
end

#### Plot out other design parameters
other_para_names = ["PFEI", "LD_cruise","LHV_cruise", "TSFC_cruise", "TSEC_cruise", "vel_cruise", "massTO", "massFuelReserved", "massFuelBurned", "massEmpty", "rangeBreguet", "thrustOneEngine_takeoff", "thrustOneEngine_cruise1"]
for idxPara = 1:length(other_para_names) #Loop through all the parameter names
    p = plot(xlabel="Range [nmi]", ylabel="$(other_para_names[idxPara])", dpi=800)
    for idxFuel in eachindex(Fuels) #Through each fuel    
        # Get constrained parameters
        value_lst = []
        for idxRange in eachindex(paraSet[idxFuel]) # Assume optparaset have the same number of valid fuel and ranges as the paraset
            push!(value_lst, getfield(paraSet[idxFuel][idxRange], Symbol(other_para_names[idxPara])))
        end
        # Plot parameters
        plot!(p, range_matrix[idxFuel], value_lst, marker=:cross, lw=2, label=Fuels[idxFuel])
    end
    savefig(p, joinpath(save_dirSub, "CompOtherPara_$(other_para_names[idxPara]).png"))
end