"""
This script compares the design parameters with the global bounds to spot any out of bounds design
"""

using TASOPT
include(__TASOPTindices__)
using Plots
using CSV
using DataFrames
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel, init_results_2Layers, fill_results!, plot_cases

#### Setup IO
# Input case names - Retrofit
model_dir  = joinpath(__TASOPTroot__,"../example/ModelSaved")
caseKeys   = "Opti_Eth_NoACT_V3_" #Parameters Comparison will be saved with an extension to this
ranges     = collect(300:100:3000)
var_names = [ #no space
    "Sweep",
    "Wing_AR",
    "Wing_Inboard_Thick_to_Chord",
    "Wing_Outboard_Thick_to_Chord",
    "Wing_Inboard_Taper_Ratio",
    "Wing_Outboard_Taper_Ratio",
    "Engine_Spanwise_Position",
    "Wing_Span_Break_Fraction",
    "Wing_Tip_Cl_Ratio",
    "Cruise_CL",
    "Cruise_Altitude",
    "Fan_PR",
    "High_Pressure_Compressor_PR",
    "Bypass_Ratio",
    "Turbine_Inlet_Temperature",
    "Vertical_Tail_AR",
    "Engine_Axial_Position"
]

# Output directory
save_dir   = joinpath(__TASOPTroot__,"../example/ModelProcessed",caseKeys*"OptParam")
mkpath(save_dir)
val_des_collected = Array{Float64}(undef, length(var_names), length(ranges))
up_bound_collected = Array{Float64}(undef, length(var_names), length(ranges))
low_bound_collected = Array{Float64}(undef, length(var_names), length(ranges))
for (i, ran) in enumerate(ranges)
    # Read in the case
    des_dir = joinpath(model_dir,caseKeys,caseKeys*"$(round(Int(ran)))_optimized_parameters.csv") #design parameters directory
    glo_dir = joinpath(model_dir,caseKeys,caseKeys*"$(round(Int(ran)))_global_bounds.csv") #global parameters directory
    
    df_des = CSV.read(des_dir, DataFrame)
    df_glo = CSV.read(glo_dir, DataFrame)

    val_des = Vector{Float64}(df_des.val_value)
    up_bound = Vector{Float64}(df_glo.bon_up_value)
    low_bound = Vector{Float64}(df_glo.bon_lo_value)
    @assert(length(var_names)==length(val_des))
    @assert(length(var_names)==length(low_bound))
    for j in eachindex(val_des)
        val_des_collected[j,i] = val_des[j]
        up_bound_collected[j,i] = up_bound[j]
        low_bound_collected[j,i] = low_bound[j]
    end
end


#### Ploting
for i in eachindex(var_names)
    p = plot(xlabel="Design Ranges (nmi)",ylabel=var_names[i],dpi=800)
    plot!(p, ranges, val_des_collected[i,:], marker=:cross, lw=2, label="optimum")
    plot!(p, ranges, up_bound_collected[i,:], marker=:xcross, lw=2, label="upper bound")
    plot!(p, ranges, low_bound_collected[i,:], marker=:xcross, lw=2, label="lower bound")
    savefig(p, joinpath(save_dir, "$(var_names[i]).png"))
end