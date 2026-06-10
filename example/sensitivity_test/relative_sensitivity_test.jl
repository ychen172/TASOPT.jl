"""
This script for a given optimal mission find the relative sensitivity of PFEI
to individiual design parameters, and plot them in descending order.
As the sensitivity to design parameters might be weak near the optimal point but
still want to evaluate the sensentivity around optimum, the design mission range is 
set to be shorter than the actual sizing range for the optimized mission.    
"""

using TASOPT
include(TASOPT.__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils","sensitivity.jl"))
using Plots
using LaTeXStrings

function sensitivityTest(modelPath,caseName,eps,saveDir,sizedRange)
    ####Load the default aircraft and setup some sizign option
    println("Reading: $(modelPath)")
    mkpath(saveDir)
    ac = quickload_aircraft(modelPath)
    ac.parm[imRange,:] .= sizedRange*1852.0 #(m)
    ac.htail.opt_sizing = TailSizing.CLmaxFwdCG #Use lift base horizontal tail sizing
    ac.htail.CL_max_fwd_CG = -0.7
    ac.vtail.opt_sizing = TailSizing.OEI #Use engine out vertical tail sizing
    ac.parg[igCLveout] = 0.5
    #### Input parameters to be size(need to be precise)
    input_params = [
        :(ac.wing.layout.ηs),                                     # Panel break eta location
        :(ac.parg[igxeng]),                                       # Engine axial location
        :(ac.parg[igyeng]),                                       # Engine vertical location
        :(ac.wing.layout.AR),                                     # Wing aspect ratio
        :(ac.para[iaCL, ipcruise1, 1]),                           # Cruise CL
        :(ac.wing.layout.sweep),                                  # Wing sweep
        :(ac.para[iaalt, ipcruise1, 1]),                          # Cruise altitude
        :(ac.wing.inboard.λ),                                     # Inboard taper ratio
        :(ac.wing.outboard.λ),                                    # Outboard taper ratio
        :(ac.wing.inboard.cross_section.thickness_to_chord),      # Inboard t/c
        :(ac.wing.outboard.cross_section.thickness_to_chord),     # Outboard t/c
        :(ac.para[iarcls, ipcruise1, 1]),                         # Break/root Cl ratio
        :(ac.para[iarclt, ipcruise1, 1]),                         # Tip/root Cl ratio
        :(ac.pare[ieTt4, ipcruise1, 1]),                          # Tt4 at cruise
        :(ac.pare[iepihc, ipcruise1, 1]),                         # HPC pressure ratio
        :(ac.pare[iepif, ipcruise1, 1]),                          # Fan pressure ratio
        :(ac.pare[iepilc, ipcruise1, 1]),                         # LPC pressure ratio
        :(ac.pare[ieBPR, ipcruise1, 1]),                          # Bypass ratio
        :(ac.htail.CL_max_fwd_CG),                                # CL max horizontal tail (probably skip)
        :(ac.parg[igCLveout]),                                    # CL max vertical tail (Skip due to duplicated effect with yEng and also bounded by physical limit)
        :(ac.parg[igcdefan]),                                     # Dead engine CD
        :(ac.para[iaMach, ipcruise1, 1]),                         # Cruise mach number
        :(ac.parm[imgamVDE1, 1]),                                 # descent_angle_top-of-descent
        :(ac.fuselage.layout.cross_section.radius),               # fuselage cross-sectional radius
        :(ac.fuselage.layout.floor_depth),                        # fuselage floor depth
        :(ac.htail.layout.AR),                                    # Horizontal tail AR
        :(ac.htail.outboard.λ),                                   # Horizontal tail taper
        :(ac.htail.layout.sweep),                                 # Horizontal tail sweep
        :(ac.vtail.layout.AR),                                    # Vertical tail AR
        :(ac.vtail.outboard.λ),                                   # Vertical tail taper
        :(ac.vtail.layout.sweep)                                  # Vertical tail sweep
    ]
    #### Input parameter names for print out only
    names_params = [
        L"\eta_{wing\:span\:break}",
        L"x_{engine}",
        L"y_{engine}",
        L"AR_{wing}",
        L"CL_{cruise}",
        L"Sweep_{wing}",
        L"Alt_{cruise}",
        L"Taper_{inner\:wing}",
        L"Taper_{outer\:wing}",
        L"Thick/C_{inner\:wing}",
        L"Thick/C_{outer\:wing}",
        L"cl_{break}/cl_{root}",
        L"cl_{tip}/cl_{root}",
        L"Tt4_{cruise}",
        L"PR_{HPC}",
        L"PR_{fan}",
        L"PR_{LPC}",
        L"BPR_{fan}",
        L"CL_{htail\:max\:fwd\:CG}",
        L"CL_{vtail\:engine\:out}",
        L"CD_{dead\:engine}",
        L"Mach_{cruise}",
        L"\gamma_{top\:of\:descent}",
        L"R_{fuselage}",
        L"Dep_{floor\:fuselage}",
        L"AR_{htail}",
        L"Taper_{htail}",
        L"Sweep_{htail}",
        L"AR_{vtail}",
        L"Taper_{vtail}",
        L"Sweep_{vtail}"
    ]

    muted_blue   = RGB(0.35, 0.55, 0.75)
    muted_green  = RGB(0.40, 0.65, 0.45)
    muted_red    = RGB(0.80, 0.45, 0.45)
    muted_orange = RGB(0.90, 0.65, 0.35)
    muted_purple = RGB(0.60, 0.50, 0.75)
    muted_gray   = RGB(0.50, 0.50, 0.50)
    muted_yellow = RGB(0.90, 0.80, 0.35)

    colors = [
        muted_gray,  # Panel break eta location :red :blue :green :orange :purple
        muted_gray,  # Engine axial location
        muted_green,  # Engine spanwise location
        muted_green,  # Wing aspect ratio
        muted_green,  # Cruise CL
        muted_green,  # Wing sweep
        muted_purple,  # Cruise altitude
        muted_gray,  # Inboard taper ratio
        muted_gray,  # Outboard taper ratio
        muted_purple,  # Inboard t/c
        muted_purple,  # Outboard t/c
        muted_gray,  # Break/root Cl ratio
        muted_gray,  # Tip/root Cl ratio
        muted_green,  # Tt4 at cruise
        muted_gray,  # HPC pressure ratio
        muted_gray,  # Fan pressure ratio
        muted_gray,  # LPC pressure ratio
        muted_gray,  # Bypass ratio
        muted_gray,  # CL max horizontal tail (probably skip)
        muted_yellow,  # CL max vertical tail (Skip due to duplicated effect with yEng and also bounded by physical limit)
        muted_gray,  # Dead engine CD
        muted_yellow,  # Cruise mach number
        muted_gray,  # descent_angle_top-of-descent
        muted_yellow,  # fuselage cross-sectional radius
        muted_gray,  # fuselage floor depth
        muted_gray,  # Horizontal tail AR
        muted_gray,  # Horizontal tail taper
        muted_gray,  # Horizontal tail sweep
        muted_gray,  # Vertical tail AR
        muted_gray,  # Vertical tail taper
        muted_gray   # Vertical tail sweep
    ]

    #### Compute sensitivities
    impactVector = get_sensitivity(input_params; model_state=ac, eps=1e-5, optimizer=false, f_out_fn=nothing, diff_scheme=:central, metric=:impact)

    #### Plot out sensitivity
    #Inpact plot
    sens = Float64.(impactVector) #From any[]
    # Sort from largest to smallest
    idx = sortperm(sens, rev=true)
    sens_sorted = sens[idx]
    names_sorted = names_params[idx]
    colors_sorted = colors[idx]
    # Bar plot
    bar(
        1:length(sens_sorted), sens_sorted*100;
        xticks = (1:length(names_sorted), names_sorted),
        color = colors_sorted,
        tick_direction = :none,
        xrotation = 90,              # vertical labels
        xtickfontsize = 11,
        yguidefontsize = 8,
        ylabel = "ΔPFEI/PFEI₀ (%) \n from $(eps*100)% ΔParameters", # xlabel = "Parameters",
        legend = false, #title = "Parameter Sensitivity (sorted, high → low)",
        dpi = 600,
        size = (1200, 400),
        yscale = :identity, #:log10 or :identity  # , #
        margin = 20Plots.mm #(Adjust this margin if label got clipped)
    )
    # Optional save
    savefig(joinpath(saveDir,caseName*"_Impact.png"))
end

####IO parameter
modelPath = joinpath(@__DIR__,"../ModelSaved/acOptim_BatJet_CT/acOptim_BatJet_CT")
modelRangeOri = [500,1000,1500,2000,2500,2900] #Range by original model at optimum
caseName = "Jet" #for figure saving
sensitivityRange = Int.(modelRangeOri .- 200) #Range intentially off-optimum of sensivitiy 
eps = 1e-5 #Use for sensitivity perturbation
saveDir = joinpath(@__DIR__,"../ModelProcessed/Sensitivity")
for (idx,range_cur) in enumerate(modelRangeOri)
    println("Run $(idx)")
    sensitivityTest(modelPath*"$(round(Int,range_cur)).jld2", caseName*"_$(round(Int,range_cur))_Run_$(round(Int,sensitivityRange[idx]))", eps, saveDir, sensitivityRange[idx])
end