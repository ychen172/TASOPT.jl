"""
This script compares parameters across multiple cases
"""

using TASOPT
include(__TASOPTindices__)
using Plots
include(joinpath(@__DIR__,"../utilities_for_postprocessing/Extract.jl"))
using .Extract: extract_acModel_compact!, init_results_2Layers, plot_cases_specified
include(joinpath(@__DIR__,"../Breguet_range_solve_offdes.jl"))
using .Breguet: Bre_off_des

#### Setup IO
# Input case names - Retrofit
model_dir  = "../ModelSaved"
caseKeys   = ["Opti_Jet_NoACT_V4_4_R1Sz_EtasEng_","Opti_Eth_NoACT_V4_3_R1Sz_EtasEng_"]
caseNames  = ["Jet Fuel"                         ,"Ethanol"                          ]
ranges     = collect(300:100:3000) 
# Output directory
save_dir      = "../ModelProcessed"
save_name     = "Aug_3_2026_Jet_Fuel_vs_Ethanol_Conventional" #sub_folder will be created
# Fields to read out
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

#### Initialization
dataset = [init_results_2Layers(length(ranges), fields) for _ in caseKeys] #[dataset[Expr][:],...]

#### Extract data for the sized missions
for (j, caseKey) in enumerate(caseKeys)
    for (i, ran) in enumerate(ranges)
        # Read in the case
        ac_dir = joinpath(model_dir,caseKey,caseKey*"$(round(Int,ran)).jld2")
        ac = quickload_aircraft(ac_dir)
        println("File, $(ac_dir), read successfully")
        # Extract data
        extract_acModel_compact!(ac, dataset[j], i)
        println("Data case $(j) at range $(ran) collected successfully")
    end
end

#### Plotting
plot_cases_specified("Design/R1 Range [nmi]", "Flight PFEI [J/J]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(parm[imPFEI,1])] for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"PFEI.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Fuel Level [%]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(parm[imVfuel,1])] .* 100.0 ./ d[:(parg[igVfmax])] for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"Fuel_Level.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Wingspan [m]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(wing.layout.span)] for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"Wingspan.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Taper Ratio of Outer Wing Panel", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(wing.outboard.λ)] for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"OuterTaperWing.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Wing Span Break Location [%]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(wing.layout.ηs)] .* 100.0 for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"SpanBreakWing.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Wing Aspect Ratio", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(wing.layout.AR)]  for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"AspectRatioWing.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Lift-to-drag Ratio at Start of Cruise", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(para[iaCL,ipcruise1,1])] ./ d[:(para[iaCD,ipcruise1,1])]  for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"LD_StartOfCruise.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Lift-to-drag Ratio at End of Cruise", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(para[iaCL,ipcruise2,1])] ./ d[:(para[iaCD,ipcruise2,1])]  for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"LD_EndOfCruise.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Engine Total Efficiency at Start of Cruise [%]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [100.0 .* (1.0 ./ (d[:(pare[ieTSFC,ipcruise1,1])] ./ gee)) .* ((cos.(d[:(para[iagamV,ipcruise1,1])]) .* d[:(pare[ieu0,ipcruise1,1])]) ./ d[:(pare[iehfuel,ipcruise1,1])])  for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"TotEffEngine_StartOfCruise.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Engine Total Efficiency at End of Cruise [%]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [100.0 * (1.0 ./ (d[:(pare[ieTSFC,ipcruise2,1])] ./ gee)) .* ((cos.(d[:(para[iagamV,ipcruise2,1])]) .* d[:(pare[ieu0,ipcruise2,1])]) ./ d[:(pare[iehfuel,ipcruise2,1])])  for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"TotEffEngine_EndOfCruise.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Aircraft Empty Weight [Ton]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [(d[:(parm[imWTO,1])] .- d[:(parm[imWfuel,1])] .- d[:(parm[imWpay,1])]) ./ gee ./ 1000.0  for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"EmptyWeight.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Fuel Weight Fraction of Total Takeoff Weight [%]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [100.0 .* (d[:(parm[imWfuel,1])] ./ d[:(parm[imWTO,1])]) for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"FuelWeightFraction.png"))
#
plot_cases_specified("Design/R1 Range [nmi]", "Wing Surface Area [m²]", 
                     [d[:(parm[imRange,1])] for d in dataset] ./ 1852.0, 
                     [d[:(wing.layout.S)] for d in dataset],
                     caseNames, 
                     joinpath(save_dir_sub,"WingSurfaceArea.png"))
#

#### Special comparative plots
#Calculate using Breguet range an estimated penalty for ethanol versus jet fuel use
PFEIEthPenalty_Breguet = zeros(length(ranges)) #[%]
iRC = 1 #Index Reference Case
iRM = 28 #Index Reference Mission
iRJ = 1 #Index of the jet case
iRE = 2 #Index of the ethanol case
for idxMiss in eachindex(ranges) #Go through each mission pair between ethanol and jet fuel case
    outJet = Bre_off_des(dataset[iRC][:(parm[imRange,1])][idxMiss], dataset[iRC][:(parm[imWpay,1])][idxMiss];
                         LD = dataset[iRC][:(para[iaCL,ipcruise1,1])][iRM]/dataset[iRC][:(para[iaCD,ipcruise1,1])][iRM],
                         eta = (1.0/(dataset[iRC][:(pare[ieTSFC,ipcruise1,1])][iRM]/gee))*((cos(dataset[iRC][:(para[iagamV,ipcruise1,1])][iRM])*dataset[iRC][:(pare[ieu0,ipcruise1,1])][iRM])/dataset[iRC][:(pare[iehfuel,ipcruise1,1])][iRM]),
                         LHV_Jkg = dataset[iRJ][:(pare[iehfuel,ipcruise1,1])][iRM], 
                         wEmp_N = dataset[iRC][:(parm[imWTO,1])][iRM]-dataset[iRC][:(parm[imWfuel,1])][iRM]-dataset[iRC][:(parm[imWpay,1])][iRM],
                         rhoFuel_kgm3 = dataset[iRJ][:(parg[igrhofuel])][iRM], frac_rese=dataset[iRC][:(parg[igfreserve])][iRM],
                         wTO_Max_N=4e10, wPay_Max_N=1e10, volFuel_Max_m3=1e10/1e5) #Not considering physical limitations
    outEth = Bre_off_des(dataset[iRC][:(parm[imRange,1])][idxMiss], dataset[iRC][:(parm[imWpay,1])][idxMiss];
                         LD = dataset[iRC][:(para[iaCL,ipcruise1,1])][iRM]/dataset[iRC][:(para[iaCD,ipcruise1,1])][iRM],
                         eta = (1.0/(dataset[iRC][:(pare[ieTSFC,ipcruise1,1])][iRM]/gee))*((cos(dataset[iRC][:(para[iagamV,ipcruise1,1])][iRM])*dataset[iRC][:(pare[ieu0,ipcruise1,1])][iRM])/dataset[iRC][:(pare[iehfuel,ipcruise1,1])][iRM]),
                         LHV_Jkg = dataset[iRE][:(pare[iehfuel,ipcruise1,1])][iRM], 
                         wEmp_N = dataset[iRC][:(parm[imWTO,1])][iRM]-dataset[iRC][:(parm[imWfuel,1])][iRM]-dataset[iRC][:(parm[imWpay,1])][iRM],
                         rhoFuel_kgm3 = dataset[iRE][:(parg[igrhofuel])][iRM], frac_rese=dataset[iRC][:(parg[igfreserve])][iRM],
                         wTO_Max_N=4e10, wPay_Max_N=1e10, volFuel_Max_m3=1e10/1e5) #Not considering physical limitations
    PFEIEthPenalty_Breguet[idxMiss] = 100.0*(outEth["PFEI_JJ_out"]-outJet["PFEI_JJ_out"])/outJet["PFEI_JJ_out"]
end

plot_cases_specified("Design/R1 Range [nmi]", "Ethanol Flight PFEI Penalty [%]", 
                     [dataset[iRC][:(parm[imRange,1])],dataset[iRC][:(parm[imRange,1])]] ./ 1852.0, 
                     [100.0 .* (dataset[iRE][:(parm[imPFEI,1])] .- dataset[iRJ][:(parm[imPFEI,1])]) ./ dataset[iRJ][:(parm[imPFEI,1])],PFEIEthPenalty_Breguet],
                     ["Optimized Parameters","Constant Efficiencies (Breguet)"], 
                     joinpath(save_dir_sub,"EthPFEIPenalty.png"))
#
