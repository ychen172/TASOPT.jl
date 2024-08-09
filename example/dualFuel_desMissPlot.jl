# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOPT. 

# 1) Load TASOPT
using TASOPT
using PyPlot
using CSV
# using Tables
include("../src/misc/index.inc")
ft_to_m = 0.3048
# From Optimization
"""
               AR                  Alt(ft)             Cl                  Λ                   λs                  λt                   hboxo                hboxs                rcls                rclt                Tt4CR               iepihc              iepif
Ethanol 1500: [11.941686657713525, 30511.776452022423, 0.5848967822299598, 27.712995630599664, 0.9417718748745083, 0.10412119292629382, 0.14811779467488953, 0.12580371789052663, 1.0545034510486877, 0.9989229403529972, 1749.8606928474742, 14.977689688881302, 1.4301673363334233] PFEI = 0.6513
Jet A   1500: [11.991949818973007, 31284.521653455868, 0.5930853805805898, 27.45888882421664, 0.958879620436785, 0.11555475645318984, 0.14367703488865205, 0.12197546585940724, 1.0331998484328888, 0.9999708411275662, 1768.1487526911849, 14.989386166042735, 1.879087280526366]    PFEI = 0.6167
EthaJet 1500: [11.209746212440592, 31248.335606979665, 0.5822564781297511, 26.83688322046656, 0.8645278884243015, 0.1381719896560668, 0.14225631183926352, 0.12078733440024264, 1.1565051727391427, 0.9998672032093319, 1756.579167574885, 14.909003005180068, 1.700769613859498]     PFEI = 0.6308 0.6410
JetEtha 1500: [11.319998631871222, 30811.23725202253, 0.5813759825551024, 27.640718319157607, 0.8772736830151879, 0.1350984828369109, 0.14855570828902398, 0.12530518480239944, 1.1631814900480104, 0.9419705949013939, 1753.3494643968827, 15.0, 1.3887491288940033]                 PFEI = 0.6485 0.6394
Ethanol 3000: [11.972782021599212, 26855.31676655534, 0.5897546909709177, 27.934121648805956, 0.9083625579263621, 0.10038536693707056, 0.148118228816289, 0.133653028913655, 1.1045845807165509, 0.9991970828125056, 1763.2463808321709, 14.991966665820446, 1.8797278394132353]      PFEI = 0.7957
Jet A   3000: [12.0, 29036.283196079105, 0.5847390857049931, 27.536926874931076, 0.9455429779743656, 0.11162958807659483, 0.14804764447122798, 0.12667338569408848, 1.052751755460775, 0.9837147694563335, 1789.0037092463112, 14.995243338018897, 1.3567628966201042]                PFEI = 0.6756
EthaJet 3000: [11.767524327190367, 28084.037784162752, 0.5818414565848145, 27.700212481634704, 0.9300673536492501, 0.107470417949662, 0.14761018071108312, 0.12914782388715706, 1.0792445068013392, 0.9980608775454582, 1797.21150483555, 14.989353182560986, 1.8997543951239186]     PFEI = 0.7045 0.7109
JetEtha 3000: [11.73825703104477, 26821.182812272666, 0.5800928093129452, 27.071979085742164, 0.8202211816967958, 0.10011141740947224, 0.14286434753223193, 0.12515531292922516, 1.2320048914143524, 0.9137484062312613, 1780.3630579099809, 14.998545567272956, 1.2207582475268322]  PFEI = 0.7679 0.7613
Ethanol 2250: [11.995790823232323, 27436.801350118487, 0.5826911328437688, 27.634381158682476, 0.9224809175036062, 0.1039364578365561, 0.1490086413142229, 0.12791090503589803, 1.0853705461074912, 0.9985612550494998, 1783.9146612129932, 14.997065344994738, 1.8989347314942426]   PFEI = 0.7076
Jet A   2250: [11.973393618436805, 31254.582327577267, 0.5860052516023335, 27.68704614674552, 0.9519924570118375, 0.10395349365694001, 0.1481992196333408, 0.1254146503963879, 1.0539986029999036, 1.0, 1763.694866386413, 14.997854217768928, 1.8588863082387639]                    PFEI = 0.6384
EthaJet 2250: [11.995940864779865, 30464.73778738585, 0.591496130187772, 27.607515749377402, 0.9471598891712532, 0.1020463692216292, 0.1462669140493254, 0.12545260382405868, 1.054347591062121, 0.9991422231224164, 1769.124077885052, 14.996258137948766, 1.8584298224056814]       PFEI = 0.6565 0.6640
JetEtha 2250: [11.997237766721248, 28290.769979357814, 0.582912373048561, 27.66710351806454, 0.9325984194554195, 0.10536834277053876, 0.14877874745298114, 0.12738253318079706, 1.0764419072918554, 0.9149268452567285, 1779.8579748353382, 14.904886391112868, 1.795734616563026]    PFEI = 0.6915 0.6848
EJ 1500Dual : [11.99451901200959, 31533.65389444799, 0.5846243009397081, 27.74958240195238, 0.9407468249931482, 0.10589384713663758, 0.14780787433897713, 0.1250055804350304, 1.0720774392606942, 0.8731648387365094, 1745.2488393511476, 14.571589150524822, 1.8887525996672818]     PFEI = 0.6366
JE 1500Dual : [12.0, 31013.979762758965, 0.5844366591095647, 27.68596643876029, 0.9396100888724693, 0.10324499406287652, 0.1482481847912658, 0.12588178783610118, 1.058410576889806, 0.9936746721704857, 1746.4961555956363, 14.898876895521129, 1.8999999999999995]                  PFEI = 0.6342
"""
savedir = "Movie/"
if !isdir(savedir)
    # If it doesn't exist, create the "optimization" directory
    mkdir(savedir)
end
# 2) Include input file for desired aircraft/
nameAircraftModel = "../src/IO/experiment_input_2250.toml"
ac = read_aircraft_model(nameAircraftModel) # MODIFY <path> appropriately
saveName = savedir*"JetAEtha2250_230"
x = [11.997237766721248, 28290.769979357814, 0.582912373048561, 27.66710351806454, 0.9325984194554195, 0.10536834277053876, 0.14877874745298114, 0.12738253318079706, 1.0764419072918554, 0.9149268452567285, 1779.8579748353382, 14.904886391112868, 1.795734616563026]
# 2.5) Change fuel type
ac.pari[iifuel] = 322471 #(JetA:25 Ethanol:32 JetAEtha29%Blend: 322429 JetAEtha71%Blend: 322471)
ac.parg[igrhofuel] = 794.504 #(JetA:817.0 Ethanol:789.0 JetAEtha29%Blend: 805.649 JetAEtha71%Blend: 794.504)
# 3) Set the parameters based on optimization result
ac.parg[igAR] = x[1] # Aspect Ratio 
ac.para[iaalt, ipcruise1, :] .=  x[2] * ft_to_m # Cruise Altitude
ac.para[iaCL, ipclimb1+1:ipdescentn-1, :] .= x[3] # CL
ac.parg[igsweep] = x[4] # Wing sweep 
ac.parg[iglambdas] = x[5] #inner_panel_taper_ratio
ac.parg[iglambdat] = x[6] #outer_panel_taper_ratio
ac.parg[ighboxo] = x[7] #root_thickness_to_chord
ac.parg[ighboxs] = x[8] #spanbreak_thickness_to_chord
ac.para[iarcls, ipclimb1+1 : ipdescentn-1, :] .= x[9]   #  rcls    break/root cl ratio = cls/clo
ac.para[iarclt, ipclimb1+1 : ipdescentn-1, :] .= x[10]   #  rclt    tip  /root cl ratio = clt/clo
ac.pare[ieTt4, ipcruise1:ipcruise2, :] .= x[11] # Tt4
ac.pare[iepihc, ipclimb1+1 : ipdescentn-1, :] .= x[12] # High Pressure Compressor Pressure Ratio
ac.pare[iepif, ipclimbn, :] .= x[13] #Fan PR 
ac.pare[iepilc, :, :] .= 3 # Low Pressure Compressure Pressure Ratio set to 3
# 3) Size aircraft
TASOPT.size_aircraft!(ac, iter =500, printiter=true)
# 4) Collect data
AltRec        = [x[2]] #ft
RanRec        = [ac.parg[igRange]./1852.0] #nmi
WMTORec       = [ac.parg[igWMTO]./(9.81*1000)] #Ton (metric)
WFuelRec      = [ac.parg[igWfuel]./(9.81*1000)] #Ton (metric)
WPayRec       = [ac.parg[igWpay]./(9.81*1000)] #Ton (metric)
PFEIRec       = [ac.parm[imPFEI,1]] #(J/J)
WTO_WTOmaxRec = [ac.parm[imWTO,1]/ac.parg[igWMTO]] #N/N
Wf_WfmaxRec   = [ac.parg[igWfuel]/ac.parg[igWfmax]] #N/N
areaWingRec   = [ac.parg[igS]] #m2
ARWingRec     = [ac.parg[igAR]] #Aspect ratio
spanWingRec   = [ac.parg[igb]] #m
diaFanRec     = [ac.parg[igdfan]] #m
FnTotCRRec    = [ac.pare[ieFe,ipcruise1,1]] #N
WEmpRec       = [WMTORec - WFuelRec - WPayRec] # Ton Metric
EFuelRec      = [PFEIRec.*WPayRec*1000*9.81.*RanRec*1852.0] #Joul
SweepRec      = [ac.parg[igsweep]] #deg
## Add additional output on cargo bay fuel volume
Wfmax = ac.parg[igWfmax]
Wf    = ac.parg[igWfuel]
# Compute the additional fuel tank volume available from the cargo bay
AFuse = ac.parg[igAfuse] #fuselage crosssection area [m2]
lShell = ac.parg[igxshell2]-ac.parg[igxshell1] #length of the cylindrical sector [m]
WCargo = AFuse*0.45*lShell*ac.parg[igrhofuel]*9.81 #[N]
# Finish additional fuel tank calculation
if ((Wf-Wfmax)>0)
    fCargoFuel = [(Wf - Wfmax)/WCargo] #Percentage volume taken in cargo space of fuel storage
else
    fCargoFuel = [0.0]
end

outputTup = (AltRec=AltRec,RanRec=RanRec,WMTORec=WMTORec,WFuelRec=WFuelRec
             ,WPayRec=WPayRec,PFEIRec=PFEIRec,WTO_WTOmaxRec=WTO_WTOmaxRec
             ,Wf_WfmaxRec=Wf_WfmaxRec,areaWingRec=areaWingRec,ARWingRec=ARWingRec
             ,spanWingRec=spanWingRec,diaFanRec=diaFanRec,FnTotCRRec=FnTotCRRec,WEmpRec=WEmpRec,EFuelRec=EFuelRec,SweepRec=SweepRec,fCargoFuel=fCargoFuel)
CSV.write(saveName*"MissDetail.csv",  outputTup, writeheader=true)
##Create a mask to mask out unreported phases
maskRep = ac.pare[ieFe,:,1].>0 #Reported Phase has non zero thrust
phases = ["ST","RO","TO","CB","B1","B2","B3","B4","B5","C1","C2","D1","D2","D3","D4","D5","Test"]
phases = phases[maskRep]
print("Reported Phases Are:", phases,"\n")
##Read out other parameters
###Aero Parameters
timeOptMiss = ac.para[iatime,maskRep,1] #second
ranOptMiss = ac.para[iaRange,maskRep,1] #meter
altOptMiss = ac.para[iaalt,maskRep,1] #meter
machOptMiss = ac.para[iaMach,maskRep,1]
weiOptMiss = ac.para[iafracW,maskRep,1]*(ac.parm[imWTO,1]/9.81) #kg
gamOptMiss = ac.para[iagamV,maskRep,1] #rad
LDROptMiss = ac.para[iaCL,maskRep,1]./ac.para[iaCD,maskRep,1] #Lift to drag ratio
###Engine Parameters
hfOptMiss = ac.pare[iehfuel,maskRep,1] #J/kg equivalent heating value
TfuelOptMiss = ac.pare[ieTfuel,maskRep,1] #K fuel temperature
Tt3OptMiss = ac.pare[ieTt3,maskRep,1] #K
Pt3OptMiss = ac.pare[iept3,maskRep,1] #Pa 
Tt4OptMiss = ac.pare[ieTt4,maskRep,1] #K
Pt4OptMiss = ac.pare[iept4,maskRep,1] #Pa
FnOptMiss = ac.pare[ieFe,maskRep,1] #N Total Thrust for each engines
mdotfOptMiss = ac.pare[iemcore,maskRep,1].*ac.pare[ieff,maskRep,1] #kg/s for eacg engines
Cpa = 0.5.*(ac.pare[iecpt3,maskRep,1].+ac.pare[iecpt4,maskRep,1])
ffbMiss   = (Cpa.*(Tt4OptMiss.-Tt3OptMiss))./(hfOptMiss.*ac.pare[ieetab,maskRep,1].-Cpa.*(Tt4OptMiss.-TfuelOptMiss))
mdot3OptMiss = mdotfOptMiss./ffbMiss #kg/s for each engines air flow into the combustor (exclude bypass cooling flow)
#Output Additional Data at the optimal mission
outputTup = (Phase=phases,Time=timeOptMiss,Range=ranOptMiss,Altitude=altOptMiss,MachNumber=machOptMiss,Weight=weiOptMiss
            ,ClimbAngle=gamOptMiss,LiftDragRatio=LDROptMiss,HeatingValue=hfOptMiss,FuelTemp=TfuelOptMiss
            ,Tt3=Tt3OptMiss,Pt3=Pt3OptMiss,Tt4=Tt4OptMiss,Pt4=Pt4OptMiss,Thrust=FnOptMiss,mdotFuel=mdotfOptMiss,mdot3=mdot3OptMiss)
CSV.write(saveName*"MissDetail2.csv",  outputTup, writeheader=true)
#Plot Plane
TASOPT.stickfig(ac, label_fs = 8)
plt.savefig(saveName*"MissDetail3.png")