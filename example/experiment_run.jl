# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOPT. 

# 1) Load TASOPT
using TASOPT
include("../src/misc/index.inc")
ft_to_m = 0.3048
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# 2) Include input file for desired aircraft/
#  load default model
# example_ac = load_default_model() # simply a synonym to read_aircraft_model()
# Alternatively you can load your desired input file 
ac = read_aircraft_model("../src/IO/experiment_input.toml") # MODIFY <path> appropriately

# 2.5) Change fuel type
# ac.pari[iifuel] = 1

# 3) Size aircraft
ac.para[iaalt, ipcruise1, 1] =  3.5e4 * ft_to_m # Cruise Altitude
time_wsize = @elapsed size_aircraft!(ac,iter=135)
#println("Time to size aircraft = $time_wsize s")

# 3.5) Read out the size of each variable
# display(size(ac.pari))
# display(size(ac.parg))
# display(size(ac.parm))
# display(size(ac.para))
# display(size(ac.pare))

# 3.75) Read out the total weight and flight range
println("flight range (nmi): " , ac.parg[igRange]./1852.0)
println("WMTO (1000 kg):" , ac.parg[igWMTO]./(9.8*1000))
println("PFEI:",ac.parm[imPFEI])
println("OPR:",ac.pare[iept3]/ac.pare[iept2])
println("LPCPR:",maximum(ac.pare[iepilc, :, 1]))
println("WTO/WTOmax:",ac.parm[imWTO,1]/ac.parg[igWMTO])
println("Wf/Wfmax:",ac.parg[igWfuel]/ac.parg[igWfmax])

# 4) Visualize outputs
# Output resulting geometry of aircraft
# summary(example_ac)
# Or individually look at certain aspects:
# Show weight buildup of the aircraft:
# TASOPT.weight_buildup(example_ac) 
# # Show aerodynamics:
# TASOPT.aero(example_ac)
# # Geometry:
# TASOPT.geometry(example_ac)

# 5) Plot figures
# using PyPlot
# TASOPT.stickfig(example_ac)
# plt.savefig("Example.png")