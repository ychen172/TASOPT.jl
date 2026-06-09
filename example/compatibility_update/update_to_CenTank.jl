"""
This script update the model saved from current branch to the CenTank branch readable format
"""

using TASOPT
using JLD2
include(__TASOPTindices__)

@kwdef mutable struct fuselage_tank_new
    """Fuel type name"""
    fueltype::String = ""
    """Fuel tank count"""
    tank_count::Int64 = 0
    """Fuel tank location"""
    placement::String = ""
    """Flag for insulation sizing"""
    sizes_insulation::Bool = false
    """Weight of fuel in one tank (N)"""
    Wfuelintank::Float64 = 0.0

    clearance_fuse::Float64 = 0.0

    """Vector with insulation layer thickness (m)"""
    t_insul::Vector{Float64} = Float64[]
    """Vector with insulation materials"""
    material_insul::Vector{ThermalInsulator} = Float64[]
    """Vector with insulation layer design indices"""
    iinsuldes::Vector{Int64} = Float64[]
    """Length of cylindrical portion of tank (m)"""
    l_cyl_inner::Float64 = 0.0
    """Length of inner tank (m)"""
    l_inner::Float64 = 0.0
    """Inner tank radius (m)"""
    Rinnertank::Float64 = 0.0
    """Vector with surface areas of insulation tank heads (m^2)"""
    Shead_insul::Vector{Float64} = Float64[]

    """Inner vessel material"""
    inner_material::StructuralAlloy = StructuralAlloy("Al-2219-T87")
    """Outer vessel material"""
    outer_material::StructuralAlloy = StructuralAlloy("Al-2219-T87")
    """Tank head aspect ratio"""
    ARtank::Float64 = 0.0
    """Angular location of inner vessel stiffeners"""
    theta_inner::Float64 = 0.0
    """Vector with angular location of outer vessel stiffeners"""
    theta_outer::Vector{Float64} = Float64[]
    """Number of intermediate stiffeners in outer vessel"""
    Ninterm::Float64 = 1.0
    
    """Venting pressure (Pa)"""
    pvent::Float64 = 0.0
    """Fill pressure (Pa)"""
    pinitial::Float64 = 0.0
    """Minimum allowable tank pressure (Pa)"""
    pmin::Float64 = 0.0
    """Departure hold time (s)"""
    t_hold_orig::Float64 = 0.0
    """Arrival hold time (s)"""
    t_hold_dest::Float64 = 0.0
    """Sea-level temperature for tank design (K)"""
    TSLtank::Vector{Float64} = []

    """Liquid fuel density (kg/m^3)"""
    rhofuel::Float64 = 0.0
    """Liquid fuel temperature in tank (K)"""
    Tfuel::Float64 = 0.0
    """Gas fuel density (kg/m^3)"""
    rhofuelgas::Float64 = 0.0
    """Fuel specific enthalpy of vaporization (J/kg)"""
    hvap::Float64 = 0.0
    """Percentage tank boiloff rate at start of cruise (%/h)"""
    boiloff_rate::Float64 = 0.0

    """Vessel additional mass fraction"""
    ftankadd::Float64 = 0.0
    """Vessel weld efficiency"""
    ew::Float64 = 0.0
    """Minimum ullage fraction"""
    ullage_frac::Float64 = 0.0
    """Heat leakage factor"""
    qfac::Float64 = 0.0
    """Pressure rise factor"""
    pfac::Float64 = 0.0

    """For addtional fuel tank"""
    ACT_eta_vol::Float64 = 0.88
    ACT_eta_wei::Float64 = 0.93
    ACT_A::Float64 = 0.0
    ACT_l::Float64 = 0.0
    ACT_W::Float64 = 0.0
    ACT_dx::Float64 = 0.0
    ACT_fuse_l_extend::Float64 = 0.0
end

function fuselage_tank_new(fuse_tank::fuselage_tank)
    return fuselage_tank_new(;
        fueltype = fuse_tank.fueltype,
        tank_count = fuse_tank.tank_count,
        placement = fuse_tank.placement,
        sizes_insulation = fuse_tank.sizes_insulation,
        Wfuelintank = fuse_tank.Wfuelintank,
        clearance_fuse = fuse_tank.clearance_fuse,
        t_insul = copy(fuse_tank.t_insul),
        material_insul = copy(fuse_tank.material_insul),
        iinsuldes = copy(fuse_tank.iinsuldes),
        l_cyl_inner = fuse_tank.l_cyl_inner,
        l_inner = fuse_tank.l_inner,
        Rinnertank = fuse_tank.Rinnertank,
        Shead_insul = copy(fuse_tank.Shead_insul),
        inner_material = fuse_tank.inner_material,
        outer_material = fuse_tank.outer_material,
        ARtank = fuse_tank.ARtank,
        theta_inner = fuse_tank.theta_inner,
        theta_outer = copy(fuse_tank.theta_outer),
        Ninterm = fuse_tank.Ninterm,
        pvent = fuse_tank.pvent,
        pinitial = fuse_tank.pinitial,
        pmin = fuse_tank.pmin,
        t_hold_orig = fuse_tank.t_hold_orig,
        t_hold_dest = fuse_tank.t_hold_dest,
        TSLtank = copy(fuse_tank.TSLtank),
        rhofuel = fuse_tank.rhofuel,
        Tfuel = fuse_tank.Tfuel,
        rhofuelgas = fuse_tank.rhofuelgas,
        hvap = fuse_tank.hvap,
        boiloff_rate = fuse_tank.boiloff_rate,
        ftankadd = fuse_tank.ftankadd,
        ew = fuse_tank.ew,
        ullage_frac = fuse_tank.ullage_frac,
        qfac = fuse_tank.qfac,
        pfac = fuse_tank.pfac
    )
end

@kwdef mutable struct aircraft_tmp{WS}
    name::String
    description::String
    options::TASOPT.Options

    parg::Vector{Float64}
    parm::Array{Float64,2}
    para::Array{Float64,3}
    pare::Array{Float64,3}

    is_sized::Vector{Bool}

    fuselage::TASOPT.Fuselage
    fuse_tank::Union{TASOPT.fuselage_tank, fuselage_tank_new}

    wing::TASOPT.Wing
    htail::TASOPT.Tail
    vtail::TASOPT.Tail
    engine::TASOPT.Engine
    landing_gear::TASOPT.LandingGear

    wake_system::WS
end

function aircraft_tmp(ac::aircraft)
    return aircraft_tmp(
        ac.name,
        ac.description,
        ac.options,
        ac.parg,
        ac.parm,
        ac.para,
        ac.pare,
        ac.is_sized,
        ac.fuselage,
        fuselage_tank_new(ac.fuse_tank),  # <-- conversion happens here
        ac.wing,
        ac.htail,
        ac.vtail,
        ac.engine,
        ac.landing_gear,
        ac.wake_system
    )
end

function update_model!(ac::aircraft_tmp)
    parg = ac.parg #Vector{Float64}
    parm = ac.parm #Array{Float64, 2}
    pare = ac.pare #Array{Float64, 3}

    push!(parg,parg[igWfmax]/parg[igrhofuel]) #igVfmax
    push!(parg,parg[igWfuel]/parg[igrhofuel]) #igVfuel

    parm = vcat(parm, fill(parg[igWfuel]/parg[igrhofuel],size(parm,2))') #imVfuel
    pare = cat(pare, fill(parg[igrhofuel], 1, size(pare,2), size(pare,3)); dims=1)#ierhofuel_driven

    ac.parg = parg
    ac.parm = parm
    ac.pare = pare
end

#
loadpath = joinpath(@__DIR__,"../ModelSaved/acOptimized_BatOptEth/acOptimized_BatOptEth300.jld2")
savepath = joinpath(@__DIR__,"../ModelSaved/testSave.jld2")
@load loadpath ac
ac = aircraft_tmp(ac)
update_model!(ac)
@save savepath ac
