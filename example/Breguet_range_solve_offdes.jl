module Breguet

    export Bre_off_des

    """
    Bre_off_des(ranges_m, wPay_N;
            LD, eta, LHV_Jkg, wEmp_N, rhoFuel_kgm3, frac_rese,
            wTO_Max_N, wPay_Max_N, volFuel_Max_m3,
            gee = 9.81)
    Output:
        Dict["range_m_out", "wPay_N_out", "wFuel_N_out", "mFuel_kg_out", "PFEI_JJ_out"]
    """
    function Bre_off_des(ranges_m::Real, wPay_N::Real;
                        LD::Real, eta::Real, LHV_Jkg::Real, wEmp_N::Real, rhoFuel_kgm3::Real, frac_rese::Real,
                        wTO_Max_N::Real, wPay_Max_N::Real, volFuel_Max_m3::Real,
                        gee::Real = 9.81)

        #### Promote to Float64 once
        ranges_m = Float64(ranges_m)
        wPay_N = Float64(wPay_N)
        LD = Float64(LD)
        eta = Float64(eta)
        LHV_Jkg = Float64(LHV_Jkg)
        wEmp_N = Float64(wEmp_N)
        rhoFuel_kgm3 = Float64(rhoFuel_kgm3)
        frac_rese = Float64(frac_rese)
        wTO_Max_N = Float64(wTO_Max_N)
        wPay_Max_N = Float64(wPay_Max_N)
        volFuel_Max_m3 = Float64(volFuel_Max_m3)
        gee = Float64(gee)

        # Check for inputs
        @assert (ranges_m>=0.0)
        @assert (wPay_N>=0.0)
        @assert (LD>0)
        @assert (eta>0)
        @assert (LHV_Jkg>0)
        @assert (wEmp_N>0)
        @assert (rhoFuel_kgm3>0)
        @assert (frac_rese>=0.0 && frac_rese < (exp((ranges_m*gee)/(LHV_Jkg*eta*LD))-1)^(-1))
        @assert (wTO_Max_N>=(wEmp_N+wPay_Max_N) && wTO_Max_N>=(wEmp_N+volFuel_Max_m3*rhoFuel_kgm3*gee))
        @assert (wPay_Max_N>0.0)
        @assert (volFuel_Max_m3>0.0)
        @assert (gee>0)
        if wPay_N == 0.0
            wPay_N = 1e-10 #(N)
        end
        if ranges_m == 0.0
            ranges_m = 1e-10 #(m)
        end
        
        #### Constants
        C1 = LHV_Jkg * eta * LD / gee
        C2 = exp(ranges_m/C1) - 1.0
        C3 = C2/(1.0 - C2*frac_rese) #wFuel = C3*(wEmp + wPay)

        #### Caps Calculation
        wPay_N_UB1 = wPay_Max_N #(N)
        wPay_N_UB2 = wTO_Max_N / (1.0 + C3 * (1.0 + frac_rese)) - wEmp_N #(N)
        wPay_N_UB3 = (volFuel_Max_m3 * gee * rhoFuel_kgm3) / (C3 * (1.0 + frac_rese)) - wEmp_N #(N)
        wPay_N_max = min(wPay_N_UB1,wPay_N_UB2,wPay_N_UB3)
        if wPay_N_max<0.0
            return Dict("range_m_out" => -1.0,
                        "wPay_N_out" => -1.0,
                        "wFuel_N_out" => -1.0,
                        "mFuel_kg_out" => -1.0,
                        "PFEI_JJ_out" => -1.0)
        end
        
        #### Test feasibility
        range_m_out  = ranges_m
        wPay_N_out   = 0.0
        if (wPay_N>wPay_N_UB3) && (wPay_N>wPay_N_UB2) && (wPay_N>wPay_N_UB1)
            wPay_N_out = wPay_N_max
        elseif (wPay_N>wPay_N_UB3) && (wPay_N>wPay_N_UB2)
            wPay_N_out = min(wPay_N_UB2,wPay_N_UB3)
        elseif (wPay_N>wPay_N_UB1) && (wPay_N>wPay_N_UB3)
            wPay_N_out = min(wPay_N_UB1,wPay_N_UB3)
        elseif (wPay_N>wPay_N_UB2) && (wPay_N>wPay_N_UB1)
            wPay_N_out = min(wPay_N_UB1,wPay_N_UB2)
        elseif (wPay_N>wPay_N_UB1)
            wPay_N_out = wPay_N_UB1
        elseif (wPay_N>wPay_N_UB2)
            wPay_N_out = wPay_N_UB2
        elseif (wPay_N>wPay_N_UB3)
            wPay_N_out = wPay_N_UB3
        else
            wPay_N_out = wPay_N
        end
        wFuel_N_out = C3*(wEmp_N + wPay_N_out)
        mFuel_kg_out = wFuel_N_out/gee
        PFEI_JJ_out = (wPay_N_out*range_m_out)/(mFuel_kg_out*LHV_Jkg) 
        
        #### Output
        return Dict("range_m_out" => range_m_out,
                    "wPay_N_out" => wPay_N_out,
                    "wFuel_N_out" => wFuel_N_out,
                    "mFuel_kg_out" => mFuel_kg_out,
                    "PFEI_JJ_out" => PFEI_JJ_out)

    end


end #end module