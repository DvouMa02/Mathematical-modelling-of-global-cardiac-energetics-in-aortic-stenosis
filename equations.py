import numpy as np
from numba import njit


def initialize_activation_functions(ctype, periode, m1, m2, Emax, Emin):
    t = np.linspace(0, periode, 10000)

    if ctype == "atrium":
        tau1 = 0.11 * periode
        tau2 = 0.18 * periode

        g1 = (t/tau1) ** m1
        g2 = (t/tau2) ** m2
        k = (Emax - Emin) / max((g1/(1+g1)) * (1/(1+g2)))

        return tau1, tau2, k

    if ctype == "ventricle":
        tau1 = 0.269 * periode
        tau2 = 0.452 * periode

        g1 = (t/tau1) ** m1
        g2 = (t/tau2) ** m2
        k = (Emax - Emin) / max((g1/(1+g1)) * (1/(1+g2)))

        return tau1, tau2, k
    


def activation_function(chamber_type, x, tau1, tau2, m1, m2, k, Emin, periode=0.0, delta=0.0):
    if chamber_type == "ventricle":
        g1 = (x / tau1) ** m1
        g2 = (x / tau2) ** m2

        Et = k * (g1 / (1 + g1)) * (1 / (1 + g2)) + Emin
        return Et

    elif chamber_type == "atrium":
        v = np.copy(x)
        v[np.where(v < 0.15*periode)] = v[np.where(v < 0.15*periode)] + periode
        v[np.where(v >= 0.15*periode)] = v[np.where(v >= 0.15*periode)] - delta
        v[np.where(v < 0)] = 0

        g1 = (v / tau1) ** m1
        g2 = (v / tau2) ** m2

        Et = k * (g1 / (1 + g1)) * (1 / (1 + g2)) + Emin
        return Et
    

def Et_Rmicro(x, m1, m2, Rmin, Rmean, periode, rt_all=False):
    if isinstance(x, float):
        x = np.linspace(0, periode, 10000)
    tau1 = 0.269 * periode
    tau2 = 0.452 * periode

    g1 = (x/tau1) ** m1
    g2 = (x/tau2) ** m2
    k = (Rmean - Rmin) / np.mean(((g1/(1+g1)) * (1/(1+g2)))**2)

    Et_R = k * (((g1 / (1 + g1)) * (1 / (1 + g2)))**2) + Rmin
    if rt_all:
        return Et_R, tau1, tau2, k
    return Et_R





@njit
def system_odes_scipy_faster(t, u, params):
    
    # unpack state variables
    Vra, Vrv, Vla, Vlv, Qtv, Qpv, Qmv, Qav, eps_tv, eps_pv, eps_mv, eps_av, Pa_syst, Pv_syst, Pa_pulm, Pv_pulm, Pa_LAD, Pv_LAD, Pa_LCX, Pv_LCX, Pa_RCA, Pv_RCA = u
    
    conversion_factor, delta, periode, tau1_ra, tau2_ra, m1_ra, m2_ra, k_ra, Emin_ra, tau1_rv, tau2_rv, m1_rv, m2_rv, k_rv, Emin_rv, tau1_la, tau2_la, m1_la, m2_la, k_la, Emin_la, tau1_lv, tau2_lv, m1_lv, m2_lv, k_lv, Emin_lv, V0ra, V0rv, V0la, V0lv, density, Rao, Ra_syst, Rv_syst, Rpa, Ra_pulm, Rv_pulm, Ca_syst, Cv_syst, Ca_pulm, Cv_pulm, Amax_tv, Amin_tv, Amax_pv, Amin_pv, Amax_mv, Amin_mv, Amax_av, Amin_av, Kvo_tv, Kvc_tv, Kvo_pv, Kvc_pv, Kvo_mv, Kvc_mv, Kvo_av, Kvc_av, Ra_LAD, Ra_LCX, Ra_RCA, Rv_LAD, Rv_LCX, Rv_RCA, Ca_LAD, Ca_LCX, Ca_RCA, Cmicro_LAD, Cmicro_LCX, Cmicro_RCA, tau1_LAD, tau2_LAD, m1_LAD, m2_LAD, k_LAD, Rmicro_min_LAD, tau1_LCX, tau2_LCX, m1_LCX, m2_LCX, k_LCX, Rmicro_min_LCX, tau1_RCA, tau2_RCA, m1_RCA, m2_RCA, k_RCA, Rmicro_min_RCA = params



    # chambers
    t_mod = t % periode
    if t_mod < 0.15 * periode:
        t_atrium = t_mod + periode -delta
    else:
        t_atrium = t_mod - delta

    if t_atrium < 0.0:
        t_atrium = 0.0

    g1_ra = (t_atrium / tau1_ra) ** m1_ra
    g2_ra = (t_atrium / tau2_ra) ** m2_ra
    Era = k_ra * (g1_ra / (1 + g1_ra)) * (1 / (1 + g2_ra)) + Emin_ra
    Era = Era * conversion_factor

    g1_rv = (t / tau1_rv) ** m1_rv
    g2_rv = (t / tau2_rv) ** m2_rv
    Erv = k_rv * (g1_rv / (1 + g1_rv)) * (1 / (1 + g2_rv)) + Emin_rv
    Erv_dot = ((Erv-Emin_rv)/max(t,1e-6))*(m1_rv/(1+g1_rv) - (m2_rv*g2_rv)/(1+g2_rv))
    Erv = Erv*conversion_factor
    Erv_dot = Erv_dot*conversion_factor

    g1_la = (t_atrium / tau1_la) ** m1_la
    g2_la = (t_atrium / tau2_la) ** m2_la
    Ela = k_la * (g1_la / (1 + g1_la)) * (1 / (1 + g2_la)) + Emin_la
    Ela = Ela*conversion_factor

    g1_lv = (t / tau1_lv) ** m1_lv
    g2_lv = (t / tau2_lv) ** m2_lv
    Elv = k_lv * (g1_lv / (1 + g1_lv)) * (1 / (1 + g2_lv)) + Emin_lv
    Elv_dot = ((Elv-Emin_lv)/max(t,1e-6))*(m1_lv/(1+g1_lv) - (m2_lv*g2_lv)/(1+g2_lv))
    Elv = Elv*conversion_factor
    Elv_dot = Elv_dot*conversion_factor
    

    Pra = Era*(Vra - V0ra)
    Prv = Erv*(Vrv - V0rv)
    Pla = Ela*(Vla - V0la)
    Plv = Elv*(Vlv - V0lv)

    Qv_pulm = (Pv_pulm - Pla)/Rv_pulm #if Pv_pulm > Pla else 0.0 ##
    Qv_syst = (Pv_syst - Pra)/Rv_syst #if Pv_syst > Pra else 0.0 ##
    Qv_LAD = (Pv_LAD - Pra)/Rv_LAD #if Pv_LAD > Pra else 0.0
    Qv_LCX = (Pv_LCX - Pra)/Rv_LCX #if Pv_LCX > Pra else 0.0
    Qv_RCA = (Pv_RCA - Pra)/Rv_RCA #if Pv_RCA > Pra else 0.0
    Qv_syst_total = Qv_syst + Qv_LAD + Qv_LCX + Qv_RCA


    Vra_dot = Qv_syst_total - Qtv
    Vrv_dot = Qtv - Qpv
    Vla_dot = Qv_pulm - Qmv
    Vlv_dot = Qmv - Qav


    # valves
    Aeff_tv = (Amax_tv - Amin_tv)*eps_tv + Amin_tv
    Aeff_tv = max(Aeff_tv, 1e-9)
    Btv = density/(2*Aeff_tv**2)
    Mtv = density/(np.sqrt(np.pi*Aeff_tv))

    Aeff_pv = (Amax_pv - Amin_pv)*eps_pv + Amin_pv
    Aeff_pv = max(Aeff_pv, 1e-9)
    Bpv = density/(2*Aeff_pv**2)
    Mpv = density/(np.sqrt(np.pi*Aeff_pv))

    Aeff_mv = (Amax_mv - Amin_mv)*eps_mv + Amin_mv
    Aeff_mv = max(Aeff_mv, 1e-9)
    Bmv = density/(2*Aeff_mv**2)
    Mmv = density/(np.sqrt(np.pi*Aeff_mv))

    Aeff_ao = (Amax_av - Amin_av)*eps_av + Amin_av
    Aeff_ao = max(Aeff_ao, 1e-9)
    Bav = density/(2*Aeff_ao**2)
    Mav = density/(np.sqrt(np.pi*Aeff_ao))

    Qtv_dot = (Pra - Prv - Btv*Qtv*abs(Qtv))/Mtv
    Qpv_dot = (Prv - Pa_pulm - Bpv*Qpv*abs(Qpv))/Mpv
    Qmv_dot = (Pla - Plv - Bmv*Qmv*abs(Qmv))/Mmv
    Pao = (Qav + Pa_syst/Rao + Pa_LAD/Ra_LAD + Pa_LCX/Ra_LCX + Pa_RCA/Ra_RCA)/(1/Rao + 1/Ra_LAD + 1/Ra_LCX + 1/Ra_RCA)
    Qav_dot = (Plv - Pao - Bav*Qav*abs(Qav))/Mav

    eps_tv_dot = (1.0-eps_tv)*Kvo_tv*(Pra - Prv) if Pra - Prv >= 0 else eps_tv*Kvc_tv*(Pra - Prv)
    eps_pv_dot = (1.0-eps_pv)*Kvo_pv*(Prv - Pa_pulm) if Prv - Pa_pulm >= 0 else eps_pv*Kvc_pv*(Prv - Pa_pulm)
    eps_mv_dot = (1.0-eps_mv)*Kvo_mv*(Pla - Plv) if Pla - Plv >= 0 else eps_mv*Kvc_mv*(Pla - Plv)
    eps_av_dot = (1.0-eps_av)*Kvo_av*(Plv - Pa_syst) if Plv - Pa_syst >= 0 else eps_av*Kvc_av*(Plv - Pa_syst)
    

    # pulmonary circulation
    Pv_pulm_dot = ((Pa_pulm-Qpv*Rpa-Pv_pulm)/Ra_pulm - Qv_pulm)/Cv_pulm
    Pa_pulm_dot = (Qpv + Ca_pulm*Rpa*Qpv_dot - (Pa_pulm-Qpv*Rpa-Pv_pulm)/Ra_pulm)/Ca_pulm

    # systemic circulation
    Pv_syst_dot = ((Pa_syst-Pv_syst)/Ra_syst - Qv_syst)/Cv_syst
    Pa_syst_dot = ((Pao - Pa_syst)/Rao - (Pa_syst-Pv_syst)/Ra_syst)/Ca_syst

    # coronary circulation
    g1_LAD = (t / tau1_LAD) ** m1_LAD
    g2_LAD = (t / tau2_LAD) ** m2_LAD
    Rmicro_LAD = k_LAD * (((g1_LAD / (1 + g1_LAD)) * (1 / (1 + g2_LAD)))**2) + Rmicro_min_LAD
    #Rmicro_LAD = Rmicro_LAD * conversion_factor
    IMPlv_dot = 0.6*(Elv_dot*(Vlv-V0lv) + Elv*Vlv_dot)
    Pa_LAD_dot = ((Pao-Pa_LAD)/Ra_LAD - (Pa_LAD-Pv_LAD)/Rmicro_LAD)/Ca_LAD
    Pv_LAD_dot = ((Pa_LAD-Pv_LAD)/Rmicro_LAD + Cmicro_LAD*IMPlv_dot - Qv_LAD)/Cmicro_LAD

    g1_LCX = (t / tau1_LCX) ** m1_LCX
    g2_LCX = (t / tau2_LCX) ** m2_LCX
    Rmicro_LCX = k_LCX * (((g1_LCX / (1 + g1_LCX)) * (1 / (1 + g2_LCX)))**2) + Rmicro_min_LCX
    #Rmicro_LCX = Rmicro_LCX * conversion_factor
    Pa_LCX_dot = ((Pao-Pa_LCX)/Ra_LCX - (Pa_LCX-Pv_LCX)/Rmicro_LCX)/Ca_LCX
    Pv_LCX_dot = ((Pa_LCX-Pv_LCX)/Rmicro_LCX + Cmicro_LCX*IMPlv_dot - Qv_LCX)/Cmicro_LCX

    g1_RCA = (t / tau1_RCA) ** m1_RCA
    g2_RCA = (t / tau2_RCA) ** m2_RCA
    Rmicro_RCA = k_RCA * (((g1_RCA / (1 + g1_RCA)) * (1 / (1 + g2_RCA)))**2) + Rmicro_min_RCA
    #Rmicro_RCA = Rmicro_RCA * conversion_factor
    IMPrv_dot = 0.6*(Erv_dot*(Vrv-V0rv) + Erv*Vrv_dot)
    Pa_RCA_dot = ((Pao-Pa_RCA)/Ra_RCA - (Pa_RCA-Pv_RCA)/Rmicro_RCA)/Ca_RCA
    Pv_RCA_dot = ((Pa_RCA-Pv_RCA)/Rmicro_RCA + Cmicro_RCA*IMPrv_dot - Qv_RCA)/Cmicro_RCA


    return np.array([Vra_dot, Vrv_dot, Vla_dot, Vlv_dot, Qtv_dot, Qpv_dot, Qmv_dot, Qav_dot, eps_tv_dot, eps_pv_dot, eps_mv_dot, eps_av_dot, Pa_syst_dot, Pv_syst_dot, Pa_pulm_dot, Pv_pulm_dot, Pa_LAD_dot, Pv_LAD_dot, Pa_LCX_dot, Pv_LCX_dot, Pa_RCA_dot, Pv_RCA_dot])

