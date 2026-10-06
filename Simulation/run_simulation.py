import numpy as np
from scipy.integrate import solve_ivp, simpson
from shapely.geometry import Polygon
from numba import njit
import pandas as pd
import multiprocessing as mp
import json
from copy import deepcopy

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





def cycle(
        Emax_ra, Emin_ra, m1_ra, m2_ra, V0ra,
        Emax_rv, Emin_rv, m1_rv, m2_rv, V0rv,
        Emax_la, Emin_la, m1_la, m2_la, V0la,
        Emax_lv, Emin_lv, m1_lv, m2_lv, V0lv,
        tau1_ra, tau2_ra, k_ra, tau1_rv, tau2_rv, k_rv,
        tau1_la, tau2_la, k_la, tau1_lv, tau2_lv, k_lv,
        tau1_LAD, tau2_LAD, k_LAD, tau1_LCX, tau2_LCX, k_LCX, tau1_RCA, tau2_RCA, k_RCA,
        Amax_tv, Amin_tv, Kvo_tv, Kvc_tv,
        Amax_pv, Amin_pv, Kvo_pv, Kvc_pv,
        Amax_mv, Amin_mv, Kvo_mv, Kvc_mv,
        Amax_av, Amin_av, Kvo_av, Kvc_av,
        HR, density, Rao, Ca_syst, Ra_syst, Cv_syst, Rv_syst, Rpa, Ca_pulm, Ra_pulm, Cv_pulm, Rv_pulm,
        Ra_LAD, Ra_LCX, Ra_RCA, Rmicro_max_LAD, Rmicro_min_LAD, Rmicro_max_LCX, Rmicro_min_LCX, Rmicro_max_RCA, Rmicro_min_RCA, Rv_LAD, Rv_LCX, Rv_RCA, Ca_LAD, Ca_LCX, Ca_RCA, Cmicro_LAD, Cmicro_LCX, Cmicro_RCA, hypertrophy_factor,
        icVra, icVrv, icVla, icVlv, icQtv, icQpv, icQmv, icQav, iceps_tv, iceps_pv, iceps_mv, iceps_av, icPa_syst, icPv_syst, icPa_pulm, icPv_pulm,
        icPa_LAD, icPv_LAD, icPa_LCX, icPv_LCX, icPa_RCA, icPv_RCA
    ):

    # conversion factor
    conversion_factor = 1333.22


    

    #-----solving the system of ODEs-----#
    periode = 60/HR
    u0 = [icVra, icVrv, icVla, icVlv, icQtv, icQpv, icQmv, icQav, iceps_tv, iceps_pv, iceps_mv, iceps_av, icPa_syst, icPv_syst, icPa_pulm, icPv_pulm, icPa_LAD, icPv_LAD, icPa_LCX, icPv_LCX, icPa_RCA, icPv_RCA]
    tspan = (0, periode)
    t_to_solve = np.linspace(tspan[0], tspan[-1], 2000)
    delta = 0.85*periode


    params_flat = np.array([
    conversion_factor, delta, periode,
    tau1_ra, tau2_ra, m1_ra, m2_ra, k_ra, Emin_ra, tau1_rv, tau2_rv, m1_rv, m2_rv, k_rv, Emin_rv,
    tau1_la, tau2_la, m1_la, m2_la, k_la, Emin_la, tau1_lv, tau2_lv, m1_lv, m2_lv, k_lv, Emin_lv,
    V0ra, V0rv, V0la, V0lv, density,
    Rao, Ra_syst, Rv_syst, Rpa, Ra_pulm, Rv_pulm, Ca_syst, Cv_syst, Ca_pulm, Cv_pulm,
    Amax_tv, Amin_tv, Amax_pv, Amin_pv, Amax_mv, Amin_mv, Amax_av, Amin_av,
    Kvo_tv, Kvc_tv, Kvo_pv, Kvc_pv, Kvo_mv, Kvc_mv, Kvo_av, Kvc_av,
    Ra_LAD, Ra_LCX, Ra_RCA, Rv_LAD, Rv_LCX, Rv_RCA, Ca_LAD, Ca_LCX, Ca_RCA, Cmicro_LAD, Cmicro_LCX, Cmicro_RCA,
    tau1_LAD, tau2_LAD, m1_lv, m2_lv, k_LAD, Rmicro_min_LAD, 
    tau1_LCX, tau2_LCX, m1_lv, m2_lv, k_LCX, Rmicro_min_LCX, 
    tau1_RCA, tau2_RCA, m1_rv, m2_rv, k_RCA, Rmicro_min_RCA
    ], dtype=np.float64)

    solution = solve_ivp(system_odes_scipy_faster, t_span=tspan, y0=u0, t_eval=t_to_solve, method="BDF", rtol=1e-6,atol=1e-9,args=(params_flat,))

    Vra = solution.y[0]
    Vrv = solution.y[1]
    Vla = solution.y[2]
    Vlv = solution.y[3]
    Qtv = solution.y[4]
    Qpv = solution.y[5]
    Qmv = solution.y[6]
    Qav = solution.y[7]
    eps_tv = solution.y[8]
    eps_pv = solution.y[9]
    eps_mv = solution.y[10]
    eps_av = solution.y[11]
    Pa_syst = solution.y[12]/conversion_factor
    Pv_syst = solution.y[13]/conversion_factor
    Pa_pulm = solution.y[14]/conversion_factor
    Pv_pulm = solution.y[15]/conversion_factor
    Pa_LAD = solution.y[16]/conversion_factor
    Pv_LAD = solution.y[17]/conversion_factor
    Pa_LCX = solution.y[18]/conversion_factor
    Pv_LCX = solution.y[19]/conversion_factor
    Pa_RCA = solution.y[20]/conversion_factor
    Pv_RCA = solution.y[21]/conversion_factor

    


    #-----returning results-----#
    return {
        "Pa_syst_max": max(Pa_syst),
        "Pa_syst_min": min(Pa_syst),
        "Vav": simpson(Qav, x=t_to_solve),
        "Vmv": simpson(Qmv, x=t_to_solve),

        "t_to_solve": t_to_solve,
        "Vra": Vra,
        "Vrv": Vrv,
        "Vla": Vla,
        "Vlv": Vlv,
        "Qtv": Qtv,
        "Qpv": Qpv,
        "Qmv": Qmv,
        "Qav": Qav,
        "Pa_syst": Pa_syst,
        "Pv_syst": Pv_syst,
        "Pa_pulm": Pa_pulm,
        "Pv_pulm": Pv_pulm,
        "Pa_LAD": Pa_LAD,
        "Pv_LAD": Pv_LAD,
        "Pa_LCX": Pa_LCX,
        "Pv_LCX": Pv_LCX,
        "Pa_RCA": Pa_RCA,
        "Pv_RCA": Pv_RCA,

        "icVra": Vra[-1],
        "icVrv": Vrv[-1],
        "icVla": Vla[-1],
        "icVlv": Vlv[-1],
        "icQtv": Qtv[-1],
        "icQpv": Qpv[-1],
        "icQmv": Qmv[-1],
        "icQav": Qav[-1],
        "iceps_tv": eps_tv[-1],
        "iceps_pv": eps_pv[-1],
        "iceps_mv": eps_mv[-1],
        "iceps_av": eps_av[-1],
        "icPa_syst": Pa_syst[-1]*conversion_factor,
        "icPv_syst": Pv_syst[-1]*conversion_factor,
        "icPa_pulm": Pa_pulm[-1]*conversion_factor,
        "icPv_pulm": Pv_pulm[-1]*conversion_factor,
        "icPa_LAD": Pa_LAD[-1]*conversion_factor,
        "icPv_LAD": Pv_LAD[-1]*conversion_factor,
        "icPa_LCX": Pa_LCX[-1]*conversion_factor,
        "icPv_LCX": Pv_LCX[-1]*conversion_factor,
        "icPa_RCA": Pa_RCA[-1]*conversion_factor,
        "icPv_RCA": Pv_RCA[-1]*conversion_factor
        }
    
        



def multiple_cycles(parameters_input, number_of_cycles, criterion=0.005):
    #-----load parameters & conversion-----#
    conversion_factor = 1333.22
    # chambers
    chambers_params = parameters_input["chambers"]
    Emax_ra = chambers_params["Emax_ra"]
    Emin_ra = chambers_params["Emin_ra"]
    m1_ra = chambers_params["m1_ra"]
    m2_ra = chambers_params["m2_ra"]
    V0ra = chambers_params["V0ra"]
    Emax_rv = chambers_params["Emax_rv"]
    Emin_rv = chambers_params["Emin_rv"]
    m1_rv = chambers_params["m1_rv"]
    m2_rv = chambers_params["m2_rv"]
    V0rv = chambers_params["V0rv"]
    Emax_la = chambers_params["Emax_la"]
    Emin_la = chambers_params["Emin_la"]
    m1_la = chambers_params["m1_la"]
    m2_la = chambers_params["m2_la"]
    V0la = chambers_params["V0la"]
    Emax_lv = chambers_params["Emax_lv"]
    Emin_lv = chambers_params["Emin_lv"]
    m1_lv = chambers_params["m1_lv"]
    m2_lv = chambers_params["m2_lv"]
    V0lv = chambers_params["V0lv"]

    # valves
    valves_params = parameters_input["valves"]
    Amax_tv = valves_params["Amax_tv"]
    Amin_tv = valves_params["Amin_tv"]
    Kvo_tv = valves_params["Kvo_tv"] 
    Kvc_tv = valves_params["Kvc_tv"] 
    Amax_pv = valves_params["Amax_pv"]
    Amin_pv = valves_params["Amin_pv"]
    Kvo_pv = valves_params["Kvo_pv"] 
    Kvc_pv = valves_params["Kvc_pv"] 
    Amax_mv = valves_params["Amax_mv"]
    Amin_mv = valves_params["Amin_mv"]
    Kvo_mv = valves_params["Kvo_mv"] 
    Kvc_mv = valves_params["Kvc_mv"] 
    Amax_av = valves_params["Amax_av"]
    Amin_av = valves_params["Amin_av"]
    Kvo_av = valves_params["Kvo_av"] 
    Kvc_av = valves_params["Kvc_av"] 

    # systemic and pulmonarycirculation
    circulation_params = parameters_input["circulation"]
    HR = circulation_params["HR"]
    density = circulation_params["density"]
    Rao = circulation_params["Rao"]*conversion_factor
    Ca_syst = circulation_params["Ca_syst"]/conversion_factor
    Ra_syst = circulation_params["Ra_syst"]*conversion_factor
    Cv_syst = circulation_params["Cv_syst"]/conversion_factor
    Rv_syst = circulation_params["Rv_syst"]*conversion_factor
    Rpa = circulation_params["Rpa"]*conversion_factor
    Ca_pulm = circulation_params["Ca_pulm"]/conversion_factor
    Ra_pulm = circulation_params["Ra_pulm"]*conversion_factor
    Cv_pulm = circulation_params["Cv_pulm"]/conversion_factor
    Rv_pulm = circulation_params["Rv_pulm"]*conversion_factor

    # coronary circulation
    coronary_params = parameters_input["CBF_model"]
    Ra_LAD = coronary_params["Ra_LAD"]*conversion_factor
    Ra_LCX = coronary_params["Ra_LCX"]*conversion_factor
    Ra_RCA = coronary_params["Ra_RCA"]*conversion_factor
    Rmicro_max_LAD = coronary_params["Rmicro_max_LAD"]*conversion_factor
    Rmicro_max_LCX = coronary_params["Rmicro_max_LCX"]*conversion_factor
    Rmicro_max_RCA = coronary_params["Rmicro_max_RCA"]*conversion_factor
    Rv_LAD = coronary_params["Rv_LAD"]*conversion_factor
    Rv_LCX = coronary_params["Rv_LCX"]*conversion_factor
    Rv_RCA = coronary_params["Rv_RCA"]*conversion_factor
    Ca_LAD = coronary_params["Ca_LAD"]/conversion_factor
    Ca_LCX = coronary_params["Ca_LCX"]/conversion_factor
    Ca_RCA = coronary_params["Ca_RCA"]/conversion_factor
    Cmicro_LAD = coronary_params["Cmicro_LAD"]/conversion_factor
    Cmicro_LCX = coronary_params["Cmicro_LCX"]/conversion_factor
    Cmicro_RCA = coronary_params["Cmicro_RCA"]/conversion_factor
    hypertrophy_factor = coronary_params["hypertrophy_factor"]

    # initial conditions
    initial_conditions = parameters_input["initial_conditions"]
    icVra = initial_conditions["icVra"]
    icVrv = initial_conditions["icVrv"]
    icVla = initial_conditions["icVla"]
    icVlv = initial_conditions["icVlv"]
    icQtv = initial_conditions["icQtv"]
    icQpv = initial_conditions["icQpv"]
    icQmv = initial_conditions["icQmv"]
    icQav = initial_conditions["icQav"]
    iceps_tv = initial_conditions["iceps_tv"]
    iceps_pv = initial_conditions["iceps_pv"]
    iceps_mv = initial_conditions["iceps_mv"]
    iceps_av = initial_conditions["iceps_av"]
    icPa_syst = initial_conditions["icPa_syst"]*conversion_factor
    icPv_syst = initial_conditions["icPv_syst"]*conversion_factor
    icPa_pulm = initial_conditions["icPa_pulm"]*conversion_factor
    icPv_pulm = initial_conditions["icPv_pulm"]*conversion_factor
    icPa_LAD = initial_conditions["icPa_LAD"]*conversion_factor
    icPv_LAD = initial_conditions["icPv_LAD"]*conversion_factor
    icPa_LCX = initial_conditions["icPa_LCX"]*conversion_factor
    icPv_LCX = initial_conditions["icPv_LCX"]*conversion_factor
    icPa_RCA = initial_conditions["icPa_RCA"]*conversion_factor
    icPv_RCA = initial_conditions["icPv_RCA"]*conversion_factor

    # set Rmicro_min before hypertrophy influence
    Rmicro_min_LAD = Rmicro_max_LAD*0.7
    Rmicro_min_LCX = Rmicro_max_LCX*0.7
    Rmicro_min_RCA = Rmicro_max_RCA*0.75

    # LV hypertrophy influence for LAD and LCX
    Rmicro_max_LAD = Rmicro_max_LAD*hypertrophy_factor
    Rmicro_max_LCX = Rmicro_max_LCX*hypertrophy_factor

    # activation function parameters
    periode = 60/HR
    tau1_ra, tau2_ra, k_ra = initialize_activation_functions("atrium", periode, m1_ra, m2_ra, Emax_ra, Emin_ra)
    tau1_rv, tau2_rv, k_rv = initialize_activation_functions("ventricle", periode, m1_rv, m2_rv, Emax_rv, Emin_rv)
    tau1_la, tau2_la, k_la = initialize_activation_functions("atrium", periode, m1_la, m2_la, Emax_la, Emin_la)
    tau1_lv, tau2_lv, k_lv = initialize_activation_functions("ventricle", periode, m1_lv, m2_lv, Emax_lv, Emin_lv)
    _, tau1_LAD, tau2_LAD, k_LAD = Et_Rmicro(0.0, m1_lv, m2_lv, Rmicro_min_LAD, Rmicro_max_LAD, periode, rt_all=True)
    _, tau1_LCX, tau2_LCX, k_LCX = Et_Rmicro(0.0, m1_lv, m2_lv, Rmicro_min_LCX, Rmicro_max_LCX, periode, rt_all=True)
    _, tau1_RCA, tau2_RCA, k_RCA = Et_Rmicro(0.0, m1_rv, m2_rv, Rmicro_min_RCA, Rmicro_max_RCA, periode, rt_all=True)

    #-----initialize-----#
    list_of_results = []
    stability_counter = 0
    prev_metrics = None

    # compute
    for i in range(number_of_cycles):

        results = cycle(
        Emax_ra, Emin_ra, m1_ra, m2_ra, V0ra,
        Emax_rv, Emin_rv, m1_rv, m2_rv, V0rv,
        Emax_la, Emin_la, m1_la, m2_la, V0la,
        Emax_lv, Emin_lv, m1_lv, m2_lv, V0lv,
        tau1_ra, tau2_ra, k_ra, tau1_rv, tau2_rv, k_rv,
        tau1_la, tau2_la, k_la, tau1_lv, tau2_lv, k_lv,
        tau1_LAD, tau2_LAD, k_LAD, tau1_LCX, tau2_LCX, k_LCX, tau1_RCA, tau2_RCA, k_RCA,
        Amax_tv, Amin_tv, Kvo_tv, Kvc_tv,
        Amax_pv, Amin_pv, Kvo_pv, Kvc_pv,
        Amax_mv, Amin_mv, Kvo_mv, Kvc_mv,
        Amax_av, Amin_av, Kvo_av, Kvc_av,
        HR, density, Rao, Ca_syst, Ra_syst, Cv_syst, Rv_syst, Rpa, Ca_pulm, Ra_pulm, Cv_pulm, Rv_pulm,
        Ra_LAD, Ra_LCX, Ra_RCA, Rmicro_max_LAD, Rmicro_min_LAD, Rmicro_max_LCX, Rmicro_min_LCX, Rmicro_max_RCA, Rmicro_min_RCA, Rv_LAD, Rv_LCX, Rv_RCA, Ca_LAD, Ca_LCX, Ca_RCA, Cmicro_LAD, Cmicro_LCX, Cmicro_RCA, hypertrophy_factor,
        icVra, icVrv, icVla, icVlv, icQtv, icQpv, icQmv, icQav, iceps_tv, iceps_pv, iceps_mv, iceps_av, icPa_syst, icPv_syst, icPa_pulm, icPv_pulm,
        icPa_LAD, icPv_LAD, icPa_LCX, icPv_LCX, icPa_RCA, icPv_RCA
        )

        # change initial values
        icVra = results["icVra"]
        icVrv = results["icVrv"]
        icVla = results["icVla"]
        icVlv = results["icVlv"]
        icQtv = results["icQtv"]
        icQpv = results["icQpv"]
        icQmv = results["icQmv"]
        icQav = results["icQav"]
        iceps_tv = results["iceps_tv"]
        iceps_pv = results["iceps_pv"]
        iceps_mv = results["iceps_mv"]
        iceps_av = results["iceps_av"]
        icPa_syst = results["icPa_syst"]
        icPv_syst = results["icPv_syst"]
        icPa_pulm = results["icPa_pulm"]
        icPv_pulm = results["icPv_pulm"]
        icPa_LAD = results["icPa_LAD"]
        icPv_LAD = results["icPv_LAD"]
        icPa_LCX = results["icPa_LCX"]
        icPv_LCX = results["icPv_LCX"]
        icPa_RCA = results["icPa_RCA"]
        icPv_RCA = results["icPv_RCA"]

        
        current_metrics = {
            "Pa_syst_max": results["Pa_syst_max"],
            "Pa_syst_min": results["Pa_syst_min"],
            "Vav": results["Vav"],
            "Vmv": results["Vmv"]
        }
    
        # check stability - 5 times
        if i > 1:
            if (
                abs(current_metrics["Pa_syst_max"] - prev_metrics["Pa_syst_max"]) < criterion
                and abs(current_metrics["Pa_syst_min"] - prev_metrics["Pa_syst_min"]) < criterion
                and abs(current_metrics["Vav"] - prev_metrics["Vav"]) < criterion
                and abs(current_metrics["Vav"] - current_metrics["Vmv"]) < criterion
                ):
                if stability_counter > 4:
                    stability_reached = True
                    break
                else:
                    stability_counter += 1
            else:
                stability_counter = 0

        prev_metrics = current_metrics

    if stability_counter <= 4:
        stability_reached = False

    # unpack
    t_to_solve = results["t_to_solve"]
    Vra = results["Vra"]
    Vrv = results["Vrv"]
    Vla = results["Vla"]
    Vlv = results["Vlv"]
    Qtv = results["Qtv"]
    Qpv = results["Qpv"]
    Qmv = results["Qmv"]
    Qav = results["Qav"]
    Pa_syst = results["Pa_syst"]
    Pv_syst = results["Pv_syst"]
    Pa_pulm = results["Pa_pulm"]
    Pv_pulm = results["Pv_pulm"]
    Pa_LAD = results["Pa_LAD"]
    Pv_LAD = results["Pv_LAD"]
    Pa_LCX = results["Pa_LCX"]
    Pv_LCX = results["Pv_LCX"]
    Pa_RCA = results["Pa_RCA"]
    Pv_RCA = results["Pv_RCA"]

    #-----calculate output values-----#

    Era = activation_function(chamber_type="atrium", x=t_to_solve, tau1=tau1_ra, tau2=tau2_ra, m1=m1_ra, m2=m2_ra, k=k_ra, Emin=Emin_ra, periode=periode, delta=0.85*periode)
    Erv = activation_function(chamber_type="ventricle", x=t_to_solve, tau1=tau1_rv, tau2=tau2_rv, m1=m1_rv, m2=m2_rv, k=k_rv, Emin=Emin_rv)
    Ela = activation_function(chamber_type="atrium", x=t_to_solve, tau1=tau1_la, tau2=tau2_la, m1=m1_la, m2=m2_la, k=k_la, Emin=Emin_la, periode=periode, delta=0.85*periode)
    Elv = activation_function(chamber_type="ventricle", x=t_to_solve, tau1=tau1_lv, tau2=tau2_lv, m1=m1_lv, m2=m2_lv, k=k_lv, Emin=Emin_lv)

    Pra = Era*(Vra - V0ra)
    Prv = Erv*(Vrv - V0rv)
    Pla = Ela*(Vla - V0la)
    Plv = Elv*(Vlv - V0lv)

    Pao = ((Qav + Pa_syst*conversion_factor/Rao + Pa_LAD*conversion_factor/Ra_LAD + Pa_LCX*conversion_factor/Ra_LCX + Pa_RCA*conversion_factor/Ra_RCA)/(1/Rao + 1/Ra_LAD + 1/Ra_LCX + 1/Ra_RCA))/conversion_factor
    Ppt = (Pa_pulm*conversion_factor + Rpa*Qpv)/conversion_factor
    
    Pao_max = max(Pao)
    Pao_min = min(Pao)
    Pao_mean = np.mean(Pao)

    Ppt_max = max(Ppt)
    Ppt_min = min(Ppt)
    Ppt_mean = np.mean(Ppt)

    Pa_syst_max = max(Pa_syst)
    Pa_syst_min = min(Pa_syst)
    Pa_syst_mean = np.mean(Pa_syst)

    Pa_pulm_max = max(Pa_pulm)
    Pa_pulm_min = min(Pa_pulm)
    Pa_pulm_mean = np.mean(Pa_pulm)

    Pv_syst_max = max(Pv_syst)
    Pv_syst_min = min(Pv_syst)
    Pv_syst_mean = np.mean(Pv_syst)

    Pv_pulm_max = max(Pv_pulm)
    Pv_pulm_min = min(Pv_pulm)
    Pv_pulm_mean = np.mean(Pv_pulm)

    Pa_LAD_max = max(Pa_LAD)
    Pa_LAD_min = min(Pa_LAD)
    Pa_LAD_mean = np.mean(Pa_LAD)
    Pa_LCX_max = max(Pa_LCX)
    Pa_LCX_min = min(Pa_LCX)
    Pa_LCX_mean = np.mean(Pa_LCX)
    Pa_RCA_max = max(Pa_RCA)
    Pa_RCA_min = min(Pa_RCA)
    Pa_RCA_mean = np.mean(Pa_RCA)

    Pv_LAD_max = max(Pv_LAD)
    Pv_LAD_min = min(Pv_LAD)
    Pv_LAD_mean = np.mean(Pv_LAD)
    Pv_LCX_max = max(Pv_LCX)
    Pv_LCX_min = min(Pv_LCX)
    Pv_LCX_mean = np.mean(Pv_LCX)
    Pv_RCA_max = max(Pv_RCA)
    Pv_RCA_min = min(Pv_RCA)
    Pv_RCA_mean = np.mean(Pv_RCA)

    Pra_max = max(Pra)
    Pra_min = min(Pra)
    Pra_mean = np.mean(Pra)
    Prv_max = max(Prv)
    Prv_min = min(Prv)
    Prv_mean = np.mean(Prv)
    Pla_max = max(Pla)
    Pla_min = min(Pla)
    Pla_mean = np.mean(Pla)
    Plv_max = max(Plv)
    Plv_min = min(Plv)
    Plv_mean = np.mean(Plv)
    EDP_rv = Prv[-1]
    EDP_lv = Plv[-1]

    Vra_max = max(Vra)
    Vra_min = min(Vra)
    Vra_mean = np.mean(Vra)
    Vrv_max = max(Vrv)
    Vrv_min = min(Vrv)
    Vrv_mean = np.mean(Vrv)
    Vla_max = max(Vla)
    Vla_min = min(Vla)
    Vla_mean = np.mean(Vla)
    Vlv_max = max(Vlv)
    Vlv_min = min(Vlv)
    Vlv_mean = np.mean(Vlv)

    Vtv = simpson(Qtv, x=t_to_solve)
    Vpv = simpson(Qpv, x=t_to_solve)
    Vmv = simpson(Qmv, x=t_to_solve)
    Vav = simpson(Qav, x=t_to_solve)
    Qtv_max = max(Qtv)
    Qpv_max = max(Qpv)
    Qmv_max = max(Qmv)
    Qav_max = max(Qav)

    start_systole_lv_index = np.where(Qav>0.1)[0][0]
    end_systole_lv_index = np.where(Qav>0.1)[0][-1]
    start_systole_rv_index = np.where(Qpv>0.1)[0][0]
    end_systole_rv_index = np.where(Qpv>0.1)[0][-1]

    SW_rv_polygon = Polygon(zip(Vrv, Prv))
    SW_rv = SW_rv_polygon.area*(133.332*10**-6)
    SW_lv_polygon = Polygon(zip(Vlv, Plv))
    SW_lv = SW_lv_polygon.area*(133.332*10**-6)
    Pes_lv = Plv[end_systole_lv_index]
    Pes_rv = Prv[end_systole_rv_index]
    PE_lv = (0.5*Pes_lv*(Vlv_min-V0lv))*(133.332*10**-6)
    PE_rv = (0.5*Pes_rv*(Vrv_min-V0rv))*(133.332*10**-6)
    PVA_lv = SW_lv + PE_lv
    PVA_rv = SW_rv + PE_rv
    VO2_lv = 2.46*PVA_lv + 0.7
    VO2_rv = 2.46*PVA_rv + 0.2
    VO2_total = VO2_lv + VO2_rv

    Rmicro_LAD = Et_Rmicro(t_to_solve, m1_lv, m2_lv, Rmicro_min_LAD, Rmicro_max_LAD, periode)
    Rmicro_LCX = Et_Rmicro(t_to_solve, m1_lv, m2_lv, Rmicro_min_LCX, Rmicro_max_LCX, periode)
    Rmicro_RCA = Et_Rmicro(t_to_solve, m1_rv, m2_rv, Rmicro_min_RCA, Rmicro_max_RCA, periode)
    Qmicro_LAD = (Pa_LAD- Pv_LAD)/(Rmicro_LAD/conversion_factor)
    Qmicro_LCX = (Pa_LCX - Pv_LCX)/(Rmicro_LCX/conversion_factor)
    Qmicro_RCA = (Pa_RCA - Pv_RCA)/(Rmicro_RCA/conversion_factor)
    Qmicro_CBF_LAD = simpson(Qmicro_LAD, x=t_to_solve)
    Qmicro_CBF_LCX = simpson(Qmicro_LCX, x=t_to_solve)
    Qmicro_CBF_RCA = simpson(Qmicro_RCA, x=t_to_solve)
    Qmicro_CBF_total = simpson(Qmicro_LAD+Qmicro_LCX+Qmicro_RCA, x=t_to_solve)
    Ein = Qmicro_CBF_total*20*0.2*0.95
    Ein_LAD = Qmicro_CBF_LAD*20*0.2*0.95
    Ein_LCX = Qmicro_CBF_LCX*20*0.2*0.95
    Ein_RCA = Qmicro_CBF_RCA*20*0.2*0.95
    
    CO_lv = Vav * HR
    CO_rv = Vpv * HR
    LVET = t_to_solve[end_systole_lv_index] - t_to_solve[start_systole_lv_index]
    RVET = t_to_solve[end_systole_rv_index] - t_to_solve[start_systole_rv_index]
    Ees_lv = Pes_lv/(Vlv_min-V0lv)
    Ees_rv = Pes_rv/(Vrv_min-V0rv)
    Ea_lv = Pes_lv/Vav
    Ea_rv = Pes_rv/Vpv
    VAC_lv = Ea_lv/Ees_lv
    VAC_rv = Ea_rv/Ees_rv

    d = {
        "stability_reached": stability_reached,
        "stability_reached_cycle": i+1,

        "Pao_max": Pao_max,
        "Pao_min": Pao_min,
        "Pao_mean": Pao_mean,
        "Ppt_max": Ppt_max,
        "Ppt_min": Ppt_min,
        "Ppt_mean": Ppt_mean,
        "Pa_syst_max": Pa_syst_max,
        "Pa_syst_min": Pa_syst_min,
        "Pa_syst_mean": Pa_syst_mean,
        "Pv_syst_max": Pv_syst_max,
        "Pv_syst_min": Pv_syst_min,
        "Pv_syst_mean": Pv_syst_mean,
        "Pa_pulm_max": Pa_pulm_max,
        "Pa_pulm_min": Pa_pulm_min,
        "Pa_pulm_mean": Pa_pulm_mean,
        "Pv_pulm_max": Pv_pulm_max,
        "Pv_pulm_min": Pv_pulm_min,
        "Pv_pulm_mean": Pv_pulm_mean,
        "Pa_LAD_max": Pa_LAD_max,
        "Pa_LAD_min": Pa_LAD_min,
        "Pa_LAD_mean": Pa_LAD_mean,
        "Pa_LCX_max": Pa_LCX_max,
        "Pa_LCX_min": Pa_LCX_min,
        "Pa_LCX_mean": Pa_LCX_mean,
        "Pa_RCA_max": Pa_RCA_max,
        "Pa_RCA_min": Pa_RCA_min,
        "Pa_RCA_mean": Pa_RCA_mean,
        "Pv_LAD_max": Pv_LAD_max,
        "Pv_LAD_min": Pv_LAD_min,
        "Pv_LAD_mean": Pv_LAD_mean,
        "Pv_LCX_max": Pv_LCX_max,
        "Pv_LCX_min": Pv_LCX_min,
        "Pv_LCX_mean": Pv_LCX_mean,
        "Pv_RCA_max": Pv_RCA_max,
        "Pv_RCA_min": Pv_RCA_min,
        "Pv_RCA_mean": Pv_RCA_mean,
        "Pra_max": Pra_max,
        "Pra_min": Pra_min,
        "Pra_mean": Pra_mean,
        "Prv_max": Prv_max,
        "Prv_min": Prv_min,
        "Prv_mean": Prv_mean,
        "Pes_rv": Pes_rv,
        "Pla_max": Pla_max,
        "Pla_min": Pla_min,
        "Pla_mean": Pla_mean,
        "Plv_max": Plv_max,
        "Plv_min": Plv_min,
        "Plv_mean": Plv_mean,
        "Pes_lv": Pes_lv,
        "EDP_rv": EDP_rv,
        "EDP_lv": EDP_lv,
        "Vra_max": Vra_max,
        "Vra_min": Vra_min,
        "Vra_mean": Vra_mean,
        "Vrv_max": Vrv_max,
        "Vrv_min": Vrv_min,
        "Vrv_mean": Vrv_mean,
        "Vla_max": Vla_max,
        "Vla_min": Vla_min,
        "Vla_mean": Vla_mean,
        "Vlv_max": Vlv_max,
        "Vlv_min": Vlv_min,
        "Vlv_mean": Vlv_mean,
        "Vtv": Vtv,
        "Vpv": Vpv,
        "Vmv": Vmv,
        "Vav": Vav,
        "Qtv_max": Qtv_max,
        "Qpv_max": Qpv_max,
        "Qmv_max": Qmv_max,
        "Qav_max": Qav_max,
        "SW_rv": SW_rv,
        "SW_lv": SW_lv,
        "PE_rv": PE_rv,
        "PE_lv": PE_lv,
        "PVA_rv": PVA_rv,
        "PVA_lv": PVA_lv,
        "VO2_rv": VO2_rv,
        "VO2_lv": VO2_lv,
        "VO2_total": VO2_total,
        "Ein": Ein,
        "Ein_LAD": Ein_LAD,
        "Ein_LCX": Ein_LCX,
        "Ein_RCA": Ein_RCA,
        "CO_rv": CO_rv,
        "CO_lv": CO_lv,
        "RVET": RVET,
        "LVET": LVET,
        "Ees_rv": Ees_rv,
        "Ees_lv": Ees_lv,
        "Ea_rv": Ea_rv,
        "Ea_lv": Ea_lv,
        "VAC_rv": VAC_rv,
        "VAC_lv": VAC_lv,
        "DeltaE": Ein-VO2_total,
        "Patient_ID": parameters_input["Patient_ID"]
    }
    z = np.array(list(d.values()))
    return z

def worker_simulation(args):
    inx, row, base_params = args
    
    parameters_input = deepcopy(base_params)

    parameters_input["chambers"]["Emax_lv"] = row["Emax_lv"]
    parameters_input["chambers"]["Emin_lv"] = row["Emin_lv"]
    
    parameters_input["valves"]["Amax_av"] = row["Amax_av"]

    parameters_input["circulation"]["HR"] = row["HR"]
    parameters_input["circulation"]["Ra_syst"] = row["Ra_syst"]
    parameters_input["circulation"]["Rv_syst"] = row["Rv_syst"]
    parameters_input["circulation"]["Cv_syst"] = row["Cv_syst"]
    
    parameters_input["CBF_model"]["Rmicro_max_LAD"] = row["Rmicro_max_LAD"]
    parameters_input["CBF_model"]["Rmicro_max_LCX"] = row["Rmicro_max_LCX"]
    
    parameters_input["initial_conditions"]["icPv_syst"] = row["icPv_syst"]

    if "k_hyp" in row.keys():
        parameters_input["CBF_model"]["hypertrophy_factor"] = row["k_hyp"]

    parameters_input["Patient_ID"] = row["Patient_ID"]

    try:
        results = multiple_cycles(parameters_input, number_of_cycles=800, criterion=0.005)
        return inx, results
    except:
        return inx, None


def main_fc(file_name_universal_input, file_name_specific_input, file_name_output):
    with open(f"{file_name_universal_input}.json", "r") as f: 
        base_params = json.load(f)

    input_df = pd.read_csv(f"{file_name_specific_input}.csv") 
    Y = np.zeros((len(input_df), 99))

    tasks = [(inx, input_df.iloc[inx], base_params) for inx in range(len(input_df))] 

    num_cores = 48 
    print(f"Starting parallel processing with {num_cores} cores for {file_name_specific_input}", flush=True)

    completed = 0
    with mp.Pool(processes=num_cores) as pool:
        for inx, results in pool.imap_unordered(worker_simulation, tasks):
            if results is not None:
                Y[inx,:] = results
            
            completed += 1
            if completed % 25000 == 0:
                print(f"Completed {completed} out of {len(input_df)} simulations.", flush=True)
                np.save(f'{file_name_output}.npy', Y)

    output = {
        "stability_reached": Y[:,0],
        "stability_reached_cycle": Y[:,1],
        "Pao_max": Y[:,2],
        "Pao_min": Y[:,3],
        "Pao_mean": Y[:,4],
        "Ppt_max": Y[:,5],
        "Ppt_min": Y[:,6],
        "Ppt_mean": Y[:,7],
        "Pa_syst_max": Y[:,8],
        "Pa_syst_min": Y[:,9],
        "Pa_syst_mean": Y[:,10],
        "Pv_syst_max": Y[:,11],
        "Pv_syst_min": Y[:,12],
        "Pv_syst_mean": Y[:,13],
        "Pa_pulm_max": Y[:,14],
        "Pa_pulm_min": Y[:,15],
        "Pa_pulm_mean": Y[:,16],
        "Pv_pulm_max": Y[:,17],
        "Pv_pulm_min": Y[:,18],
        "Pv_pulm_mean": Y[:,19],
        "Pa_LAD_max": Y[:,20],
        "Pa_LAD_min": Y[:,21],
        "Pa_LAD_mean": Y[:,22],
        "Pa_LCX_max": Y[:,23],
        "Pa_LCX_min": Y[:,24],
        "Pa_LCX_mean": Y[:,25],
        "Pa_RCA_max": Y[:,26],
        "Pa_RCA_min": Y[:,27],
        "Pa_RCA_mean": Y[:,28],
        "Pv_LAD_max": Y[:,29],
        "Pv_LAD_min": Y[:,30],
        "Pv_LAD_mean": Y[:,31],
        "Pv_LCX_max": Y[:,32],
        "Pv_LCX_min": Y[:,33],
        "Pv_LCX_mean": Y[:,34],
        "Pv_RCA_max": Y[:,35],
        "Pv_RCA_min": Y[:,36],
        "Pv_RCA_mean": Y[:,37],
        "Pra_max": Y[:,38],
        "Pra_min": Y[:,39],
        "Pra_mean": Y[:,40],
        "Prv_max": Y[:,41],
        "Prv_min": Y[:,42],
        "Prv_mean": Y[:,43],
        "Pes_rv": Y[:,44],
        "Pla_max": Y[:,45],
        "Pla_min": Y[:,46],
        "Pla_mean": Y[:,47],
        "Plv_max": Y[:,48],
        "Plv_min": Y[:,49],
        "Plv_mean": Y[:,50],
        "Pes_lv": Y[:,51],
        "EDP_rv": Y[:,52],
        "EDP_lv": Y[:,53],
        "Vra_max": Y[:,54],
        "Vra_min": Y[:,55],
        "Vra_mean": Y[:,56],
        "Vrv_max": Y[:,57],
        "Vrv_min": Y[:,58],
        "Vrv_mean": Y[:,59],
        "Vla_max": Y[:,60],
        "Vla_min": Y[:,61],
        "Vla_mean": Y[:,62],
        "Vlv_max": Y[:,63],
        "Vlv_min": Y[:,64],
        "Vlv_mean": Y[:,65],
        "Vtv": Y[:,66],
        "Vpv": Y[:,67],
        "Vmv": Y[:,68],
        "Vav": Y[:,69],
        "Qtv_max": Y[:,70],
        "Qpv_max": Y[:,71],
        "Qmv_max": Y[:,72],
        "Qav_max": Y[:,73],
        "SW_rv": Y[:,74],
        "SW_lv": Y[:,75],
        "PE_rv": Y[:,76],
        "PE_lv": Y[:,77],
        "PVA_rv": Y[:,78],
        "PVA_lv": Y[:,79],
        "VO2_rv": Y[:,80],
        "VO2_lv": Y[:,81],
        "VO2_total": Y[:,82],
        "Ein": Y[:,83],
        "Ein_LAD": Y[:,84],
        "Ein_LCX": Y[:,85],
        "Ein_RCA": Y[:,86],
        "CO_rv": Y[:,87],
        "CO_lv": Y[:,88],
        "RVET": Y[:,89],
        "LVET": Y[:,90],
        "Ees_rv": Y[:,91],
        "Ees_lv": Y[:,92],
        "Ea_rv": Y[:,93],
        "Ea_lv": Y[:,94],
        "VAC_rv": Y[:,95],
        "VAC_lv": Y[:,96],
        "DeltaE": Y[:,97],
        "Patient_ID": Y[:,98]
    }
    output_df = pd.DataFrame(output)
    output_df.to_csv(f'{file_name_output}.csv', index=False)
    print(f"Computation completed for {file_name_specific_input} saved as {file_name_output}.csv", flush=True)





