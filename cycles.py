import numpy as np
from equations import *
from scipy.integrate import solve_ivp, simpson
from shapely.geometry import Polygon

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

    Era = activation_function(chamber_type="atrium", x=t_to_solve, tau1=tau1_ra, tau2=tau2_ra, m1=m1_ra, m2=m2_ra, k=k_ra, Emin=Emin_ra, periode=periode, delta=delta)
    Erv = activation_function(chamber_type="ventricle", x=t_to_solve, tau1=tau1_rv, tau2=tau2_rv, m1=m1_rv, m2=m2_rv, k=k_rv, Emin=Emin_rv)
    Ela = activation_function(chamber_type="atrium", x=t_to_solve, tau1=tau1_la, tau2=tau2_la, m1=m1_la, m2=m2_la, k=k_la, Emin=Emin_la, periode=periode, delta=delta)
    Elv = activation_function(chamber_type="ventricle", x=t_to_solve, tau1=tau1_lv, tau2=tau2_lv, m1=m1_lv, m2=m2_lv, k=k_lv, Emin=Emin_lv)

    Pra = Era*(Vra - V0ra)
    Prv = Erv*(Vrv - V0rv)
    Pla = Ela*(Vla - V0la)
    Plv = Elv*(Vlv - V0lv)

    #-----calculate output values-----#
    Pao = ((Qav + Pa_syst*conversion_factor/Rao + Pa_LAD*conversion_factor/Ra_LAD + Pa_LCX*conversion_factor/Ra_LCX + Pa_RCA*conversion_factor/Ra_RCA)/(1/Rao + 1/Ra_LAD + 1/Ra_LCX + 1/Ra_RCA))/conversion_factor
    



    STK_syst = max(Pa_syst)
    DTK_syst = min(Pa_syst)
    SV_lv = simpson(Qav, x=t_to_solve)
    VR_lv = simpson(Qmv, x=t_to_solve)


    STK_pulm = max(Pa_pulm)
    DTK_pulm = min(Pa_pulm)
    SV_rv = simpson(Qpv, x=t_to_solve)
    VR_rv = simpson(Qmv, x=t_to_solve)


    Qv_pulm = (Pv_pulm - Pla)/Rv_pulm
    Qv_syst = (Pv_syst - Pra)/Rv_syst
    Qv_LAD = (Pv_LAD - Pra)/(Rv_LAD/conversion_factor) #if Pv_LAD > Pra else 0.0
    Qv_LCX = (Pv_LCX - Pra)/(Rv_LCX/conversion_factor) #if Pv_LCX > Pra else 0.0
    Qv_RCA = (Pv_RCA - Pra)/(Rv_RCA/conversion_factor) #if Pv_RCA > Pra else 0.0
    Qv_syst_total = Qv_syst + Qv_LAD + Qv_LCX + Qv_RCA

    SV_ra = simpson(Qtv, x=t_to_solve)
    VR_ra = simpson(Qv_syst_total, x=t_to_solve)
    SV_la = simpson(Qmv, x=t_to_solve)
    VR_la = simpson(Qv_pulm, x=t_to_solve)

    ESV_lv = min(Vlv)
    ESV_rv = min(Vrv)



    #-----CBF model values-----#
    Rmicro_LAD = Et_Rmicro(t_to_solve, m1_lv, m2_lv, Rmicro_min_LAD, Rmicro_max_LAD, periode)
    Rmicro_LCX = Et_Rmicro(t_to_solve, m1_lv, m2_lv, Rmicro_min_LCX, Rmicro_max_LCX, periode)
    Rmicro_RCA = Et_Rmicro(t_to_solve, m1_rv, m2_rv, Rmicro_min_RCA, Rmicro_max_RCA, periode)

    Qa_LAD = (Pao - Pa_LAD)/(Ra_LAD/conversion_factor)
    Qa_LCX = (Pao - Pa_LCX)/(Ra_LCX/conversion_factor)
    Qa_RCA = (Pao - Pa_RCA)/(Ra_RCA/conversion_factor)

    Qmicro_LAD = (Pa_LAD- Pv_LAD)/(Rmicro_LAD/conversion_factor)
    Qmicro_LCX = (Pa_LCX - Pv_LCX)/(Rmicro_LCX/conversion_factor)
    Qmicro_RCA = (Pa_RCA - Pv_RCA)/(Rmicro_RCA/conversion_factor)

    Qa_CBF_total = simpson(Qa_LAD+Qa_LCX+Qa_RCA, x=t_to_solve)
    Qa_CBF_total = simpson(Qa_LAD, x=t_to_solve) + simpson(Qa_LCX, x=t_to_solve) + simpson(Qa_RCA, x=t_to_solve)
    Qmicro_CBF_total = simpson(Qmicro_LAD+Qmicro_LCX+Qmicro_RCA, x=t_to_solve)
    Qmicro_CBF_LAD = simpson(Qmicro_LAD, x=t_to_solve)
    Qmicro_CBF_LCX = simpson(Qmicro_LCX, x=t_to_solve)
    Qmicro_CBF_RCA = simpson(Qmicro_RCA, x=t_to_solve)
    Qv_CBF_total = simpson(Qv_LAD+Qv_LCX+Qv_RCA, x=t_to_solve)


    #-----energies-----#
    SW_rv_polygon = Polygon(zip(Vrv, Prv))
    SW_rv = SW_rv_polygon.area*(133.332*10**-6)
    SW_lv_polygon = Polygon(zip(Vlv, Plv))
    SW_lv = SW_lv_polygon.area*(133.332*10**-6)

    start_systole_lv_index = np.where(Qav>0.1)[0][0]
    end_systole_lv_index = np.where(Qav>0.1)[0][-1]
    start_systole_rv_index = np.where(Qpv>0.1)[0][0]
    end_systole_rv_index = np.where(Qpv>0.1)[0][-1]

    Pes_lv = Plv[end_systole_lv_index]
    Pes_rv = Prv[end_systole_rv_index]

    PE_lv = (0.5*Pes_lv*(ESV_lv-V0lv))*(133.332*10**-6)
    PE_rv = (0.5*Pes_rv*(ESV_rv-V0rv))*(133.332*10**-6)

    PVA_lv = SW_lv + PE_lv
    PVA_rv = SW_rv + PE_rv

    VO2_lv = 2.46*PVA_lv + 0.7
    VO2_rv = 2.46*PVA_rv + 0.2
    VO2_total = VO2_lv + VO2_rv

    Ein = Qmicro_CBF_total*20*0.2*0.8 
    Ein_LAD = Qmicro_CBF_LAD*20*0.2*0.8 
    Ein_LCX = Qmicro_CBF_LCX*20*0.2*0.8 
    Ein_RCA = Qmicro_CBF_RCA*20*0.2*0.8 

    LVET = t_to_solve[end_systole_lv_index] - t_to_solve[start_systole_lv_index]
    RVET = t_to_solve[end_systole_rv_index] - t_to_solve[start_systole_rv_index]
    Ees_lv = Pes_lv/(ESV_lv-V0lv)
    Ees_rv = Pes_rv/(ESV_rv-V0rv)

    Ea_lv = Pes_lv/SV_lv
    Ea_rv = Pes_rv/SV_rv
    VAC_lv = Ea_lv/Ees_lv
    VAC_rv = Ea_rv/Ees_rv

    PG_lv_mean = np.mean(Plv[start_systole_lv_index:end_systole_lv_index]-Pa_syst[start_systole_lv_index:end_systole_lv_index])


    #-----returning results-----#
    return {
        "t": t_to_solve,
        "chambers": {
            "Vra": Vra,
            "Vrv": Vrv,
            "Vla": Vla,
            "Vlv": Vlv,
            "Era": Era,
            "Erv": Erv,
            "Ela": Ela,
            "Elv": Elv,
            "Pra": Pra,
            "Prv": Prv,
            "Pla": Pla,
            "Plv": Plv
        },
        "valves": {
            "Qtv": Qtv,
            "Qpv": Qpv,
            "Qmv": Qmv,
            "Qav": Qav,
            "eps_tv": eps_tv,
            "eps_pv": eps_pv,
            "eps_mv": eps_mv,
            "eps_av": eps_av
        },
        "circulation": {
            "Pao": Pao,
            "Pa_syst": Pa_syst,
            "Pv_syst": Pv_syst,
            "Pa_pulm": Pa_pulm,
            "Pv_pulm": Pv_pulm
        },
        "coronary_circulation": {
            "Pa_LAD": Pa_LAD,
            "Pv_LAD": Pv_LAD,
            "Pa_LCX": Pa_LCX,
            "Pv_LCX": Pv_LCX,
            "Pa_RCA": Pa_RCA,
            "Pv_RCA": Pv_RCA,
            "Qa_LAD": Qa_LAD,
            "Qa_LCX": Qa_LCX,
            "Qa_RCA": Qa_RCA,
            "Qmicro_LAD": Qmicro_LAD,
            "Qmicro_LCX": Qmicro_LCX,
            "Qmicro_RCA": Qmicro_RCA,
            "Qv_LAD": Qv_LAD,
            "Qv_LCX": Qv_LCX,
            "Qv_RCA": Qv_RCA,
            "Qa_CBF_total": Qa_CBF_total,
            "Qmicro_CBF_total": Qmicro_CBF_total,
            "Qv_CBF_total": Qv_CBF_total
        },
        "output_values": {
            "STK_syst": STK_syst,
            "DTK_syst": DTK_syst,
            "MAP_syst": np.mean(Pa_syst),
            "SV_lv": SV_lv,
            "VR_lv": VR_lv,
            "EDV_lv": max(Vlv),
            "ESV_lv": min(Vlv),
            "Plvmax": max(Plv),
            "STK_pulm": STK_pulm,
            "DTK_pulm": DTK_pulm,
            "MAP_pulm": np.mean(Pa_pulm),
            "SV_rv": SV_rv,
            "VR_rv": VR_rv,
            "EDV_rv": max(Vrv),
            "Prvmax": max(Prv),
            "SV_la": SV_la,
            "VR_la": VR_la,
            "SV_ra": SV_ra,
            "VR_ra": VR_ra,
            "Ein": Ein,
            "Ein_LAD": Ein_LAD,
            "Ein_LCX": Ein_LCX,
            "Ein_RCA": Ein_RCA,
            "Eout": VO2_total,
            "Eout_lv": VO2_lv,
            "Eout_rv": VO2_rv,
            "SW_lv": SW_lv,
            "PE_lv": PE_lv,
            "SW_rv": SW_rv,
            "PE_rv": PE_rv,
            "maxVla": max(Vla),
            "minVla": min(Vla),
            "maxVra": max(Vra),
            "minVra": min(Vra),
            "LVET": LVET,
            "Pes_lv": Pes_lv,
            "Ees_lv": Ees_lv,
            "EF_lv": SV_lv/max(Vlv)*100,
            "EF_rv": SV_rv/max(Vrv)*100,
            "Plvmin": min(Plv),
            "EDP_lv": Plv[-1],
            "Pv_syst_max": max(Pv_syst),
            "Pv_syst_min": min(Pv_syst),
            "Pv_syst_mean": np.mean(Pv_syst),
            "Prvmin": min(Prv),
            "EDP_rv": Prv[-1],
            "Pv_pulm_max": max(Pv_pulm),
            "Pv_pulm_min": min(Pv_pulm),
            "Pv_pulm_mean": np.mean(Pv_pulm),
            "Plamax": max(Pla),
            "Plamin": min(Pla),
            "Pramax": max(Pra),
            "Pramin": min(Pra),
            "RVET": RVET,
            "Ees_rv": Ees_rv,
            "CO": SV_lv*HR,
            "Ea_lv": Ea_lv,
            "Ea_rv": Ea_rv,
            "VAC_lv": VAC_lv,
            "VAC_rv": VAC_rv,
            "Pes_rv": Pes_rv,
            "Qmicro_CBF_LAD": Qmicro_CBF_LAD,
            "Qmicro_CBF_LCX": Qmicro_CBF_LCX,
            "PG_lv_mean": PG_lv_mean
        }
    }


def multiple_cycles(parameters_input, number_of_cycles, check_stability=False, criterion=0.005, store_all=True):
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
    Rmicro_min_LAD = coronary_params["Rmicro_min_LAD"]*conversion_factor
    Rmicro_max_LCX = coronary_params["Rmicro_max_LCX"]*conversion_factor
    Rmicro_min_LCX = coronary_params["Rmicro_min_LCX"]*conversion_factor
    Rmicro_max_RCA = coronary_params["Rmicro_max_RCA"]*conversion_factor
    Rmicro_min_RCA = coronary_params["Rmicro_min_RCA"]*conversion_factor
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

    # compute
    for i in range(number_of_cycles):
        print(f"Current loop: {i+1}", end='\r')

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
        icVra = results["chambers"]["Vra"][-1]
        icVrv = results["chambers"]["Vrv"][-1]
        icVla = results["chambers"]["Vla"][-1]
        icVlv = results["chambers"]["Vlv"][-1]
        icQtv = results["valves"]["Qtv"][-1]
        icQpv = results["valves"]["Qpv"][-1]
        icQmv = results["valves"]["Qmv"][-1]
        icQav = results["valves"]["Qav"][-1]
        iceps_tv = results["valves"]["eps_tv"][-1]
        iceps_pv = results["valves"]["eps_pv"][-1]
        iceps_mv = results["valves"]["eps_mv"][-1]
        iceps_av = results["valves"]["eps_av"][-1]
        icPa_syst = results["circulation"]["Pa_syst"][-1]*conversion_factor
        icPv_syst = results["circulation"]["Pv_syst"][-1]*conversion_factor
        icPa_pulm = results["circulation"]["Pa_pulm"][-1]*conversion_factor
        icPv_pulm = results["circulation"]["Pv_pulm"][-1]*conversion_factor
        icPa_LAD = results["coronary_circulation"]["Pa_LAD"][-1]*conversion_factor
        icPv_LAD = results["coronary_circulation"]["Pv_LAD"][-1]*conversion_factor
        icPa_LCX = results["coronary_circulation"]["Pa_LCX"][-1]*conversion_factor
        icPv_LCX = results["coronary_circulation"]["Pv_LCX"][-1]*conversion_factor
        icPa_RCA = results["coronary_circulation"]["Pa_RCA"][-1]*conversion_factor
        icPv_RCA = results["coronary_circulation"]["Pv_RCA"][-1]*conversion_factor

        
        list_of_results.append(results)
    
        # check stability - 5 times
        if i > 1 and check_stability:
            if (
                abs(list_of_results[i]["output_values"]["STK_syst"] - list_of_results[i-1]["output_values"]["STK_syst"]) < criterion
                and abs(list_of_results[i]["output_values"]["DTK_syst"] - list_of_results[i-1]["output_values"]["DTK_syst"]) < criterion
                and abs(list_of_results[i]["output_values"]["SV_lv"] - list_of_results[i-1]["output_values"]["SV_lv"]) < criterion
                and abs(list_of_results[i]["output_values"]["SV_lv"] - list_of_results[i]["output_values"]["VR_lv"]) < criterion
                ):
                if stability_counter > 4:
                    print(f"Stability reached at cycle {i+1}.")
                    break
                else:
                    stability_counter += 1
            else:
                stability_counter = 0

    if stability_counter <= 4 and check_stability:
        print("Stability not reached.")
        print(parameters_input)

    if store_all:
        return list_of_results
    else:
        return results, i+1

if __name__ == "__main__":
    parameters_input = {
        "chambers": {
            "Emax_ra": ...,
            "Emin_ra": ...,
            "m1_ra": ...,
            "m2_ra": ...,
            "V0ra": ...,
            "Emax_rv": ...,
            "Emin_rv": ...,
            "m1_rv": ...,
            "m2_rv": ...,
            "V0rv": ...,
            "Emax_la": ...,
            "Emin_la": ...,
            "m1_la": ...,
            "m2_la": ...,
            "V0la": ...,
            "Emax_lv": ...,
            "Emin_lv": ...,
            "m1_lv": ...,
            "m2_lv": ...,
            "V0lv": ...
        },
        "valves": {
            "Amax_tv": ...,
            "Amin_tv": ...,
            "Kvo_tv": ...,
            "Kvc_tv": ...,
            "Amax_pv": ...,
            "Amin_pv": ...,
            "Kvo_pv": ...,
            "Kvc_pv": ...,
            "Amax_mv": ...,
            "Amin_mv": ...,
            "Kvo_mv": ...,
            "Kvc_mv": ...,
            "Amax_av": ...,
            "Amin_av": ...,
            "Kvo_av": ...,
            "Kvc_av": ...,
            
        },
        "circulation": {
            "HR": ...,
            "density": ...,
            "Rao": ...,
            "Ca_syst": ...,
            "Ra_syst": ...,
            "Cv_syst": ...,
            "Rv_syst": ...,
            "Rpa": ...,
            "Ca_pulm": ...,
            "Ra_pulm": ...,
            "Cv_pulm": ...,
            "Rv_pulm": ...,
        },
        "CBF_model": {
            "Ra_LAD": ...,
            "Ra_LCX": ...,
            "Ra_RCA": ...,
            "Rmicro_max_LAD": ...,
            "Rmicro_min_LAD": ...,
            "Rmicro_max_LCX": ...,
            "Rmicro_min_LCX": ...,
            "Rmicro_max_RCA": ...,
            "Rmicro_min_RCA": ...,
            "Rv_LAD": ...,
            "Rv_LCX": ...,
            "Rv_RCA": ...,
            "Ca_LAD": ...,
            "Ca_LCX": ...,
            "Ca_RCA": ...,
            "Cmicro_LAD": ...,
            "Cmicro_LCX": ...,
            "Cmicro_RCA": ...,
            "hypertrophy_factor":...
            
        },
        "initial_conditions": {
            "icVra": ...,
            "icVrv": ...,
            "icVla": ...,
            "icVlv": ...,
            "icQtv": ...,
            "icQpv": ...,
            "icQmv": ...,
            "icQav": ...,
            "iceps_tv": ...,
            "iceps_pv": ...,
            "iceps_mv": ...,
            "iceps_av": ...,
            "icPa_syst": ...,
            "icPv_syst": ...,
            "icPa_pulm": ...,
            "icPv_pulm": ...,
            "icPa_LAD": ...,
            "icPv_LAD": ...,
            "icPa_LCX": ...,
            "icPv_LCX": ...,
            "icPa_RCA": ...,
            "icPv_RCA": ...,
        }
    }

    list_of_results = multiple_cycles(parameters_input=parameters_input, number_of_cycles=10, check_stability=True, criterion=0.005, store_all=True)