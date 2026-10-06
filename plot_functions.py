import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def complete_plot_results(results, save=False, CBFplotall=False, output_path=None):
    # unpack results
    t = results["t"]
    Vra = results["chambers"]["Vra"]
    Vrv = results["chambers"]["Vrv"]
    Vla = results["chambers"]["Vla"]
    Vlv = results["chambers"]["Vlv"]
    Era = results["chambers"]["Era"]
    Erv = results["chambers"]["Erv"]
    Ela = results["chambers"]["Ela"]
    Elv = results["chambers"]["Elv"]
    Pra = results["chambers"]["Pra"]
    Prv = results["chambers"]["Prv"]
    Pla = results["chambers"]["Pla"]
    Plv = results["chambers"]["Plv"]
    Qtv = results["valves"]["Qtv"]
    Qpv = results["valves"]["Qpv"]
    Qmv = results["valves"]["Qmv"]
    Qav = results["valves"]["Qav"]
    eps_tv = results["valves"]["eps_tv"]
    eps_pv = results["valves"]["eps_pv"]
    eps_mv = results["valves"]["eps_mv"]
    eps_av = results["valves"]["eps_av"]
    Pa_syst = results["circulation"]["Pa_syst"]
    Pv_syst = results["circulation"]["Pv_syst"]
    Pa_pulm = results["circulation"]["Pa_pulm"]
    Pv_pulm = results["circulation"]["Pv_pulm"]
    Qa_LAD = results["coronary_circulation"]["Qa_LAD"]
    Qa_LCX = results["coronary_circulation"]["Qa_LCX"]
    Qa_RCA = results["coronary_circulation"]["Qa_RCA"]
    Qmicro_LAD = results["coronary_circulation"]["Qmicro_LAD"]
    Qmicro_LCX = results["coronary_circulation"]["Qmicro_LCX"]
    Qmicro_RCA = results["coronary_circulation"]["Qmicro_RCA"]
    Qv_LAD = results["coronary_circulation"]["Qv_LAD"]
    Qv_LCX = results["coronary_circulation"]["Qv_LCX"]
    Qv_RCA = results["coronary_circulation"]["Qv_RCA"]
    Pao = results["circulation"]["Pao"]

    centimeter = 1/2.54
    fig = plt.figure(figsize=(20*centimeter, 28*centimeter), layout="constrained")
    subfigs = fig.subfigures(nrows=3, ncols=1)

    #-----LEFT HEART-----#
    # plot results
    fig_leftheart = subfigs[0]
    fig_leftheart.suptitle("Left Heart", fontsize=12, fontweight='bold')
    axs_leftheart = fig_leftheart.subplots(nrows=2, ncols=3)

    # elastances
    axs_leftheart[0, 0].plot(t, Ela, label=r"$E_{la}$", linestyle="-")
    axs_leftheart[0, 0].plot(t, Elv, label=r"$E_{lv}$", linestyle="-")
    axs_leftheart[0, 0].set_xlabel("Time [s]")
    axs_leftheart[0, 0].set_ylabel("Elastance [mmHg/mL]")
    axs_leftheart[0, 0].set_title("Elastances")
    axs_leftheart[0, 0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # volumes
    axs_leftheart[0, 1].plot(t, Vla, label=r"$V_{la}$", linestyle="-")
    axs_leftheart[0, 1].plot(t, Vlv, label=r"$V_{lv}$", linestyle="-")
    axs_leftheart[0, 1].set_xlabel("Time [s]")
    axs_leftheart[0, 1].set_ylabel("Volume [mL]")
    axs_leftheart[0, 1].set_title("Volumes")
    axs_leftheart[0, 1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # flows
    axs_leftheart[0, 2].plot(t, Qmv, label=r"$Q_{mv}$", linestyle="-")
    axs_leftheart[0, 2].plot(t, Qav, label=r"$Q_{av}$", linestyle="-")
    axs_leftheart[0, 2].set_xlabel("Time [s]")
    axs_leftheart[0, 2].set_ylabel("Flow [mL/s]")
    axs_leftheart[0, 2].set_title("Flows")
    axs_leftheart[0, 2].legend(loc='upper right', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # pressures
    axs_leftheart[1, 0].plot(t, Pla, label=r"$P_{la}$", linestyle="-", color="r")
    axs_leftheart[1, 0].plot(t, Plv, label=r"$P_{lv}$", linestyle="-", color="g")
    axs_leftheart[1, 0].plot(t, Pa_syst, label=r"$P_{a}^{syst}$", linestyle="-", color="b")
    axs_leftheart[1, 0].plot(t, Pv_syst, label=r"$P_{v}^{syst}$", linestyle="-", color="brown")
    axs_leftheart[1, 0].plot(t, Pao, label=r"$P_{ao}$", linestyle="-", color="black")
    axs_leftheart[1, 0].set_xlabel("Time [s]")
    axs_leftheart[1, 0].set_ylabel("Pressure [mmHg]")
    axs_leftheart[1, 0].set_title("Pressures")
    axs_leftheart[1, 0].legend(loc='upper right', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # pV loop
    axs_leftheart[1, 1].plot(Vla, Pla, label=r"$LA$", linestyle="-")
    axs_leftheart[1, 1].plot(Vlv, Plv, label=r"$LV$", linestyle="-")
    axs_leftheart[1, 1].set_xlabel("Volume [mL]")
    axs_leftheart[1, 1].set_ylabel("Pressure [mmHg]")
    axs_leftheart[1, 1].set_title("pV Loops")
    axs_leftheart[1, 1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # valve states
    axs_leftheart[1, 2].plot(t, eps_mv, label=r"$\epsilon_{mv}$", linestyle="-")
    axs_leftheart[1, 2].plot(t, eps_av, label=r"$\epsilon_{av}$", linestyle="-")
    axs_leftheart[1, 2].set_xlabel("Time [s]")
    axs_leftheart[1, 2].set_ylabel("Valve State [-]")
    axs_leftheart[1, 2].set_title("Valve States")
    axs_leftheart[1, 2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)


    #-----Right HEART-----#
    # plot results
    fig_rightheart = subfigs[1]
    fig_rightheart.suptitle("Right Heart", fontsize=12, fontweight='bold')
    axs_rightheart = fig_rightheart.subplots(nrows=2, ncols=3)
    # elastances
    axs_rightheart[0, 0].plot(t, Era, label=r"$E_{ra}$", linestyle="-")
    axs_rightheart[0, 0].plot(t, Erv, label=r"$E_{rv}$", linestyle="-")
    axs_rightheart[0, 0].set_xlabel("Time [s]")
    axs_rightheart[0, 0].set_ylabel("Elastance [mmHg/mL]")
    axs_rightheart[0, 0].set_title("Elastances")
    axs_rightheart[0, 0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # volumes
    axs_rightheart[0, 1].plot(t, Vra, label=r"$V_{ra}$", linestyle="-")
    axs_rightheart[0, 1].plot(t, Vrv, label=r"$V_{rv}$", linestyle="-")
    axs_rightheart[0, 1].set_xlabel("Time [s]")
    axs_rightheart[0, 1].set_ylabel("Volume [mL]")
    axs_rightheart[0, 1].set_title("Volumes")
    axs_rightheart[0, 1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # flows
    axs_rightheart[0, 2].plot(t, Qtv, label=r"$Q_{tv}$", linestyle="-")
    axs_rightheart[0, 2].plot(t, Qpv, label=r"$Q_{pv}$", linestyle="-")
    axs_rightheart[0, 2].set_xlabel("Time [s]")
    axs_rightheart[0, 2].set_ylabel("Flow [mL/s]")
    axs_rightheart[0, 2].set_title("Flows")
    axs_rightheart[0, 2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # pressures
    axs_rightheart[1, 0].plot(t, Pra, label=r"$P_{ra}$", linestyle="-", color="r")
    axs_rightheart[1, 0].plot(t, Prv, label=r"$P_{rv}$", linestyle="-", color="g")
    axs_rightheart[1, 0].plot(t, Pa_pulm, label=r"$P_{a}^{pulm}$", linestyle="-", color="b")
    axs_rightheart[1, 0].plot(t, Pv_pulm, label=r"$P_{v}^{pulm}$", linestyle="-", color="brown")
    axs_rightheart[1, 0].set_xlabel("Time [s]")
    axs_rightheart[1, 0].set_ylabel("Pressure [mmHg]")
    axs_rightheart[1, 0].set_title("Pressures")
    axs_rightheart[1, 0].legend(loc='upper right', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # pV loop
    axs_rightheart[1, 1].plot(Vra, Pra, label=r"$RA$", linestyle="-")
    axs_rightheart[1, 1].plot(Vrv, Prv, label=r"$RV$", linestyle="-")
    axs_rightheart[1, 1].set_xlabel("Volume [mL]")
    axs_rightheart[1, 1].set_ylabel("Pressure [mmHg]")
    axs_rightheart[1, 1].set_title("pV Loop")
    axs_rightheart[1, 1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    # valve states
    axs_rightheart[1, 2].plot(t, eps_tv, label=r"$\epsilon_{tv}$", linestyle="-")
    axs_rightheart[1, 2].plot(t, eps_pv, label=r"$\epsilon_{pv}$", linestyle="-")
    axs_rightheart[1, 2].set_xlabel("Time [s]")
    axs_rightheart[1, 2].set_ylabel("Valve State [-]")
    axs_rightheart[1, 2].set_title("Valve States")
    axs_rightheart[1, 2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    #-----CBF model-----#
    # plot results
    fig_CBF = subfigs[2]
    fig_CBF.suptitle("Coronary Blood Flow", fontsize=12, fontweight='bold')
    if CBFplotall:
        axs_CBF = fig_CBF.subplots(nrows=3, ncols=3)
        # flows
        axs_CBF[0,0].plot(t, Qa_LAD, label=r"$Q_{a}^{LAD}$", linestyle="-")
        axs_CBF[0,0].set_xlabel("Time [s]")
        axs_CBF[0,0].set_ylabel("Flow [mL/s]")
        axs_CBF[0,0].set_title(r"$Q_{a}^{LAD}$")
        axs_CBF[0,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[0,1].plot(t, Qa_LCX, label=r"$Q_{a}^{LCX}$", linestyle="-")
        axs_CBF[0,1].set_xlabel("Time [s]")
        axs_CBF[0,1].set_ylabel("Flow [mL/s]")
        axs_CBF[0,1].set_title(r"$Q_{a}^{LCX}$")
        axs_CBF[0,1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[0,2].plot(t, Qa_RCA, label=r"$Q_{a}^{RCA}$", linestyle="-")
        axs_CBF[0,2].set_xlabel("Time [s]")
        axs_CBF[0,2].set_ylabel("Flow [mL/s]")
        axs_CBF[0,2].set_title(r"$Q_{a}^{RCA}$")
        axs_CBF[0,2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[1,0].plot(t, Qmicro_LAD, label=r"$Q_{micro}^{LAD}$", linestyle="-")
        axs_CBF[1,0].set_xlabel("Time [s]")
        axs_CBF[1,0].set_ylabel("Flow [mL/s]")
        axs_CBF[1,0].set_title(r"$Q_{micro}^{LAD}$")
        axs_CBF[1,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[1,1].plot(t, Qmicro_LCX, label=r"$Q_{micro}^{LCX}$", linestyle="-")
        axs_CBF[1,1].set_xlabel("Time [s]")
        axs_CBF[1,1].set_ylabel("Flow [mL/s]")
        axs_CBF[1,1].set_title(r"$Q_{micro}^{LCX}$")
        axs_CBF[1,1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[1,2].plot(t, Qmicro_RCA, label=r"$Q_{micro}^{RCA}$", linestyle="-")
        axs_CBF[1,2].set_xlabel("Time [s]")
        axs_CBF[1,2].set_ylabel("Flow [mL/s]")
        axs_CBF[1,2].set_title(r"$Q_{micro}^{RCA}$")
        axs_CBF[1,2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[2,0].plot(t, Qv_LAD, label=r"$Q_{v}^{LAD}$", linestyle="-")
        axs_CBF[2,0].set_xlabel("Time [s]")
        axs_CBF[2,0].set_ylabel("Flow [mL/s]")
        axs_CBF[2,0].set_title(r"$Q_{v}^{LAD}$")
        axs_CBF[2,0].legend(loc='upper right', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[2,1].plot(t, Qv_LCX, label=r"$Q_{v}^{LCX}$", linestyle="-")
        axs_CBF[2,1].set_xlabel("Time [s]")
        axs_CBF[2,1].set_ylabel("Flow [mL/s]")
        axs_CBF[2,1].set_title(r"$Q_{v}^{LCX}$")
        axs_CBF[2,1].legend(loc='upper right', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[2,2].plot(t, Qv_RCA, label=r"$Q_{v}^{RCA}$", linestyle="-")
        axs_CBF[2,2].set_xlabel("Time [s]")
        axs_CBF[2,2].set_ylabel("Flow [mL/s]")
        axs_CBF[2,2].set_title(r"$Q_{v}^{RCA}$")
        axs_CBF[2,2].legend(loc='upper right', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)
    else:
        axs_CBF = fig_CBF.subplots(nrows=1, ncols=3)
        axs_CBF[0].plot(t, Qmicro_LAD, label=r"$Q_{micro}^{LAD}$", linestyle="-")
        axs_CBF[0].set_xlabel("Time [s]")
        axs_CBF[0].set_ylabel("Flow [mL/s]")
        axs_CBF[0].set_title(r"$Q_{micro}^{LAD}$")
        axs_CBF[0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[1].plot(t, Qmicro_LCX, label=r"$Q_{micro}^{LCX}$", linestyle="-")
        axs_CBF[1].set_xlabel("Time [s]")
        axs_CBF[1].set_ylabel("Flow [mL/s]")
        axs_CBF[1].set_title(r"$Q_{micro}^{LCX}$")
        axs_CBF[1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

        axs_CBF[2].plot(t, Qmicro_RCA, label=r"$Q_{micro}^{RCA}$", linestyle="-")
        axs_CBF[2].set_xlabel("Time [s]")
        axs_CBF[2].set_ylabel("Flow [mL/s]")
        axs_CBF[2].set_title(r"$Q_{micro}^{RCA}$")
        axs_CBF[2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False)

    

    # print("Qa_CBF_total: ", results["coronary_circulation"]["Qa_CBF_total"])
    # print("Qmicro_CBF_total: ", results["coronary_circulation"]["Qmicro_CBF_total"])
    # print("Qv_CBF_total: ", results["coronary_circulation"]["Qv_CBF_total"])

    if save:
        plt.savefig(output_path)
    
    plt.show()


def SVA_plot_from_csv(path_noAS, path_AS, save_fig=False, complete_plot=False, show_all=False):
    # loading data
    df_noAS = pd.read_csv(path_noAS)

    STK_syst_list_noAS = df_noAS['STK_syst_list'].tolist()
    DTK_syst_list_noAS = df_noAS['DTK_syst_list'].tolist()
    MAP_syst_list_noAS = df_noAS['MAP_syst_list'].tolist()
    Plvmax_list_noAS = df_noAS['Plvmax_list'].tolist()
    SV_lv_list_noAS = df_noAS['SV_lv_list'].tolist()
    EDV_lv_list_noAS = df_noAS['EDV_lv_list'].tolist()
    STK_pulm_list_noAS = df_noAS['STK_pulm_list'].tolist()
    DTK_pulm_list_noAS = df_noAS['DTK_pulm_list'].tolist()
    MAP_pulm_list_noAS = df_noAS['MAP_pulm_list'].tolist()
    Prvmax_list_noAS = df_noAS['Prvmax_list'].tolist()
    SV_rv_list_noAS = df_noAS['SV_rv_list'].tolist()
    EDV_rv_list_noAS = df_noAS['EDV_rv_list'].tolist()
    maxVla_list_noAS = df_noAS['maxVla_list'].tolist()
    minVla_list_noAS = df_noAS['minVla_list'].tolist()
    maxVra_list_noAS = df_noAS['maxVra_list'].tolist()
    minVra_list_noAS = df_noAS['minVra_list'].tolist()
    Eout_total_list_noAS = df_noAS['Eout_total_list'].tolist()
    Eout_lv_list_noAS = df_noAS['Eout_lv_list'].tolist()
    Eout_rv_list_noAS = df_noAS['Eout_rv_list'].tolist()
    SW_lv_list_noAS = df_noAS['SW_lv_list'].tolist()
    PE_lv_list_noAS = df_noAS['PE_lv_list'].tolist()
    SW_rv_list_noAS = df_noAS['SW_rv_list'].tolist()
    PE_rv_list_noAS = df_noAS['PE_rv_list'].tolist()
    Ein_total_list_noAS = df_noAS['Ein_total_list'].tolist()
    Ein_LAD_list_noAS = df_noAS['Ein_LAD_list'].tolist()
    Ein_LCX_list_noAS = df_noAS['Ein_LCX_list'].tolist()
    Ein_RCA_list_noAS = df_noAS['Ein_RCA_list'].tolist()
    Ees_lv_list_noAS = df_noAS['Ees_lv_list'].tolist()
    Pes_lv_list_noAS = df_noAS['Pes_lv_list'].tolist()
    LVET_list_noAS = df_noAS['LVET_list'].tolist()
    EF_lv_list_noAS = df_noAS["EF_lv_list"].tolist()
    EF_rv_list_noAS = df_noAS["EF_rv_list"].tolist()
    Plvmin_list_noAS = df_noAS["Plvmin_list"].tolist()
    EDP_lv_list_noAS = df_noAS["EDP_lv_list"].tolist()
    Pv_syst_max_list_noAS = df_noAS["Pv_syst_max_list"].tolist()
    Pv_syst_min_list_noAS = df_noAS["Pv_syst_min_list"].tolist()
    Pv_syst_mean_list_noAS = df_noAS["Pv_syst_mean_list"].tolist()
    Prvmin_list_noAS = df_noAS["Prvmin_list"].tolist()
    EDP_rv_list_noAS = df_noAS["EDP_rv_list"].tolist()
    Pv_pulm_max_list_noAS = df_noAS["Pv_pulm_max_list"].tolist()
    Pv_pulm_min_list_noAS = df_noAS["Pv_pulm_min_list"].tolist()
    Pv_pulm_mean_list_noAS = df_noAS["Pv_pulm_mean_list"].tolist()
    Plamax_list_noAS = df_noAS["Plamax_list"].tolist()
    Pramax_list_noAS = df_noAS["Pramax_list"].tolist()
    Plamin_list_noAS = df_noAS["Plamin_list"].tolist()
    Pramin_list_noAS = df_noAS["Pramin_list"].tolist()
    RVET_list_noAS = df_noAS["RVET_list"].tolist()
    Ees_rv_list_noAS = df_noAS["Ees_rv_list"].tolist()
    CO_list_noAS = df_noAS["CO_list"].tolist()
    Ea_lv_list_noAS = df_noAS["Ea_lv_list"].tolist()
    Ea_rv_list_noAS = df_noAS["Ea_lv_list"].tolist()
    VAC_lv_list_noAS = df_noAS["VAC_lv_list"].tolist()
    VAC_rv_list_noAS = df_noAS["VAC_rv_list"].tolist()

    df_AS = pd.read_csv(path_AS)

    STK_syst_list_AS = df_AS['STK_syst_list'].tolist()
    DTK_syst_list_AS = df_AS['DTK_syst_list'].tolist()
    MAP_syst_list_AS = df_AS['MAP_syst_list'].tolist()
    Plvmax_list_AS = df_AS['Plvmax_list'].tolist()
    SV_lv_list_AS = df_AS['SV_lv_list'].tolist()
    EDV_lv_list_AS = df_AS['EDV_lv_list'].tolist()
    STK_pulm_list_AS = df_AS['STK_pulm_list'].tolist()
    DTK_pulm_list_AS = df_AS['DTK_pulm_list'].tolist()
    MAP_pulm_list_AS = df_AS['MAP_pulm_list'].tolist()
    Prvmax_list_AS = df_AS['Prvmax_list'].tolist()
    SV_rv_list_AS = df_AS['SV_rv_list'].tolist()
    EDV_rv_list_AS = df_AS['EDV_rv_list'].tolist()
    maxVla_list_AS = df_AS['maxVla_list'].tolist()
    minVla_list_AS = df_AS['minVla_list'].tolist()
    maxVra_list_AS = df_AS['maxVra_list'].tolist()
    minVra_list_AS = df_AS['minVra_list'].tolist()
    Eout_total_list_AS = df_AS['Eout_total_list'].tolist()
    Eout_lv_list_AS = df_AS['Eout_lv_list'].tolist()
    Eout_rv_list_AS = df_AS['Eout_rv_list'].tolist()
    SW_lv_list_AS = df_AS['SW_lv_list'].tolist()
    PE_lv_list_AS = df_AS['PE_lv_list'].tolist()
    SW_rv_list_AS = df_AS['SW_rv_list'].tolist()
    PE_rv_list_AS = df_AS['PE_rv_list'].tolist()
    Ein_total_list_AS = df_AS['Ein_total_list'].tolist()
    Ein_LAD_list_AS = df_AS['Ein_LAD_list'].tolist()
    Ein_LCX_list_AS = df_AS['Ein_LCX_list'].tolist()
    Ein_RCA_list_AS = df_AS['Ein_RCA_list'].tolist()
    Ees_lv_list_AS = df_AS['Ees_lv_list'].tolist()
    Pes_lv_list_AS = df_AS['Pes_lv_list'].tolist()
    LVET_list_AS = df_AS['LVET_list'].tolist()
    EF_lv_list_AS = df_AS["EF_lv_list"].tolist()
    EF_rv_list_AS = df_AS["EF_rv_list"].tolist()
    Plvmin_list_AS = df_AS["Plvmin_list"].tolist()
    EDP_lv_list_AS = df_AS["EDP_lv_list"].tolist()
    Pv_syst_max_list_AS = df_AS["Pv_syst_max_list"].tolist()
    Pv_syst_min_list_AS = df_AS["Pv_syst_min_list"].tolist()
    Pv_syst_mean_list_AS = df_AS["Pv_syst_mean_list"].tolist()
    Prvmin_list_AS = df_AS["Prvmin_list"].tolist()
    EDP_rv_list_AS = df_AS["EDP_rv_list"].tolist()
    Pv_pulm_max_list_AS = df_AS["Pv_pulm_max_list"].tolist()
    Pv_pulm_min_list_AS = df_AS["Pv_pulm_min_list"].tolist()
    Pv_pulm_mean_list_AS = df_AS["Pv_pulm_mean_list"].tolist()
    Plamax_list_AS = df_AS["Plamax_list"].tolist()
    Pramax_list_AS = df_AS["Pramax_list"].tolist()
    Plamin_list_AS = df_AS["Plamin_list"].tolist()
    Pramin_list_AS = df_AS["Pramin_list"].tolist()
    RVET_list_AS = df_AS["RVET_list"].tolist()
    Ees_rv_list_AS = df_AS["Ees_rv_list"].tolist()
    CO_list_AS = df_AS["CO_list"].tolist()
    Ea_lv_list_AS = df_AS["Ea_lv_list"].tolist()
    Ea_rv_list_AS = df_AS["Ea_lv_list"].tolist()
    VAC_lv_list_AS = df_AS["VAC_lv_list"].tolist()
    VAC_rv_list_AS = df_AS["VAC_rv_list"].tolist()

    var = df_noAS["var"].tolist()
    var_name = df_noAS["var_name"][0]
    var_unit = df_noAS["var_unit"][0]

    # plotting
    
        
    centimeter = 1/2.54

    fig, axs = plt.subplots(2, 3, figsize=(20*centimeter, 11*centimeter), layout='constrained', sharex=True)

    # axs[0,0].plot(var, y, label=r"$SBP^{syst}$", color="blue", linestyle='-')
    # axs[0,0].plot(var, y, label=r"$DBP^{syst}$", color="brown", linestyle='-')
    axs[0,0].plot(var, MAP_syst_list_noAS, label=r"$MAP_{syst}$", color="red", linestyle='-')
    axs[0,0].plot(var, Plvmax_list_noAS, label=r"$P_{lv}^{max}$", color="green", linestyle='-')
    ax_twin00 = axs[0,0].twinx()
    ax_twin00.plot(var, EDV_lv_list_noAS, label=r"$EDV_{lv}$", color="blue", linestyle='--')
    ax_twin00.plot(var, SV_lv_list_noAS, label=r"$SV_{lv}$", color="purple", linestyle='--')
    axs[0,0].set_ylabel("Normal Aortic Valve" + "\n" + "\n" + r"Pressure [mmHg]", 
                fontsize=8, 
                multialignment='center') 
    axs[0,0].set_title("Systemic Circulation", fontsize=10)
    lines1, labels1 = axs[0,0].get_legend_handles_labels()
    lines2, labels2 = ax_twin00.get_legend_handles_labels()
    axs[0,0].legend(lines1 + lines2, labels1 + labels2, loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

    # axs[0,1].plot(var, y, label=r"$SBP^{pulm}$", color="blue", linestyle='-')
    # axs[0,1].plot(var, y, label=r"$DBP^{pulm}$", color="brown", linestyle='-')
    axs[0,1].plot(var, MAP_pulm_list_noAS, label=r"$MAP_{pulm}$", color="red", linestyle='-')
    axs[0,1].plot(var, Prvmax_list_noAS, label=r"$P_{rv}^{max}$", color="green", linestyle='-')
    ax_twin01 = axs[0,1].twinx()
    ax_twin01.plot(var, EDV_rv_list_noAS, label=r"$EDV_{rv}$", color="blue", linestyle='--')
    ax_twin01.plot(var, SV_rv_list_noAS, label=r"$SV_{rv}$", color="purple", linestyle='--')
    #axs[0,1].set_ylabel("Pressure [mmHg]")
    ax_twin01.set_ylabel(f"Volume [ml]")
    axs[0,1].set_title("Pulmonary Circulation")
    lines1, labels1 = axs[0,1].get_legend_handles_labels()
    lines2, labels2 = ax_twin01.get_legend_handles_labels()
    axs[0,1].legend(lines1 + lines2, labels1 + labels2, loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'


    axs[0,2].plot(var, Ein_total_list_noAS, label=r"$E_{in}$", color="blue", linestyle='-')
    axs[0,2].plot(var, Eout_total_list_noAS, label=r"$E_{out}$", color="y", linestyle='-')
    axs[0,2].plot(var, np.array(Ein_total_list_noAS)-np.array(Eout_total_list_noAS), label=r"$\Delta E$", color="red")
    axs[0,2].set_ylabel("Energy [J/beat]")
    axs[0,2].set_title("Energies")
    axs[0,2].legend(loc='upper left', ncol=3, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'


    # axs[1,0].plot(var, y, label=r"$SBP^{syst}$", color="blue", linestyle='-')
    # axs[1,0].plot(var, y, label=r"$DBP^{syst}$", color="brown", linestyle='-')
    axs[1,0].plot(var, MAP_syst_list_AS, label=r"$MAP_{syst}$", color="red", linestyle='-')
    axs[1,0].plot(var, Plvmax_list_AS, label=r"$P_{lv}^{max}$", color="green", linestyle='-')
    ax_twin10 = axs[1,0].twinx()
    ax_twin10.plot(var, EDV_lv_list_AS, label=r"$EDV_{lv}$", color="blue", linestyle='--')
    ax_twin10.plot(var, SV_lv_list_AS, label=r"$SV_{lv}$", color="purple", linestyle='--')
    axs[1,0].set_xlabel(f"{var_name} {var_unit}")
    axs[1,0].set_ylabel("Severe Aortic Stenosis" + "\n" + "\n" + r"Pressure [mmHg]", 
                fontsize=8, 
                multialignment='center')
    lines1, labels1 = axs[1,0].get_legend_handles_labels()
    lines2, labels2 = ax_twin10.get_legend_handles_labels()
    axs[1,0].legend(lines1 + lines2, labels1 + labels2, loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

    # axs[1,1].plot(var, y, label=r"$SBP^{pulm}$", color="blue", linestyle='-')
    # axs[1,1].plot(var, y, label=r"$DBP^{pulm}$", color="brown", linestyle='-')
    axs[1,1].plot(var, MAP_pulm_list_AS, label=r"$MAP_{pulm}$", color="red", linestyle='-')
    axs[1,1].plot(var, Prvmax_list_AS, label=r"$P_{rv}^{max}$", color="green", linestyle='-')
    ax_twin11 = axs[1,1].twinx()
    ax_twin11.plot(var, EDV_rv_list_AS, label=r"$EDV_{rv}$", color="blue", linestyle='--')
    ax_twin11.plot(var, SV_rv_list_AS, label=r"$SV_{rv}$", color="purple", linestyle='--')
    axs[1,1].set_xlabel(f"{var_name} {var_unit}")
    #axs[1,1].set_ylabel("Pressure [mmHg]")
    ax_twin11.set_ylabel(f"Volume [ml]")
    lines1, labels1 = axs[1,1].get_legend_handles_labels()
    lines2, labels2 = ax_twin11.get_legend_handles_labels()
    axs[1,1].legend(lines1 + lines2, labels1 + labels2, loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

    axs[1,2].plot(var, Ein_total_list_AS, label=r"$E_{in}$", color="blue", linestyle='-')
    axs[1,2].plot(var, Eout_total_list_AS, label=r"$E_{out}$", color="y", linestyle='-')
    axs[1,2].plot(var, np.array(Ein_total_list_AS)-np.array(Eout_total_list_AS), label=r"$\Delta E$", color="red")
    axs[1,2].set_xlabel(f"{var_name} {var_unit}")
    axs[1,2].set_ylabel("Energy [J]")
    axs[1,2].legend(loc='upper left', ncol=3, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'


    if save_fig:
        plt.savefig(f"SVA_figures/final_{var_name}.png", dpi=400) 
    plt.show()


    if complete_plot:
        fig, axs = plt.subplots(2, 4, figsize=(16, 7), layout='constrained', sharex=True)

        axs[0,0].plot(var, STK_syst_list_noAS, label=r"$SBP^{syst}$", color="blue", linestyle='-')
        axs[0,0].plot(var, DTK_syst_list_noAS, label=r"$DBP^{syst}$", color="brown", linestyle='-')
        axs[0,0].plot(var, MAP_syst_list_noAS, label=r"$MAP^{syst}$", color="red", linestyle='-')
        axs[0,0].plot(var, Pv_syst_mean_list_noAS, label=r"$CVP$", color="pink", linestyle='-')
        axs[0,0].plot(var, Plvmax_list_noAS, label=r"$P_{lv}^{max}$", color="green", linestyle='-')
        axs[0,0].set_ylabel("Normal Aortic Valve" + "\n" + "\n" + r"Pressure [mmHg]", 
                    fontsize=12, 
                    multialignment='center') 
        axs[0,0].set_title("Pressures - Systemic Circulation")
        axs[0,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,1].plot(var, STK_pulm_list_noAS, label=r"$SBP^{pulm}$", color="blue", linestyle='-')
        axs[0,1].plot(var, DTK_pulm_list_noAS, label=r"$DBP^{pulm}$", color="brown", linestyle='-')
        axs[0,1].plot(var, MAP_pulm_list_noAS, label=r"$MAP^{pulm}$", color="red", linestyle='-')
        axs[0,1].plot(var, Pv_pulm_mean_list_noAS, label=r"$PVP$", color="pink", linestyle='-')
        axs[0,1].plot(var, Prvmax_list_noAS, label=r"$P_{rv,max}$", color="green", linestyle='-')
        axs[0,1].set_ylabel("Pressure [mmHg]")
        axs[0,1].set_title("Pressures - Pulmonary Circulation")
        axs[0,1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,2].plot(var, EDV_lv_list_noAS, label=r"$EDV_{lv}$", color="green", linestyle='-')
        axs[0,2].plot(var, EDV_rv_list_noAS, label=r"$EDV_{rv}$", color="brown", linestyle='-')
        axs[0,2].plot(var, SV_lv_list_noAS, label=r"$SV_{lv}$", color="blue", linestyle='-')
        ax_twin1 = axs[0,2].twinx()
        ax_twin1.plot(var, EF_lv_list_noAS, label=r"$EF_{lv}$", color="green", linestyle="--")
        ax_twin1.plot(var, EF_rv_list_noAS, label=r"$EF_{rv}$", color="brown", linestyle="--")
        axs[0,2].set_ylabel("Volume [ml]")
        ax_twin1.set_ylabel(r"EF [\%]")
        axs[0,2].set_title("Volumes")
        lines1, labels1 = axs[0,2].get_legend_handles_labels()
        lines2, labels2 = ax_twin1.get_legend_handles_labels()
        axs[0,2].legend(lines1 + lines2, labels1 + labels2, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,3].plot(var, Ein_total_list_noAS, label=r"$E_{in}$", color="green", linestyle='-')
        axs[0,3].plot(var, Eout_total_list_noAS, label=r"$E_{out}$", color="red", linestyle='-')
        axs[0,3].plot(var, np.array(Ein_total_list_noAS)-np.array(Eout_total_list_noAS), label=r"$\Delta E$", color="blue")
        # ax_twin03 = axs[0,3].twinx()
        # ax_twin03.plot(var, VAC_lv_list_noAS, label=r"$VAC$", color="black", linestyle="--")
        axs[0,3].set_ylabel("Energy [J/beat]")
        # ax_twin03.set_ylabel("VAC [-]")
        axs[0,3].set_title("Energies")
        # lines03, labels03 = axs[0,3].get_legend_handles_labels()
        # lines03t, labels03t = ax_twin03.get_legend_handles_labels()
        # axs[0,3].legend(lines03 + lines03t, labels03 + labels03t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'
        axs[0,3].legend(loc='upper left', ncol=3, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        
        axs[1,0].plot(var, STK_syst_list_AS, label=r"$SBP^{syst}$", color="blue", linestyle='-')
        axs[1,0].plot(var, DTK_syst_list_AS, label=r"$DBP^{syst}$", color="brown", linestyle='-')
        axs[1,0].plot(var, MAP_syst_list_AS, label=r"$MAP^{syst}$", color="red", linestyle='-')
        axs[1,0].plot(var, Pv_syst_mean_list_AS, label=r"$CVP$", color="pink", linestyle='-')
        axs[1,0].plot(var, Plvmax_list_AS, label=r"$P_{lv,max}$", color="green", linestyle='-')
        axs[1,0].set_xlabel(f"{var_name} {var_unit}")
        axs[1,0].set_ylabel("Severe Aortic Stenosis" + "\n" + "\n" + r"Pressure [mmHg]", 
                    fontsize=12, 
                    multialignment='center')
        axs[1,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,1].plot(var, STK_pulm_list_AS, label=r"$SBP^{pulm}$", color="blue", linestyle='-')
        axs[1,1].plot(var, DTK_pulm_list_AS, label=r"$DBP^{pulm}$", color="brown", linestyle='-')
        axs[1,1].plot(var, MAP_pulm_list_AS, label=r"$MAP^{pulm}$", color="red", linestyle='-')
        axs[1,1].plot(var, Pv_pulm_mean_list_AS, label=r"$PVP$", color="pink", linestyle='-')
        axs[1,1].plot(var, Prvmax_list_AS, label=r"$P_{rv,max}$", color="green", linestyle='-')
        axs[1,1].set_xlabel(f"{var_name} {var_unit}")
        axs[1,1].set_ylabel("Pressure [mmHg]")
        axs[1,1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,2].plot(var, EDV_lv_list_AS, label=r"$EDV_{lv}$", color="green", linestyle='-')
        axs[1,2].plot(var, EDV_rv_list_AS, label=r"$EDV_{rv}$", color="brown", linestyle='-')
        axs[1,2].plot(var, SV_lv_list_AS, label=r"$SV_{lv}$", color="blue", linestyle='-')
        ax_twin2 = axs[1,2].twinx()
        ax_twin2.plot(var, EF_lv_list_AS, label=r"$EF_{lv}$", color="green", linestyle="--")
        ax_twin2.plot(var, EF_rv_list_AS, label=r"$EF_{rv}$", color="brown", linestyle="--")
        axs[1,2].set_xlabel(f"{var_name} {var_unit}")
        axs[1,2].set_ylabel("Volume [ml]")
        ax_twin2.set_ylabel("VAC [-]")
        ax_twin2.set_ylabel(r"EF [\%]")
        lines3, labels3 = axs[1,2].get_legend_handles_labels()
        lines4, labels4 = ax_twin2.get_legend_handles_labels()
        axs[1,2].legend(lines3+lines4, labels3+labels4, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'


        axs[1,3].plot(var, Ein_total_list_AS, label=r"$E_{in}$", color="green", linestyle='-')
        axs[1,3].plot(var, Eout_total_list_AS, label=r"$E_{out}$", color="red", linestyle='-')
        axs[1,3].plot(var, np.array(Ein_total_list_AS)-np.array(Eout_total_list_AS), label=r"$\Delta E$", color="blue")
        # ax_twin13 = axs[1,3].twinx()
        # ax_twin13.plot(var, VAC_lv_list_AS, label=r"$VAC$", color="black", linestyle="--")
        axs[1,3].set_xlabel(f"{var_name} {var_unit}")
        axs[1,3].set_ylabel("Energy [J]")
        # ax_twin13.set_ylabel("VAC [-]")
        # lines13, labels13 = axs[1,3].get_legend_handles_labels()
        # lines13t, labels13t = ax_twin13.get_legend_handles_labels()
        # axs[1,3].legend(lines13 + lines13t, labels13 + labels13t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'
        axs[1,3].legend(loc='upper left', ncol=3, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        
        if save_fig:
            plt.savefig(f"SVA_figures/SVA_{var_name}_1.png") 
        if show_all:
            plt.show()

        fig, axs = plt.subplots(2, 4, figsize=(15, 7), layout='constrained', sharex=True)

        axs[0,0].plot(var, Plvmax_list_noAS, label=r"$P_{lv}^{max}$", color="blue", linestyle='-')
        axs[0,0].plot(var, Plvmin_list_noAS, label=r"$P_{lv}^{min}$", color="brown", linestyle='-')
        axs[0,0].plot(var, Prvmax_list_noAS, label=r"$P_{rv}^{max}$", color="red", linestyle='-')
        axs[0,0].plot(var, Prvmin_list_noAS, label=r"$P_{rv}^{min}$", color="green", linestyle='-')
        axs[0,0].plot(var, EDP_lv_list_noAS, label=r"$EDP_{lv}$", color="black", linestyle='-')
        axs[0,0].plot(var, EDP_rv_list_noAS, label=r"$EDP_{rv}$", color="brown", linestyle='-')
        axs[0,0].set_ylabel("Normal Aortic Valve" + "\n" + "\n" + r"Pressure [mmHg]", 
                    fontsize=12, 
                    multialignment='center') #r"$A_{max,av}^{eff} = 4.0$" +  r"$ cm^2$"
        axs[0,0].set_title("Ventricles - Pressures")
        axs[0,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,1].plot(var, EDV_lv_list_noAS, label=r"$EDV_{lv}$", color="green", linestyle='-')
        axs[0,1].plot(var, EDV_rv_list_noAS, label=r"$EDV_{rv}$", color="brown", linestyle='-')
        axs[0,1].plot(var, SV_lv_list_noAS, label=r"$SV_{lv}$", color="blue", linestyle='-')
        ax_twin01 = axs[0,1].twinx()
        ax_twin01.plot(var, EF_lv_list_noAS, label=r"$EF_{lv}$", color="green", linestyle="--")
        ax_twin01.plot(var, EF_rv_list_noAS, label=r"$EF_{rv}$", color="brown", linestyle="--")
        axs[0,1].set_ylabel("Volume [ml]")
        ax_twin01.set_ylabel(r"EF [\%]")
        axs[0,1].set_title("Ventricles - Volumes")
        lines01, labels01 = axs[0,1].get_legend_handles_labels()
        lines01t, labels01t = ax_twin01.get_legend_handles_labels()
        axs[0,1].legend(lines01+lines01t, labels01+labels01t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,2].plot(var, Plamax_list_noAS, label=r"$P_{la}^{max}$", color="green", linestyle='-')
        axs[0,2].plot(var, Plamin_list_noAS, label=r"$P_{la}^{min}$", color="red", linestyle='-')
        axs[0,2].plot(var, Pramax_list_noAS, label=r"$P_{ra}^{max}$", color="blue", linestyle='-')
        axs[0,2].plot(var, Pramin_list_noAS, label=r"$P_{ra}^{min}$", color="brown", linestyle='-')
        axs[0,2].set_ylabel("Pressure [mmHg]")
        axs[0,2].set_title("Atria - Pressures")
        axs[0,2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,3].plot(var, maxVla_list_noAS, label=r"$V_{la}^{max}$", color="green", linestyle='-')
        axs[0,3].plot(var, minVla_list_noAS, label=r"$V_{la}^{min}$", color="red", linestyle='-')
        axs[0,3].plot(var, maxVra_list_noAS, label=r"$V_{ra}^{max}$", color="blue", linestyle='-')
        axs[0,3].plot(var, minVra_list_noAS, label=r"$V_{ra}^{min}$", color="brown", linestyle='-')
        axs[0,3].set_ylabel("Volume [ml]")
        axs[0,3].set_title("Atria - Volumes")
        axs[0,3].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,0].plot(var, Plvmax_list_AS, label=r"$P_{lv}^{max}$", color="blue", linestyle='-')
        axs[1,0].plot(var, Plvmin_list_AS, label=r"$P_{lv}^{min}$", color="brown", linestyle='-')
        axs[1,0].plot(var, Prvmax_list_AS, label=r"$P_{rv}^{max}$", color="red", linestyle='-')
        axs[1,0].plot(var, Prvmin_list_AS, label=r"$P_{rv}^{min}$", color="green", linestyle='-')
        axs[1,0].plot(var, EDP_lv_list_AS, label=r"$EDP_{lv}$", color="black", linestyle='-')
        axs[1,0].plot(var, EDP_rv_list_AS, label=r"$EDP_{rv}$", color="brown", linestyle='-')
        axs[1,0].set_ylabel("Severe Aortic Stenosis" + "\n" + "\n" + r"Pressure [mmHg]", 
                    fontsize=12, 
                    multialignment='center')
        axs[1,0].set_xlabel(f"{var_name} {var_unit}")
        axs[1,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,1].plot(var, EDV_lv_list_AS, label=r"$EDV_{lv}$", color="green", linestyle='-')
        axs[1,1].plot(var, EDV_rv_list_AS, label=r"$EDV_{rv}$", color="brown", linestyle='-')
        axs[1,1].plot(var, SV_lv_list_AS, label=r"$SV_{lv}$", color="blue", linestyle='-')
        ax_twin11 = axs[1,1].twinx()
        ax_twin11.plot(var, EF_lv_list_AS, label=r"$EF_{lv}$", color="green", linestyle="--")
        ax_twin11.plot(var, EF_rv_list_AS, label=r"$EF_{rv}$", color="brown", linestyle="--")
        axs[1,1].set_ylabel("Volume [ml]")
        ax_twin11.set_ylabel(r"EF [\%]")
        lines11, labels11 = axs[1,1].get_legend_handles_labels()
        lines11t, labels11t = ax_twin11.get_legend_handles_labels()
        axs[1,1].legend(lines11+lines11t, labels11+labels11t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,2].plot(var, Plamax_list_AS, label=r"$P_{la}^{max}$", color="green", linestyle='-')
        axs[1,2].plot(var, Plamin_list_AS, label=r"$P_{la}^{min}$", color="red", linestyle='-')
        axs[1,2].plot(var, Pramax_list_AS, label=r"$P_{ra}^{max}$", color="blue", linestyle='-')
        axs[1,2].plot(var, Pramin_list_AS, label=r"$P_{ra}^{min}$", color="brown", linestyle='-')
        axs[1,2].set_ylabel("Pressure [mmHg]")
        axs[1,2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,3].plot(var, maxVla_list_AS, label=r"$V_{la}^{max}$", color="green", linestyle='-')
        axs[1,3].plot(var, minVla_list_AS, label=r"$V_{la}^{min}$", color="red", linestyle='-')
        axs[1,3].plot(var, maxVra_list_AS, label=r"$V_{ra}^{max}$", color="blue", linestyle='-')
        axs[1,3].plot(var, minVra_list_AS, label=r"$V_{ra}^{min}$", color="brown", linestyle='-')
        axs[1,3].set_ylabel("Volume [ml]")
        axs[1,3].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        if save_fig:
            plt.savefig(f"SVA_figures/SVA_{var_name}_2.png") 
        if show_all:
            plt.show()

        fig, axs = plt.subplots(2, 4, figsize=(15, 7), layout='constrained', sharex=True)

        axs[0,0].plot(var, Ein_total_list_noAS, label=r"$E_{in}$", color="blue", linestyle='-')
        axs[0,0].plot(var, Ein_LAD_list_noAS, label=r"$E_{in}^{LAD}$", color="brown", linestyle='-')
        axs[0,0].plot(var, Ein_LCX_list_noAS, label=r"$E_{in}^{LCX}$", color="red", linestyle='-')
        axs[0,0].plot(var, Ein_RCA_list_noAS, label=r"$E_{in}^{RCA}$", color="green", linestyle='-')
        axs[0,0].set_ylabel("Normal Aortic Valve" + "\n" + "\n" + r"Energy [J]", 
                    fontsize=12, 
                    multialignment='center') #r"$A_{max,av}^{eff} = 4.0$" +  r"$ cm^2$"
        axs[0,0].set_title(r"$E_{in} - detailed$")
        axs[0,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,1].plot(var, Eout_total_list_noAS, label=r"$E_{out}$", color="green", linestyle='-')
        axs[0,1].plot(var, Eout_lv_list_noAS, label=r"$E_{out}^{lv}$", color="brown", linestyle='-')
        axs[0,1].plot(var, Eout_rv_list_noAS, label=r"$E_{out}^{rv}$", color="blue", linestyle='-')
        axs[0,1].set_ylabel("Energy [J]")
        axs[0,1].set_title(r"$E_{out} - detailed$")
        axs[0,1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,2].plot(var, SW_lv_list_noAS, label=r"$SW_{lv}$", color="green", linestyle='-')
        axs[0,2].plot(var, SW_rv_list_noAS, label=r"$SW_{rv}$", color="red", linestyle='-')
        axs[0,2].plot(var, PE_lv_list_noAS, label=r"$PE_{lv}$", color="blue", linestyle='-')
        axs[0,2].plot(var, PE_rv_list_noAS, label=r"$PE_{rv}$", color="brown", linestyle='-')
        axs[0,2].set_ylabel("Energy [J]")
        axs[0,2].set_title("SW, PE")
        axs[0,2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,3].plot(var, CO_list_noAS, label=r"$CO$", color="green", linestyle='-')
        ax_twin03 = axs[0,3].twinx()
        ax_twin03.plot(var, LVET_list_noAS, label=r"$LVET$", color="red", linestyle="--")
        ax_twin03.plot(var, RVET_list_noAS, label=r"$RVET$", color="blue", linestyle="--")
        axs[0,3].set_ylabel("Cardiac Output [ml/min]")
        ax_twin03.set_ylabel("Time [s]")
        axs[0,3].set_title("Cardiac Output, Ejection Time")
        lines03, labels03 = axs[0,3].get_legend_handles_labels()
        lines03t, labels03t = ax_twin03.get_legend_handles_labels()
        axs[0,3].legend(lines03+lines03t, labels03+labels03t, loc='upper left', ncol=1, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,0].plot(var, Ein_total_list_AS, label=r"$E_{in}$", color="blue", linestyle='-')
        axs[1,0].plot(var, Ein_LAD_list_AS, label=r"$E_{in}^{LAD}$", color="brown", linestyle='-')
        axs[1,0].plot(var, Ein_LCX_list_AS, label=r"$E_{in}^{LCX}$", color="red", linestyle='-')
        axs[1,0].plot(var, Ein_RCA_list_AS, label=r"$E_{in}^{RCA}$", color="green", linestyle='-')
        axs[1,0].set_ylabel("Severe Aortic Stenosis" + "\n" + "\n" + r"Energy [J]", 
                    fontsize=12, 
                    multialignment='center') #r"$A_{max,av}^{eff} = 4.0$" +  r"$ cm^2$"
        axs[1,0].set_xlabel(f"{var_name} {var_unit}")
        axs[1,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,1].plot(var, Eout_total_list_AS, label=r"$E_{out}$", color="green", linestyle='-')
        axs[1,1].plot(var, Eout_lv_list_AS, label=r"$E_{out}^{lv}$", color="brown", linestyle='-')
        axs[1,1].plot(var, Eout_rv_list_AS, label=r"$E_{out}^{rv}$", color="blue", linestyle='-')
        axs[1,1].set_ylabel("Energy [J]")
        axs[1,1].set_xlabel(f"{var_name} {var_unit}")
        axs[1,1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,2].plot(var, SW_lv_list_AS, label=r"$SW_{lv}$", color="green", linestyle='-')
        axs[1,2].plot(var, SW_rv_list_AS, label=r"$SW_{rv}$", color="red", linestyle='-')
        axs[1,2].plot(var, PE_lv_list_AS, label=r"$PE_{lv}$", color="blue", linestyle='-')
        axs[1,2].plot(var, PE_rv_list_AS, label=r"$PE_{rv}$", color="brown", linestyle='-')
        axs[1,2].set_ylabel("Energy [J]")
        axs[1,2].set_xlabel(f"{var_name} {var_unit}")
        axs[1,2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,3].plot(var, CO_list_AS, label=r"$CO$", color="green", linestyle='-')
        ax_twin13 = axs[1,3].twinx()
        ax_twin13.plot(var, LVET_list_AS, label=r"$LVET$", color="red", linestyle="--")
        ax_twin13.plot(var, RVET_list_AS, label=r"$RVET$", color="blue", linestyle="--")
        axs[1,3].set_ylabel("Cardiac Output [ml/min]")
        ax_twin13.set_ylabel("Time [s]")
        axs[1,3].set_title("Cardiac Output, Ejection Time")
        lines13, labels13 = axs[1,3].get_legend_handles_labels()
        lines13t, labels13t = ax_twin13.get_legend_handles_labels()
        axs[1,3].legend(lines13+lines13t, labels13+labels13t, loc='upper left', ncol=1, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        if save_fig:
            plt.savefig(f"SVA_figures/SVA_{var_name}_3.png") 
        if show_all:
            plt.show()

        fig, axs = plt.subplots(2, 4, figsize=(15, 7), layout='constrained', sharex=True)

        axs[0,0].plot(var, Pv_syst_max_list_noAS, label=r"$P_{v,max}^{syst}$", color="blue", linestyle='-')
        axs[0,0].plot(var, Pv_syst_min_list_noAS, label=r"$P_{v,min}^{syst}$", color="brown", linestyle='-')
        axs[0,0].plot(var, Pv_syst_mean_list_noAS, label=r"$CVP$", color="red", linestyle='-')
        axs[0,0].set_ylabel("Normal Aortic Valve" + "\n" + "\n" + r"Pressure [mmHg]", 
                    fontsize=12, 
                    multialignment='center') #r"$A_{max,av}^{eff} = 4.0$" +  r"$ cm^2$"
        axs[0,0].set_title(r"Venous Systemic Circulation")
        axs[0,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,1].plot(var, Pv_pulm_max_list_noAS, label=r"$P_{v,max}^{pulm}$", color="blue", linestyle='-')
        axs[0,1].plot(var, Pv_pulm_min_list_noAS, label=r"$P_{v,min}^{pulm}$", color="brown", linestyle='-')
        axs[0,1].plot(var, Pv_pulm_mean_list_noAS, label=r"$PVP$", color="red", linestyle='-')
        axs[0,1].set_ylabel("Pressure [mmHg]") 
        axs[0,1].set_title(r"Venous Pulmonary Circulation")
        axs[0,1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,2].plot(var, Ees_lv_list_noAS, label=r"$E_{es}$", color="red", linestyle='-')
        axs[0,2].plot(var, Ea_lv_list_noAS, label=r"$E_{a}$", color="blue", linestyle='-')
        ax_twin02 = axs[0,2].twinx()
        ax_twin02.plot(var, VAC_lv_list_noAS, label=r"$VAC$", color="green", linestyle="--")
        axs[0,2].set_ylabel("Elastance [mmHg/ml]")
        ax_twin02.set_ylabel("VAC [-]")
        axs[0,2].set_title("Left Ventriculo-Arterial Coupling")
        lines02, labels02 = axs[0,2].get_legend_handles_labels()
        lines02t, labels02t = ax_twin02.get_legend_handles_labels()
        axs[0,2].legend(lines02+lines02t, labels02+labels02t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[0,3].plot(var, Ees_rv_list_noAS, label=r"$E_{es}$", color="red", linestyle='-')
        axs[0,3].plot(var, Ea_rv_list_noAS, label=r"$E_{a}$", color="blue", linestyle='-')
        ax_twin03 = axs[0,3].twinx()
        ax_twin03.plot(var, VAC_rv_list_noAS, label=r"$VAC$", color="green", linestyle="--")
        axs[0,3].set_ylabel("Elastance [mmHg/ml]]")
        ax_twin03.set_ylabel("VAC [-]")
        axs[0,3].set_title("Right Ventriculo-Arterial Coupling")
        lines03, labels03 = axs[0,3].get_legend_handles_labels()
        lines03t, labels03t = ax_twin03.get_legend_handles_labels()
        axs[0,3].legend(lines03+lines03t, labels03+labels03t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,0].plot(var, Pv_syst_max_list_AS, label=r"$P_{v,max}^{syst}$", color="blue", linestyle='-')
        axs[1,0].plot(var, Pv_syst_min_list_AS, label=r"$P_{v,min}^{syst}$", color="brown", linestyle='-')
        axs[1,0].plot(var, Pv_syst_mean_list_AS, label=r"$CVP$", color="red", linestyle='-')
        axs[1,0].set_ylabel("Severe Aortic Stenosis" + "\n" + "\n" + r"Pressure [mmHg]", 
                    fontsize=12, 
                    multialignment='center') #r"$A_{max,av}^{eff} = 4.0$" +  r"$ cm^2$"
        axs[1,0].set_xlabel(f"{var_name} {var_unit}")
        axs[1,0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,1].plot(var, Pv_pulm_max_list_AS, label=r"$P_{v,max}^{pulm}$", color="blue", linestyle='-')
        axs[1,1].plot(var, Pv_pulm_min_list_AS, label=r"$P_{v,min}^{pulm}$", color="brown", linestyle='-')
        axs[1,1].plot(var, Pv_pulm_mean_list_AS, label=r"$PVP$", color="red", linestyle='-')
        axs[1,1].set_ylabel("Pressure [mmHg]") 
        axs[1,1].set_xlabel(f"{var_name} {var_unit}")
        axs[1,1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,2].plot(var, Ees_lv_list_AS, label=r"$E_{es}$", color="red", linestyle='-')
        axs[1,2].plot(var, Ea_lv_list_AS, label=r"$E_{a}$", color="blue", linestyle='-')
        ax_twin12 = axs[1,2].twinx()
        ax_twin12.plot(var, VAC_lv_list_AS, label=r"$VAC$", color="green", linestyle="--")
        axs[1,2].set_ylabel("Elastance [mmHg/ml]]")
        ax_twin12.set_ylabel("VAC [-]")
        axs[1,2].set_xlabel(f"{var_name} {var_unit}")
        lines12, labels12 = axs[1,2].get_legend_handles_labels()
        lines12t, labels12t = ax_twin12.get_legend_handles_labels()
        axs[1,2].legend(lines12+lines12t, labels12+labels12t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1,3].plot(var, Ees_rv_list_AS, label=r"$E_{es}$", color="red", linestyle='-')
        axs[1,3].plot(var, Ea_rv_list_AS, label=r"$E_{a}$", color="blue", linestyle='-')
        ax_twin13 = axs[1,3].twinx()
        ax_twin13.plot(var, VAC_rv_list_AS, label=r"$VAC$", color="green", linestyle="--")
        axs[1,3].set_ylabel("Elastance [mmHg/ml]]")
        ax_twin03.set_ylabel("VAC [-]")
        axs[1,3].set_xlabel(f"{var_name} {var_unit}")
        lines13, labels13 = axs[1,3].get_legend_handles_labels()
        lines13t, labels13t = ax_twin13.get_legend_handles_labels()
        axs[1,3].legend(lines13+lines13t, labels13+labels13t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        if save_fig:
            plt.savefig(f"SVA_figures/SVA_{var_name}_4.png") 
        if show_all:
            plt.show()



def SVA_plot_from_csv_forSao(path_noAS, save_fig=False, complete_plot=False):
    # loading data
    df_noAS = pd.read_csv(path_noAS)

    STK_syst_list_noAS = df_noAS['STK_syst_list'].tolist()
    DTK_syst_list_noAS = df_noAS['DTK_syst_list'].tolist()
    MAP_syst_list_noAS = df_noAS['MAP_syst_list'].tolist()
    Plvmax_list_noAS = df_noAS['Plvmax_list'].tolist()
    SV_lv_list_noAS = df_noAS['SV_lv_list'].tolist()
    EDV_lv_list_noAS = df_noAS['EDV_lv_list'].tolist()
    STK_pulm_list_noAS = df_noAS['STK_pulm_list'].tolist()
    DTK_pulm_list_noAS = df_noAS['DTK_pulm_list'].tolist()
    MAP_pulm_list_noAS = df_noAS['MAP_pulm_list'].tolist()
    Prvmax_list_noAS = df_noAS['Prvmax_list'].tolist()
    SV_rv_list_noAS = df_noAS['SV_rv_list'].tolist()
    EDV_rv_list_noAS = df_noAS['EDV_rv_list'].tolist()
    maxVla_list_noAS = df_noAS['maxVla_list'].tolist()
    minVla_list_noAS = df_noAS['minVla_list'].tolist()
    maxVra_list_noAS = df_noAS['maxVra_list'].tolist()
    minVra_list_noAS = df_noAS['minVra_list'].tolist()
    Eout_total_list_noAS = df_noAS['Eout_total_list'].tolist()
    Eout_lv_list_noAS = df_noAS['Eout_lv_list'].tolist()
    Eout_rv_list_noAS = df_noAS['Eout_rv_list'].tolist()
    SW_lv_list_noAS = df_noAS['SW_lv_list'].tolist()
    PE_lv_list_noAS = df_noAS['PE_lv_list'].tolist()
    SW_rv_list_noAS = df_noAS['SW_rv_list'].tolist()
    PE_rv_list_noAS = df_noAS['PE_rv_list'].tolist()
    Ein_total_list_noAS = df_noAS['Ein_total_list'].tolist()
    Ein_LAD_list_noAS = df_noAS['Ein_LAD_list'].tolist()
    Ein_LCX_list_noAS = df_noAS['Ein_LCX_list'].tolist()
    Ein_RCA_list_noAS = df_noAS['Ein_RCA_list'].tolist()
    Ees_lv_list_noAS = df_noAS['Ees_lv_list'].tolist()
    Pes_lv_list_noAS = df_noAS['Pes_lv_list'].tolist()
    LVET_list_noAS = df_noAS['LVET_list'].tolist()
    EF_lv_list_noAS = df_noAS["EF_lv_list"].tolist()
    EF_rv_list_noAS = df_noAS["EF_rv_list"].tolist()
    Plvmin_list_noAS = df_noAS["Plvmin_list"].tolist()
    EDP_lv_list_noAS = df_noAS["EDP_lv_list"].tolist()
    Pv_syst_max_list_noAS = df_noAS["Pv_syst_max_list"].tolist()
    Pv_syst_min_list_noAS = df_noAS["Pv_syst_min_list"].tolist()
    Pv_syst_mean_list_noAS = df_noAS["Pv_syst_mean_list"].tolist()
    Prvmin_list_noAS = df_noAS["Prvmin_list"].tolist()
    EDP_rv_list_noAS = df_noAS["EDP_rv_list"].tolist()
    Pv_pulm_max_list_noAS = df_noAS["Pv_pulm_max_list"].tolist()
    Pv_pulm_min_list_noAS = df_noAS["Pv_pulm_min_list"].tolist()
    Pv_pulm_mean_list_noAS = df_noAS["Pv_pulm_mean_list"].tolist()
    Plamax_list_noAS = df_noAS["Plamax_list"].tolist()
    Pramax_list_noAS = df_noAS["Pramax_list"].tolist()
    Plamin_list_noAS = df_noAS["Plamin_list"].tolist()
    Pramin_list_noAS = df_noAS["Pramin_list"].tolist()
    RVET_list_noAS = df_noAS["RVET_list"].tolist()
    Ees_rv_list_noAS = df_noAS["Ees_rv_list"].tolist()
    CO_list_noAS = df_noAS["CO_list"].tolist()
    Ea_lv_list_noAS = df_noAS["Ea_lv_list"].tolist()
    Ea_rv_list_noAS = df_noAS["Ea_lv_list"].tolist()
    VAC_lv_list_noAS = df_noAS["VAC_lv_list"].tolist()
    VAC_rv_list_noAS = df_noAS["VAC_rv_list"].tolist()

    var = df_noAS["var"].tolist()
    var_name = df_noAS["var_name"][0]
    var_unit = df_noAS["var_unit"][0]

    


    

    centimeter = 1/2.54

    fig, axs = plt.subplots(1, 3, figsize=(20*centimeter, 5.5*centimeter), layout='constrained', sharex=True)

    # axs[0,0].plot(var, y, label=r"$SBP^{syst}$", color="blue", linestyle='-')
    # axs[0,0].plot(var, y, label=r"$DBP^{syst}$", color="brown", linestyle='-')
    axs[0].plot(var, MAP_syst_list_noAS, label=r"$MAP_{syst}$", color="red", linestyle='-')
    axs[0].plot(var, Plvmax_list_noAS, label=r"$P_{lv}^{max}$", color="green", linestyle='-')
    ax_twin00 = axs[0].twinx()
    ax_twin00.plot(var, EDV_lv_list_noAS, label=r"$EDV_{lv}$", color="blue", linestyle='--')
    ax_twin00.plot(var, SV_lv_list_noAS, label=r"$SV_{lv}$", color="purple", linestyle='--')
    axs[0].set_ylabel(r"[pressure]")
    axs[0].set_xlabel(f"{var_name} {var_unit}") 
    axs[0].set_title("Systemic Circulation", fontsize=10)
    lines1, labels1 = axs[0].get_legend_handles_labels()
    lines2, labels2 = ax_twin00.get_legend_handles_labels()
    axs[0].legend(lines1 + lines2, labels1 + labels2, loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

    # axs[0,1].plot(var, y, label=r"$SBP^{pulm}$", color="blue", linestyle='-')
    # axs[0,1].plot(var, y, label=r"$DBP^{pulm}$", color="brown", linestyle='-')
    axs[1].plot(var, MAP_pulm_list_noAS, label=r"$MAP_{pulm}$", color="red", linestyle='-')
    axs[1].plot(var, Prvmax_list_noAS, label=r"$P_{rv}^{max}$", color="green", linestyle='-')
    ax_twin01 = axs[1].twinx()
    ax_twin01.plot(var, EDV_rv_list_noAS, label=r"$EDV_{rv}$", color="blue", linestyle='--')
    ax_twin01.plot(var, SV_rv_list_noAS, label=r"$SV_{rv}$", color="purple", linestyle='--')
    #axs[0,1].set_ylabel("Pressure [mmHg]")
    ax_twin01.set_ylabel(f"Volume [ml]")
    axs[1].set_xlabel(f"{var_name} {var_unit}")
    axs[1].set_title("Pulmonary Circulation")
    lines1, labels1 = axs[1].get_legend_handles_labels()
    lines2, labels2 = ax_twin01.get_legend_handles_labels()
    axs[1].legend(lines1 + lines2, labels1 + labels2, loc='upper left', ncol=2, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'


    axs[2].plot(var, Ein_total_list_noAS, label=r"$E_{in}$", color="blue", linestyle='-')
    axs[2].plot(var, Eout_total_list_noAS, label=r"$E_{out}$", color="y", linestyle='-')
    axs[2].plot(var, np.array(Ein_total_list_noAS)-np.array(Eout_total_list_noAS), label=r"$\Delta E$", color="red")
    axs[2].set_ylabel("Energy [J/beat]")
    axs[2].set_xlabel(f"{var_name} {var_unit}")
    axs[2].set_title("Energies")
    axs[2].legend(loc='upper left', ncol=3, frameon=True, columnspacing=0.8, handlelength=0.8, fontsize=6, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

    if save_fig:
        plt.savefig(f"SVA_figures/final_{var_name}.png", dpi=400) 
    plt.show()

    if complete_plot:
        fig, axs = plt.subplots(1, 4, figsize=(15, 4), layout='constrained', sharex=True)

        axs[0].plot(var, STK_syst_list_noAS, label=r"$SBP^{syst}$", color="blue", linestyle='-')
        axs[0].plot(var, DTK_syst_list_noAS, label=r"$DBP^{syst}$", color="brown", linestyle='-')
        axs[0].plot(var, MAP_syst_list_noAS, label=r"$MAP^{syst}$", color="red", linestyle='-')
        axs[0].plot(var, Pv_syst_mean_list_noAS, label=r"$CVP$", color="pink", linestyle='-')
        axs[0].plot(var, Plvmax_list_noAS, label=r"$P_{lv}^{max}$", color="green", linestyle='-')
        axs[0].set_ylabel(r"Pressure [mmHg]") 
        axs[0].set_title("Pressures - Systemic Circulation")
        axs[0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1].plot(var, STK_pulm_list_noAS, label=r"$SBP^{pulm}$", color="blue", linestyle='-')
        axs[1].plot(var, DTK_pulm_list_noAS, label=r"$DBP^{pulm}$", color="brown", linestyle='-')
        axs[1].plot(var, MAP_pulm_list_noAS, label=r"$MAP^{pulm}$", color="red", linestyle='-')
        axs[1].plot(var, Pv_pulm_mean_list_noAS, label=r"$PVP$", color="pink", linestyle='-')
        axs[1].plot(var, Prvmax_list_noAS, label=r"$P_{rv,max}$", color="green", linestyle='-')
        axs[1].set_ylabel("Pressure [mmHg]")
        axs[1].set_title("Pressures - Pulmonary Circulation")
        axs[1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[2].plot(var, EDV_lv_list_noAS, label=r"$EDV_{lv}$", color="green", linestyle='-')
        axs[2].plot(var, EDV_rv_list_noAS, label=r"$EDV_{rv}$", color="brown", linestyle='-')
        axs[2].plot(var, SV_lv_list_noAS, label=r"$SV_{lv}$", color="blue", linestyle='-')
        ax_twin1 = axs[2].twinx()
        ax_twin1.plot(var, EF_lv_list_noAS, label=r"$EF_{lv}$", color="green", linestyle="--")
        ax_twin1.plot(var, EF_rv_list_noAS, label=r"$EF_{rv}$", color="brown", linestyle="--")
        axs[2].set_ylabel("Volume [ml]")
        ax_twin1.set_ylabel(r"EF [\%]")
        axs[2].set_title("Volumes")
        lines1, labels1 = axs[2].get_legend_handles_labels()
        lines2, labels2 = ax_twin1.get_legend_handles_labels()
        axs[2].legend(lines1 + lines2, labels1 + labels2, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[3].plot(var, Ein_total_list_noAS, label=r"$E_{in}$", color="green", linestyle='-')
        axs[3].plot(var, Eout_total_list_noAS, label=r"$E_{out}$", color="red", linestyle='-')
        axs[3].plot(var, np.array(Ein_total_list_noAS)-np.array(Eout_total_list_noAS), label=r"$\Delta E$", color="blue")
        # ax_twin03 = axs[3].twinx()
        # ax_twin03.plot(var, VAC_lv_list_noAS, label=r"$VAC$", color="black", linestyle="--")
        axs[3].set_ylabel("Energy [J/beat]")
        # ax_twin03.set_ylabel("VAC [-]")
        axs[3].set_title("Energies")
        # lines03, labels03 = axs[3].get_legend_handles_labels()
        #lines03t, labels03t = ax_twin03.get_legend_handles_labels()
        # axs[3].legend(lines03 + lines03t, labels03 + labels03t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'
        axs[3].legend(loc='upper left', ncol=3, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        
        
        if save_fig:
            plt.savefig(f"SVA_figures/SVA_{var_name}_1.png") 
        plt.show()

        fig, axs = plt.subplots(1, 4, figsize=(15, 4), layout='constrained', sharex=True)

        axs[0].plot(var, Plvmax_list_noAS, label=r"$P_{lv}^{max}$", color="blue", linestyle='-')
        axs[0].plot(var, Plvmin_list_noAS, label=r"$P_{lv}^{min}$", color="brown", linestyle='-')
        axs[0].plot(var, Prvmax_list_noAS, label=r"$P_{rv}^{max}$", color="red", linestyle='-')
        axs[0].plot(var, Prvmin_list_noAS, label=r"$P_{rv}^{min}$", color="green", linestyle='-')
        axs[0].plot(var, EDP_lv_list_noAS, label=r"$EDP_{lv}$", color="black", linestyle='-')
        axs[0].plot(var, EDP_rv_list_noAS, label=r"$EDP_{rv}$", color="brown", linestyle='-')
        axs[0].set_ylabel(r"Pressure [mmHg]") #r"$A_{max,av}^{eff} = 4.0$" +  r"$ cm^2$"
        axs[0].set_title("Ventricles - Pressures")
        axs[0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1].plot(var, EDV_lv_list_noAS, label=r"$EDV_{lv}$", color="green", linestyle='-')
        axs[1].plot(var, EDV_rv_list_noAS, label=r"$EDV_{rv}$", color="brown", linestyle='-')
        axs[1].plot(var, SV_lv_list_noAS, label=r"$SV_{lv}$", color="blue", linestyle='-')
        ax_twin01 = axs[1].twinx()
        ax_twin01.plot(var, EF_lv_list_noAS, label=r"$EF_{lv}$", color="green", linestyle="--")
        ax_twin01.plot(var, EF_rv_list_noAS, label=r"$EF_{rv}$", color="brown", linestyle="--")
        axs[1].set_ylabel("Volume [ml]")
        ax_twin01.set_ylabel(r"EF [\%]")
        axs[1].set_title("Ventricles - Volumes")
        lines01, labels01 = axs[1].get_legend_handles_labels()
        lines01t, labels01t = ax_twin01.get_legend_handles_labels()
        axs[1].legend(lines01+lines01t, labels01+labels01t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[2].plot(var, Plamax_list_noAS, label=r"$P_{la}^{max}$", color="green", linestyle='-')
        axs[2].plot(var, Plamin_list_noAS, label=r"$P_{la}^{min}$", color="red", linestyle='-')
        axs[2].plot(var, Pramax_list_noAS, label=r"$P_{ra}^{max}$", color="blue", linestyle='-')
        axs[2].plot(var, Pramin_list_noAS, label=r"$P_{ra}^{min}$", color="brown", linestyle='-')
        axs[2].set_ylabel("Pressure [mmHg]")
        axs[2].set_title("Atria - Pressures")
        axs[2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[3].plot(var, maxVla_list_noAS, label=r"$V_{la}^{max}$", color="green", linestyle='-')
        axs[3].plot(var, minVla_list_noAS, label=r"$V_{la}^{min}$", color="red", linestyle='-')
        axs[3].plot(var, maxVra_list_noAS, label=r"$V_{ra}^{max}$", color="blue", linestyle='-')
        axs[3].plot(var, minVra_list_noAS, label=r"$V_{ra}^{min}$", color="brown", linestyle='-')
        axs[3].set_ylabel("Volume [ml]")
        axs[3].set_title("Atria - Volumes")
        axs[3].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        
        if save_fig:
            plt.savefig(f"SVA_figures/SVA_{var_name}_2.png") 
        plt.show()

        fig, axs = plt.subplots(1, 4, figsize=(15, 4), layout='constrained', sharex=True)

        axs[0].plot(var, Ein_total_list_noAS, label=r"$E_{in}$", color="blue", linestyle='-')
        axs[0].plot(var, Ein_LAD_list_noAS, label=r"$E_{in}^{LAD}$", color="brown", linestyle='-')
        axs[0].plot(var, Ein_LCX_list_noAS, label=r"$E_{in}^{LCX}$", color="red", linestyle='-')
        axs[0].plot(var, Ein_RCA_list_noAS, label=r"$E_{in}^{RCA}$", color="green", linestyle='-')
        axs[0].set_ylabel(r"Energy [J]") #r"$A_{max,av}^{eff} = 4.0$" +  r"$ cm^2$"
        axs[0].set_title(r"$E_{in} - detailed$")
        axs[0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1].plot(var, Eout_total_list_noAS, label=r"$E_{out}$", color="green", linestyle='-')
        axs[1].plot(var, Eout_lv_list_noAS, label=r"$E_{out}^{lv}$", color="brown", linestyle='-')
        axs[1].plot(var, Eout_rv_list_noAS, label=r"$E_{out}^{rv}$", color="blue", linestyle='-')
        axs[1].set_ylabel("Energy [J]")
        axs[1].set_title(r"$E_{out} - detailed$")
        axs[1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[2].plot(var, SW_lv_list_noAS, label=r"$SW_{lv}$", color="green", linestyle='-')
        axs[2].plot(var, SW_rv_list_noAS, label=r"$SW_{rv}$", color="red", linestyle='-')
        axs[2].plot(var, PE_lv_list_noAS, label=r"$PE_{lv}$", color="blue", linestyle='-')
        axs[2].plot(var, PE_rv_list_noAS, label=r"$PE_{rv}$", color="brown", linestyle='-')
        axs[2].set_ylabel("Energy [J]")
        axs[2].set_title("SW, PE")
        axs[2].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[3].plot(var, CO_list_noAS, label=r"$CO$", color="green", linestyle='-')
        ax_twin03 = axs[3].twinx()
        ax_twin03.plot(var, LVET_list_noAS, label=r"$LVET$", color="red", linestyle="--")
        ax_twin03.plot(var, RVET_list_noAS, label=r"$RVET$", color="blue", linestyle="--")
        axs[3].set_ylabel("Cardiac Output [ml/min]")
        ax_twin03.set_ylabel("Time [s]")
        axs[3].set_title("Cardiac Output, Ejection Time")
        lines03, labels03 = axs[3].get_legend_handles_labels()
        lines03t, labels03t = ax_twin03.get_legend_handles_labels()
        axs[3].legend(lines03+lines03t, labels03+labels03t, loc='upper left', ncol=1, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        
        if save_fig:
            plt.savefig(f"SVA_figures/SVA_{var_name}_3.png") 
        plt.show()

        fig, axs = plt.subplots(1, 4, figsize=(15, 4), layout='constrained', sharex=True)

        axs[0].plot(var, Pv_syst_max_list_noAS, label=r"$P_{v,max}^{syst}$", color="blue", linestyle='-')
        axs[0].plot(var, Pv_syst_min_list_noAS, label=r"$P_{v,min}^{syst}$", color="brown", linestyle='-')
        axs[0].plot(var, Pv_syst_mean_list_noAS, label=r"$CVP$", color="red", linestyle='-')
        axs[0].set_ylabel("Normal Aortic Valve" + "\n" + "\n" + r"Pressure [mmHg]", 
                    fontsize=12, 
                    multialignment='center') #r"$A_{max,av}^{eff} = 4.0$" +  r"$ cm^2$"
        axs[0].set_title(r"Venous Systemic Circulation")
        axs[0].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[1].plot(var, Pv_pulm_max_list_noAS, label=r"$P_{v,max}^{pulm}$", color="blue", linestyle='-')
        axs[1].plot(var, Pv_pulm_min_list_noAS, label=r"$P_{v,min}^{pulm}$", color="brown", linestyle='-')
        axs[1].plot(var, Pv_pulm_mean_list_noAS, label=r"$PVP$", color="red", linestyle='-')
        axs[1].set_ylabel("Pressure [mmHg]") 
        axs[1].set_title(r"Venous Pulmonary Circulation")
        axs[1].legend(loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[2].plot(var, Ees_lv_list_noAS, label=r"$E_{es}$", color="red", linestyle='-')
        axs[2].plot(var, Ea_lv_list_noAS, label=r"$E_{a}$", color="blue", linestyle='-')
        ax_twin02 = axs[2].twinx()
        ax_twin02.plot(var, VAC_lv_list_noAS, label=r"$VAC$", color="green", linestyle="--")
        axs[2].set_ylabel("Elastance [mmHg/ml]")
        ax_twin02.set_ylabel("VAC [-]")
        axs[2].set_title("Left Ventriculo-Arterial Coupling")
        lines02, labels02 = axs[2].get_legend_handles_labels()
        lines02t, labels02t = ax_twin02.get_legend_handles_labels()
        axs[2].legend(lines02+lines02t, labels02+labels02t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        axs[3].plot(var, Ees_rv_list_noAS, label=r"$E_{es}$", color="red", linestyle='-')
        axs[3].plot(var, Ea_rv_list_noAS, label=r"$E_{a}$", color="blue", linestyle='-')
        ax_twin03 = axs[3].twinx()
        ax_twin03.plot(var, VAC_rv_list_noAS, label=r"$VAC$", color="green", linestyle="--")
        axs[3].set_ylabel("Elastance [mmHg/ml]]")
        ax_twin03.set_ylabel("VAC [-]")
        axs[3].set_title("Right Ventriculo-Arterial Coupling")
        lines03, labels03 = axs[3].get_legend_handles_labels()
        lines03t, labels03t = ax_twin03.get_legend_handles_labels()
        axs[3].legend(lines03+lines03t, labels03+labels03t, loc='upper left', ncol=2, frameon=True, columnspacing=1.0, handlelength=1.0, fontsize=8, edgecolor="black", facecolor="white", fancybox=True, framealpha=0.9, shadow=False) # loc='upper left'

        
        if save_fig:
            plt.savefig(f"SVA_figures/SVA_{var_name}_4.png") 
        plt.show()