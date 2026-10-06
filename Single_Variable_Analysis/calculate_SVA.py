import sys
from pathlib import Path
project_root = Path(__file__).resolve().parent.parent
sys.path.append(str(project_root))
from cycles import *
from plot_functions import *
import json
import numpy as np
from copy import deepcopy
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','ieee'])
# plt.rcParams['axes.spines.top'] = False
# plt.rcParams['axes.spines.right'] = False
plt.rcParams['xtick.top'] = False
plt.rcParams['ytick.right'] = False
plt.rcParams.update({
    "font.family": "serif",   
    "font.serif": ["Times"],  
    "font.size":12})          
plt.rcParams['patch.linewidth'] = 0.5

def save_to_csv_SVA(list_of_results_SVA, var_list, file_name):
    STK_syst_list = []
    DTK_syst_list = []
    MAP_syst_list = []
    Plvmax_list = []
    SV_lv_list = []
    EDV_lv_list = []

    STK_pulm_list = []
    DTK_pulm_list = []
    MAP_pulm_list = []
    Prvmax_list = []
    SV_rv_list = []
    EDV_rv_list = []

    maxVla_list = []
    minVla_list = []
    maxVra_list = []
    minVra_list = []

    Eout_total_list = []
    Eout_lv_list = []
    Eout_rv_list = []
    SW_lv_list = []
    PE_lv_list = []
    SW_rv_list = []
    PE_rv_list = []
    Ein_total_list = []
    Ein_LAD_list = []
    Ein_LCX_list = []
    Ein_RCA_list = []

    Ees_lv_list = []
    Pes_lv_list = []

    LVET_list = []

    EF_lv_list = []
    EF_rv_list = []

    Plvmin_list = []
    EDP_lv_list = []
    Pv_syst_max_list = []
    Pv_syst_min_list = []
    Pv_syst_mean_list = []
    Prvmin_list = []
    EDP_rv_list = []
    Pv_pulm_max_list = []
    Pv_pulm_min_list = []
    Pv_pulm_mean_list = []
    Plamax_list = []
    Pramax_list = []
    Plamin_list = []
    Pramin_list = []
    RVET_list = []
    Ees_rv_list = []
    CO_list = []
    Ea_lv_list = []
    Ea_rv_list = []
    VAC_lv_list = []
    VAC_rv_list = []

    var, var_name, var_unit = var_list

    for res in list_of_results_SVA:
        STK_syst_list.append(res["output_values"]["STK_syst"])
        DTK_syst_list.append(res["output_values"]["DTK_syst"])
        MAP_syst_list.append(res["output_values"]["MAP_syst"])
        Plvmax_list.append(res["output_values"]["Plvmax"])
        Plvmin_list.append(res["output_values"]["Plvmin"])
        EDP_lv_list.append(res["output_values"]["EDP_lv"])
        SV_lv_list.append(res["output_values"]["SV_lv"])
        EDV_lv_list.append(res["output_values"]["EDV_lv"])
        Pv_syst_max_list.append(res["output_values"]["Pv_syst_max"])
        Pv_syst_min_list.append(res["output_values"]["Pv_syst_min"]) 
        Pv_syst_mean_list.append(res["output_values"]["Pv_syst_mean"]) 


        STK_pulm_list.append(res["output_values"]["STK_pulm"])
        DTK_pulm_list.append(res["output_values"]["DTK_pulm"])
        MAP_pulm_list.append(res["output_values"]["MAP_pulm"])
        Prvmax_list.append(res["output_values"]["Prvmax"])
        Prvmin_list.append(res["output_values"]["Prvmin"]) 
        EDP_rv_list.append(res["output_values"]["EDP_rv"]) 
        SV_rv_list.append(res["output_values"]["SV_rv"])
        EDV_rv_list.append(res["output_values"]["EDV_rv"])
        Pv_pulm_max_list.append(res["output_values"]["Pv_pulm_max"]) 
        Pv_pulm_min_list.append(res["output_values"]["Pv_pulm_min"]) 
        Pv_pulm_mean_list.append(res["output_values"]["Pv_pulm_mean"]) 

        maxVla_list.append(res["output_values"]["maxVla"])
        minVla_list.append(res["output_values"]["minVla"])
        maxVra_list.append(res["output_values"]["maxVra"])
        minVra_list.append(res["output_values"]["minVra"])
        Plamax_list.append(res["output_values"]["Plamax"]) 
        Pramax_list.append(res["output_values"]["Pramax"]) 
        Plamin_list.append(res["output_values"]["Plamin"]) 
        Pramin_list.append(res["output_values"]["Pramin"]) 

        Eout_total_list.append(res["output_values"]["Eout"])
        Eout_lv_list.append(res["output_values"]["Eout_lv"])
        Eout_rv_list.append(res["output_values"]["Eout_rv"])
        SW_lv_list.append(res["output_values"]["SW_lv"])
        SW_rv_list.append(res["output_values"]["SW_rv"])
        PE_lv_list.append(res["output_values"]["PE_lv"])
        PE_rv_list.append(res["output_values"]["PE_rv"])
        Ein_total_list.append(res["output_values"]["Ein"])
        Ein_LAD_list.append(res["output_values"]["Ein_LAD"])
        Ein_LCX_list.append(res["output_values"]["Ein_LCX"])
        Ein_RCA_list.append(res["output_values"]["Ein_RCA"])

        LVET_list.append(res["output_values"]["LVET"])
        RVET_list.append(res["output_values"]["RVET"]) 
        Ees_lv_list.append(res["output_values"]["Ees_lv"])
        Ees_rv_list.append(res["output_values"]["Ees_rv"]) 
        Pes_lv_list.append(res["output_values"]["Pes_lv"])

        EF_lv_list.append(res["output_values"]["EF_lv"])
        EF_rv_list.append(res["output_values"]["EF_rv"])
        CO_list.append(res["output_values"]["CO"])
        Ea_lv_list.append(res["output_values"]["Ea_lv"])
        Ea_rv_list.append(res["output_values"]["Ea_rv"])
        VAC_lv_list.append(res["output_values"]["VAC_lv"])
        VAC_rv_list.append(res["output_values"]["VAC_rv"])

        

    data = {
        "STK_syst_list": STK_syst_list,
        "DTK_syst_list": DTK_syst_list,
        "MAP_syst_list": MAP_syst_list,
        "Plvmax_list": Plvmax_list,
        "SV_lv_list": SV_lv_list,
        "EDV_lv_list": EDV_lv_list,
        "STK_pulm_list": STK_pulm_list,
        "DTK_pulm_list": DTK_pulm_list,
        "MAP_pulm_list": MAP_pulm_list,
        "Prvmax_list": Prvmax_list,
        "SV_rv_list": SV_rv_list,
        "EDV_rv_list": EDV_rv_list,
        "maxVla_list": maxVla_list,
        "minVla_list": minVla_list,
        "maxVra_list": maxVra_list,
        "minVra_list": minVra_list,
        "Eout_total_list": Eout_total_list,
        "Eout_lv_list": Eout_lv_list,
        "Eout_rv_list": Eout_rv_list,
        "SW_lv_list": SW_lv_list,
        "PE_lv_list": PE_lv_list,
        "SW_rv_list": SW_rv_list,
        "PE_rv_list": PE_rv_list,
        "Ein_total_list": Ein_total_list,
        "Ein_LAD_list": Ein_LAD_list,
        "Ein_LCX_list": Ein_LCX_list,
        "Ein_RCA_list": Ein_RCA_list,
        "Ees_lv_list": Ees_lv_list,
        "Pes_lv_list": Pes_lv_list,
        "LVET_list": LVET_list,
        "EF_lv_list": EF_lv_list,
        "EF_rv_list": EF_rv_list,
        "var": var,
        "var_name": var_name,
        "var_unit": var_unit,
        "Plvmin_list": Plvmin_list,
        "EDP_lv_list": EDP_lv_list,
        "Pv_syst_max_list": Pv_syst_max_list,
        "Pv_syst_min_list": Pv_syst_min_list,
        "Pv_syst_mean_list": Pv_syst_mean_list,
        "Prvmin_list": Prvmin_list,
        "EDP_rv_list": EDP_rv_list,
        "Pv_pulm_max_list": Pv_pulm_max_list,
        "Pv_pulm_min_list": Pv_pulm_min_list,
        "Pv_pulm_mean_list": Pv_pulm_mean_list,
        "Plamax_list": Plamax_list,
        "Pramax_list": Pramax_list,
        "Plamin_list": Plamin_list,
        "Pramin_list": Pramin_list,
        "RVET_list": RVET_list,
        "Ees_rv_list": Ees_rv_list,
        "CO_list": CO_list,
        "Ea_lv_list": Ea_lv_list,
        "Ea_rv_list": Ea_rv_list,
        "VAC_lv_list": VAC_lv_list,
        "VAC_rv_list": VAC_rv_list
        }

    df = pd.DataFrame(data)
    df.to_csv(file_name, index=False)


def run_single_simulation(base_params, path1, path2, j_val):
    params = deepcopy(base_params)
    params[path1][path2] = j_val
    res, icycle = multiple_cycles(params, 500, True, store_all=False)
    return res, icycle

def perform_SVA_multiprocessing(var_list=None,
                var_path=None,
                plot=True,
                save_fig=True,
                parameters_input_path="../input_values/parameters_input.json",
                n_jobs=-1): 
    
    var, var_name, var_unit = var_list
    var_path1, var_path2 = var_path

    # --- BLOCK 1: noAS ---
    with open(parameters_input_path, "r") as f:
        parameters_input = json.load(f)

    print(f"Running noAS simulations on {n_jobs} cores...")
    results_noAS = Parallel(n_jobs=n_jobs)(
        delayed(run_single_simulation)(parameters_input, var_path1, var_path2, j) 
        for j in var
    )
    
    list_of_results_noAS, stability_cyclenumber_noAS = zip(*results_noAS)
    
    save_to_csv_SVA(list(list_of_results_noAS), var_list, f"SVA_data/{var_name}_noAS.csv")


    # --- BLOCK 2: AS ---
    with open(parameters_input_path, "r") as f:
        parameters_input_AS = json.load(f)
    parameters_input_AS["valves"]["Amax_av"] = 0.5

    print(f"Running AS simulations on {n_jobs} cores...")
    results_AS = Parallel(n_jobs=n_jobs)(
        delayed(run_single_simulation)(parameters_input_AS, var_path1, var_path2, j) 
        for j in var
    )
    
    list_of_results_AS, stability_cyclenumber_AS = zip(*results_AS)
    save_to_csv_SVA(list(list_of_results_AS), var_list, f"SVA_data/{var_name}_AS.csv")
    if plot:
        SVA_plot_from_csv(path_noAS=f"SVA_data/{var_name}_noAS.csv",
                        path_AS=f"SVA_data/{var_name}_AS.csv",
                        save_fig=save_fig,
                        complete_plot=True,
                        )

        plt.plot(var, stability_cyclenumber_noAS, label="noAS")
        plt.plot(var, stability_cyclenumber_AS, label="AS")
        plt.legend()
        plt.show()

def perform_SVA_multiprocessing_forSao(
                var_list=None,
                var_path=None,
                save_fig=False,
                plot=True,
                parameters_input_path="../input_values/parameters_input.json",
                n_jobs=-1): 
    
    var, var_name, var_unit = var_list
    var_path1, var_path2 = var_path

    # --- BLOCK 1: noAS ---
    with open(parameters_input_path, "r") as f:
        parameters_input = json.load(f)

    print(f"Running noAS simulations on {n_jobs} cores...")
    results_noAS = Parallel(n_jobs=n_jobs)(
        delayed(run_single_simulation)(parameters_input, var_path1, var_path2, j) 
        for j in var
    )
    
    list_of_results_noAS, stability_cyclenumber_noAS = zip(*results_noAS)
    
    save_to_csv_SVA(list(list_of_results_noAS), var_list, f"SVA_data/{var_name}_noAS.csv")


    
    if plot:
        SVA_plot_from_csv_forSao(path_noAS=f"SVA_data/{var_name}_noAS.csv",
                        save_fig=save_fig,
                        complete_plot=True)

        plt.plot(var, stability_cyclenumber_noAS, label="noAS")
        plt.legend()
        plt.show()