import numpy as np
import pandas as pd
from scipy.stats import qmc
from run_simulation import main_fc
import json



# ================ Population set and analysis functions =========================================
def create_population(Npatients: int, file_name_output: str, poptype: str) -> None:
    """
    creates population noAS/AS
    """

    if poptype == "noAS":
        param_bounds = {
        "Emax_lv": [2.0, 4.0],
        "Emin_lv": [0.01, 0.1],
        "Amax_av": [3.0, 4.0],
        "HR": [60, 90],
        "Ra_syst": [0.8, 1.5],
        "Cv_syst": [50.0, 90.0],
        "Rv_syst": [0.05, 0.1],
        "Rmicro_max_LAD": [50.0, 70.0],
        "Rmicro_max_LCX": [60.0, 80.0],
        "icPv_syst": [2.0, 20.0]
        }
    elif poptype == "AS":
        param_bounds = {
        "Emax_lv": [2.0, 6.0],
        "Emin_lv": [0.1, 0.5],
        "Amax_av": [0.2, 1.0],
        "HR": [60, 90],
        "Ra_syst": [0.8, 1.5],
        "Cv_syst": [50.0, 90.0],
        "Rv_syst": [0.05, 0.1],
        "Rmicro_max_LAD": [50.0, 70.0],
        "Rmicro_max_LCX": [60.0, 80.0],
        "icPv_syst": [2.0, 20.0]
        }

    param_names = list(param_bounds.keys())
    param_names = list(param_bounds.keys())
    l_bounds = [param_bounds[key][0] for key in param_names]
    u_bounds = [param_bounds[key][1] for key in param_names]
    num_params = len(param_names)

    sampler = qmc.LatinHypercube(d=num_params, seed=42)
    sample_unit_hypercube = sampler.random(n=Npatients)

    patient_data_scaled = qmc.scale(sample_unit_hypercube, l_bounds, u_bounds)
    df = pd.DataFrame(patient_data_scaled, columns=param_names)
    df.insert(0, 'Patient_ID', range(1, Npatients + 1))

    df.to_csv(f"{file_name_output}.csv", index=False)

    print(f"Created  population {poptype}, saved as {file_name_output}.csv", flush=True)




def select_healthy(file_name_input_output: str, file_name_input_input: str, file_name_output_input: str, file_name_output_output: str) -> None:
    """
    selects healthy individuals
    """
    df = pd.read_csv(f"{file_name_input_output}.csv")

    selected = df[
        # stability
            (df["stability_reached"] == 1) &
        # energetics
            (df["DeltaE"] > 0) &
        # pressures in circulation
            (df["Pa_syst_min"] > 70) &
            (df["Pa_syst_max"] < 160) &
            (df["Ppt_mean"]<=20) &
        # right ventricle
            (df["Prv_max"]<30)&
            (df["EDP_rv"]<8)&
            (df["Vrv_max"]<220)&
        # left ventricle
            (df["EDP_lv"]<15)&
            (df["Vlv_max"]<200)&
            ((df["Vav"]/df["Vlv_max"])>0.4)
            ]
    
    selected.to_csv(f"{file_name_output_output}.csv", index=False)

    print(f"Selected output healthy subpopulation from {file_name_input_output}.csv, saved as {file_name_output_output}.csv", flush=True)

    selected_indexes = selected["Patient_ID"].tolist()

    df_in = pd.read_csv(f"{file_name_input_input}.csv")

    selected_in = df_in[df_in["Patient_ID"].isin(selected_indexes)]

    selected_in.to_csv(f"{file_name_output_input}.csv", index=False)

    print(f"Selected input healthy subpopulation from {file_name_input_input}.csv, saved as {file_name_output_input}.csv", flush=True)


def set_simulation_parameters(sim_type: int, file_name_input: str, file_name_output: str) -> None:
    """
    sets simulation conditions:
    1) isolated reduction in Rasyst
    2) reduction in Rasyst + tachycardia + increase in Emax + coronary dilatation
    3) reduction in Rasyst + tachycardia + decrease in Emax + coronary dilatation
    4) reduction in Rasyst + bradycardia + decrease in Emax + coronary dilatation
    """

    df = pd.read_csv(f"{file_name_input}.csv")

    num_patients = df.shape[0]

    if sim_type == 1:

        param_bounds = {
        "Ra_syst": [0.1, 0.5]
        }


    elif sim_type == 2:
        param_bounds = {
        "Ra_syst": [0.1, 0.5],
        "Emax_lv": [4.0, 6.0],
        "HR": [90.0, 160.0],
        "Rmicro_max_LAD": [12.5, 50.0],
        "Rmicro_max_LCX": [15.0, 60.0]
        }

    elif sim_type == 3:
        param_bounds = {
        "Ra_syst": [0.1, 0.5],
        "Emax_lv": [1.0, 2.0],
        "HR": [90.0, 160.0],
        "Rmicro_max_LAD": [12.5, 50.0],
        "Rmicro_max_LCX": [15.0, 60.0]
        }

    elif sim_type == 4:
        param_bounds = {
        "Ra_syst": [0.1, 0.5],
        "Emax_lv": [1.0, 2.0],
        "HR": [40.0, 60.0],
        "Rmicro_max_LAD": [12.5, 50.0],
        "Rmicro_max_LCX": [15.0, 60.0]
        }

    param_names = list(param_bounds.keys())
    l_bounds = [param_bounds[key][0] for key in param_names]
    u_bounds = [param_bounds[key][1] for key in param_names]
    num_params = len(param_names)

    sampler = qmc.LatinHypercube(d=num_params, seed=42)
    sample_unit_hypercube = sampler.random(n=num_patients)

    patient_data_scaled = qmc.scale(sample_unit_hypercube, l_bounds, u_bounds)
    df_to_insert = pd.DataFrame(patient_data_scaled, columns=param_names)

    for key in param_names:
        df[key] = df_to_insert[key]

    df.to_csv(f"{file_name_output}.csv", index=False)

    print(f"Created simulation type {sim_type} dataset for {file_name_input}, saved as {file_name_output}.csv", flush=True)


if __name__ == "__main__":
    # ===== noAS =====
    # generate universal input for other parameters
    parameters_input_to_save = {
        "chambers": {
            "Emax_ra": 0.25,
            "Emin_ra": 0.04,
            "m1_ra": 1.32,
            "m2_ra": 13.1,
            "V0ra": 17,
            "Emax_rv": 0.6,
            "Emin_rv": 0.04,
            "m1_rv": 1.32,
            "m2_rv": 27.4,
            "V0rv": 55,
            "Emax_la": 0.17,
            "Emin_la": 0.08,
            "m1_la": 1.32,
            "m2_la": 13.1,
            "V0la": 3,
            "Emax_lv": 3.0,
            "Emin_lv": 0.08,
            "m1_lv": 1.32,
            "m2_lv": 27.4,
            "V0lv": 10
        },
        "valves": {
            "Amax_tv": 5.0,
            "Amin_tv": 1e-6,
            "Kvo_tv": 0.03,
            "Kvc_tv": 0.03,
            "Amax_pv": 5.0,
            "Amin_pv": 1e-6,
            "Kvo_pv": 0.02,
            "Kvc_pv": 0.02,
            "Amax_mv": 5.0,
            "Amin_mv": 1e-6,
            "Kvo_mv": 0.03,
            "Kvc_mv": 0.03,
            "Amax_av": 4.0,
            "Amin_av": 1e-6,
            "Kvo_av": 0.012,
            "Kvc_av": 0.012,
            
        },
        "circulation": {
            "HR": 70.0,
            "density": 1.06,
            "Rao": 0.0398,
            "Ca_syst": 1.33,
            "Ra_syst": 1.2, 
            "Cv_syst": 70.0,
            "Rv_syst": 0.075,
            "Rpa": 0.01,
            "Ca_pulm": 3.55,
            "Ra_pulm": 0.1,
            "Cv_pulm": 20.5,
            "Rv_pulm": 0.006
        },
        "CBF_model": {
            "Ra_LAD": 0.7299,
            "Ra_LCX": 0.8759,
            "Ra_RCA": 0.9290,
            "Rmicro_max_LAD": 61.3140,
            "Rmicro_min_LAD": 61.3140*0.7,
            "Rmicro_max_LCX": 73.6768,
            "Rmicro_min_LCX": 73.6768*0.7,
            "Rmicro_max_RCA": 78.0360,
            "Rmicro_min_RCA": 78.0360*0.75,
            "Rv_LAD": 1.09489,
            "Rv_LCX": 1.31387,
            "Rv_RCA": 1.39350,
            "Ca_LAD": 0.01698,
            "Ca_LCX": 0.01698,
            "Ca_RCA": 0.01698,
            "Cmicro_LAD": 0.66215,
            "Cmicro_LCX": 0.66215,
            "Cmicro_RCA": 0.66215,
            "hypertrophy_factor": 1.0
        },
        "initial_conditions": {
            "icVra": 40,
            "icVrv": 180,
            "icVla": 27,
            "icVlv": 135,
            "icQtv": 0.0,
            "icQpv": 0.0,
            "icQmv": 0.0,
            "icQav": 0.0,
            "iceps_tv": 0.0,
            "iceps_pv": 0.0,
            "iceps_mv": 0.0,
            "iceps_av": 0.0,
            "icPa_syst": 60.0,
            "icPv_syst": 14.0,
            "icPa_pulm": 15.0,
            "icPv_pulm": 10.0,
            "icPa_LAD": 50.0,
            "icPv_LAD": 14.0,
            "icPa_LCX": 50.0,
            "icPv_LCX": 14.0,
            "icPa_RCA": 50.0,
            "icPv_RCA": 14.0
        }
    }

    with open("universal_input_noAS.json", "w") as f:
        json.dump(parameters_input_to_save, f)


    # generate dataset for noAS
    create_population(Npatients=200000, file_name_output="input_allpopulation_noAS", poptype="noAS")

    # computation for no AS + selection
    main_fc(file_name_universal_input="universal_input_noAS", file_name_specific_input="input_allpopulation_noAS", file_name_output="output_allpopulation_noAS")
    select_healthy(file_name_input_output="output_allpopulation_noAS",
                   file_name_input_input="input_allpopulation_noAS",
                   file_name_output_output="output_selpopulation_noAS",
                   file_name_output_input="input_selpopulation_noAS")

    # create and run simulations for noAS
    set_simulation_parameters(sim_type=1, file_name_input="input_selpopulation_noAS", file_name_output="input_selpopulation_sim1_noAS")
    main_fc(file_name_universal_input="universal_input_noAS", file_name_specific_input="input_selpopulation_sim1_noAS", file_name_output="output_selpopulation_sim1_noAS")

    set_simulation_parameters(sim_type=2, file_name_input="input_selpopulation_noAS", file_name_output="input_selpopulation_sim2_noAS")
    main_fc(file_name_universal_input="universal_input_noAS", file_name_specific_input="input_selpopulation_sim2_noAS", file_name_output="output_selpopulation_sim2_noAS")

    set_simulation_parameters(sim_type=3, file_name_input="input_selpopulation_noAS", file_name_output="input_selpopulation_sim3_noAS")
    main_fc(file_name_universal_input="universal_input_noAS", file_name_specific_input="input_selpopulation_sim3_noAS", file_name_output="output_selpopulation_sim3_noAS")

    set_simulation_parameters(sim_type=4, file_name_input="input_selpopulation_noAS", file_name_output="input_selpopulation_sim4_noAS")
    main_fc(file_name_universal_input="universal_input_noAS", file_name_specific_input="input_selpopulation_sim4_noAS", file_name_output="output_selpopulation_sim4_noAS")



