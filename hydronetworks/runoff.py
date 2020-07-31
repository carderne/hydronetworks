"""Functions for mode advances runoff calculations.
NOT CURRENTLY WORKING PROPERLY.
"""

import pandas as pd

SEC_P_YEAR = 3600 * 8760
SEC_P_MONTH = SEC_P_YEAR / 12


def rainfall_runoff(
    nodes_df,
    kc_path,
    precip_effectiveness=0.5,
    runoff_to_gw_fraction=0.5,
    runoff_calibration_accuracy=0.2,
):
    """
    Calculate the runoff from the rainfall using the rainfall-runoff method
    Calibrate against gscd annual values
    """

    runoff_df = pd.DataFrame(
        columns=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], index=nodes_df.index
    )
    k_c = pd.read_csv(kc_path, index_col=0)

    # formula from http://www.weap21.org/webhelp/hydrology.htm
    # for each point, calculate the runoff for each month
    # then compare to gscd annual value and modify params
    nodes_df["precip_effective"] = precip_effectiveness  # default starting values
    nodes_df["runoff_to_gw_fraction"] = runoff_to_gw_fraction
    for index, row in runoff_df.iterrows():
        print("calibrate {}".format(index))
        counter = 0
        while True:
            for col in runoff_df:
                et_ref = 1

                precip = nodes_df.loc[index, "{}{}".format("precip", col)] * SEC_P_MONTH
                precip_avail_for_et = precip * nodes_df.loc[index, "precip_effective"]
                et_potential = et_ref * k_c.loc[col, nodes_df.loc[index, "land_type"]]
                runoff = max(0, precip_avail_for_et - et_potential) + precip * (
                    1 - nodes_df.loc[index, "precip_effective"]
                )
                runoff_to_surface = runoff * (
                    1 - nodes_df.loc[index, "runoff_to_gw_fraction"]
                )

                runoff_df.loc[index, col] = runoff_to_surface

            ref_value = nodes_df.loc[index, "gscd"]
            calc_value = runoff_df.loc[index].mean() * 12  # to compare with annual gscd

            if counter > 10:
                break
            elif (
                abs((calc_value - ref_value) / ref_value) < runoff_calibration_accuracy
            ):
                break
            else:
                counter += 1
                nodes_df.loc[index, "precip_effective"] *= sorted(
                    [0.5, 1 + (calc_value - ref_value) / ref_value / 10, 1.5]
                )[1]
                nodes_df.loc[index, "runoff_to_gw_fraction"] *= sorted(
                    [0.5, 1 + (calc_value - ref_value) / ref_value / 10, 1.5]
                )[1]

    return nodes_df, runoff_df


def discharge(
    nodes_df,
    network,
    network_df,
    flowacc,
    runoff_df,
    water_loss_factor=0.00001,
    max_water_loss_fraction=0.2,
):
    area = flowacc.res[0] * flowacc.res[1]

    # calculate self discharge for each node
    # Q [m3/s] =
    #   runoff [mm/year] * flowacc [number] * area [m2] * 1e-3 [mm/m]
    #   / (8760*3600 [s/year])
    discharge_df = pd.DataFrame(columns=runoff_df.columns, index=runoff_df.index)
    for index, row in discharge_df.iterrows():
        print("discharge {}".format(index))
        for col in discharge_df:
            ru = runoff_df.loc[index, col]
            fa = nodes_df.loc[index, "flow_acc"]
            discharge_df.loc[index, col] = ru * fa * area * 1e-3 / SEC_P_YEAR

    nodes_df["discharge_local"] = (
        nodes_df["gscd"] * nodes_df["flow_acc_local"] * area * 1e-3 / SEC_P_YEAR
    )
    nodes_df["discharge_accum"] = nodes_df["discharge_local"]

    # start from Shreve order 1, contribute all discharges from upstream to downstream
    for so in range(1, max(network.T[7]) + 1):
        for arc in network:
            if arc[7] == so:
                print("contribute {}".format(arc[0]))
                nodes_df.loc[arc[6], "discharge_accum"] += nodes_df.loc[
                    arc[5], "discharge_accum"
                ] * (1 - max(max_water_loss_fraction, water_loss_factor * arc[8]))
                for col in discharge_df:
                    discharge_df.loc[arc[6], col] += discharge_df.loc[arc[5], col] * (
                        1 - max(max_water_loss_fraction, water_loss_factor * arc[8])
                    )

    nodes_df["discharge_max"] = -99
    nodes_df["discharge_mean"] = -99
    nodes_df["discharge_min"] = -99
    for index, node in nodes_df.iterrows():
        nodes_df.loc[index, "discharge_max"] = discharge_df.loc[index].max()
        nodes_df.loc[index, "discharge_mean"] = discharge_df.loc[index].mean()
        nodes_df.loc[index, "discharge_min"] = discharge_df.loc[index].min()

    # add the discharge from the upstream node of each arc to that arc
    network_df["discharge_accum"] = -99
    network_df["discharge_max"] = -99
    network_df["discharge_mean"] = -99
    network_df["discharge_min"] = -99
    for index, arc in network_df.iterrows():
        network_df.loc[index, "discharge_accum"] = nodes_df["discharge_accum"][
            arc["nif"]
        ]
        network_df.loc[index, "discharge_max"] = nodes_df["discharge_max"][arc["nif"]]
        network_df.loc[index, "discharge_mean"] = nodes_df["discharge_mean"][arc["nif"]]
        network_df.loc[index, "discharge_min"] = nodes_df["discharge_min"][arc["nif"]]
