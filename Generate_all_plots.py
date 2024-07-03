import ROOT as ROOT
import os

from configurations import (
    file_names,
    low_voltage,
    BV,
    channels,
    configurations,
)
from Analysis_functions import (
    Draw_amp,
    Get_deltaT_vs_amplitud_2D_map,
    Get_distribution_centers,
    Do_polonomial_fit,
    Get_deltaT_vs_amplitud_Corrected,
    Get_time_resolution,
    Draw_bias_scan,
)

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

for channel in channels:

    print("Working in channel: ", channel)
    sigma_config, error_config = {}, {}

    for i, configs in enumerate(configurations):

        print("Running over configs: ", configs)

        file_path = "./data/"
        plot_path = (
            "./plots/"
            + BV[configs]
            + "V/LV"
            + low_voltage[configs].replace(".", "p")
            + "/"
        )
        os.makedirs(plot_path, exist_ok=True)

        Draw_amp(
            file_path,
            configs,
            channel,
            "ch%i_amplitude_channel" % (channel),
            plot_path,
        )

        deltaT_vs_amplitud = Get_deltaT_vs_amplitud_2D_map(
            file_path,
            configs,
            channel,
            "ch%i_deltaT_vs_amplitud" % (channel),
            plot_path,
        )

        time_walk = Get_distribution_centers(
            deltaT_vs_amplitud,
            channel,
            "ch%i_deltaT_vs_amplitud_TimeWalk" % (channel),
            plot_path,
        )
        time_walk_fit = Do_polonomial_fit(
            deltaT_vs_amplitud,
            time_walk,
            configs,
            channel,
            "4",  # degree of the polynomial
            "ch%i_deltaT_vs_amplitud_TimeWalk_fitted" % (channel),
            plot_path,
        )

        deltaT_vs_amplitud_corrected = Get_deltaT_vs_amplitud_Corrected(
            file_path,
            time_walk_fit,
            configs,
            channel,
            "ch%i_deltaT_vs_amplitud_TimeWalk_fitted_corrected" % (channel),
            plot_path,
        )

        sigma, sigma_error = Get_time_resolution(
            deltaT_vs_amplitud,
            configs,
            channel,
            "ch%i_Time_resolution" % (channel),
            plot_path,
        )
        sigma_corrected, sigma_error_corrected = Get_time_resolution(
            deltaT_vs_amplitud_corrected,
            configs,
            channel,
            "ch%i_Time_resolution_corrected" % (channel),
            plot_path,
        )

        if channel == 3:
            sigma_config[configs[0]] = sigma_corrected
            error_config[configs[0]] = sigma_error_corrected
        else:
            sigma_config[configs[0]] = sigma
            error_config[configs[0]] = sigma_error

    print("Ready with all the configs")
    # Do low voltage Bias S5can
    for bias in set(BV.values()):

        plot_path = "./plots/" + bias + "V/"
        lv = []
        sigma = []
        error = []
        # for configs, lv in low_voltage.items():
        for configs, bv in BV.items():
            if bias == bv:
                lv.append(float(low_voltage[configs[0]]))
                sigma.append(sigma_config[configs[0]])
                error.append(error_config[configs[0]])

        Draw_bias_scan(
            lv,
            sigma,
            error,
            channel,
            "ch%i_board_low_voltage_scan" % (channel),
            plot_path,
            use_low_voltage=True,
        )

    # Do High voltage Bias Scan
    for voltage in set(low_voltage.values()):

        plot_path = "./plots/"
        bias = []
        sigma = []
        error = []
        # for configs, lv in low_voltage.items():
        for configs, lv in low_voltage.items():
            if voltage == lv:
                bias.append(float(BV[configs][0]))
                sigma.append(sigma_config[configs][0])
                error.append(error_config[configs][0])

        Draw_bias_scan(
            bias,
            sigma,
            error,
            channel,
            "ch%i_bias_voltage_scan_lv_" % (channel) + voltage.replace(".", "p"),
            plot_path,
            use_low_voltage=False,
        )

    print("Ready with scan plots")
