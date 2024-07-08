import ROOT
from ROOT import gROOT, kError
import os

from configurations import (
    low_voltage,
    BV,
    channels,
    configurations,
)
from Analysis_functions import (
    draw_amp,
    get_deltaT_vs_amplitud_2D_map,
    get_distribution_centers,
    do_polonomial_fit,
    get_deltaT_vs_amplitud_corrected,
    get_time_resolution,
    draw_bias_scan,
)

gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = kError

for channel in channels:
    print("Working in channel: ", channel)
    sigma_config, error_config = {}, {}

    for i, configs in enumerate(configurations):
        print("Running over configs: ", configs)

        file_path = "./data/"
        plot_path = (
            "./plots/"
            + BV[configs[0]]
            + "V/LV"
            + low_voltage[configs[0]].replace(".", "p")
            + "/"
        )
        os.makedirs(plot_path, exist_ok=True)

        draw_amp(
            file_path,
            configs,
            channel,
            "ch%i_amplitude_channel" % (channel),
            plot_path,
        )

        deltaT_vs_amplitud = get_deltaT_vs_amplitud_2D_map(
            file_path,
            configs,
            channel,
            "ch%i_deltaT_vs_amplitud" % (channel),
            plot_path,
        )

        time_walk = get_distribution_centers(
            deltaT_vs_amplitud,
            channel,
            "ch%i_deltaT_vs_amplitud_TimeWalk" % (channel),
            plot_path,
        )
        time_walk_fit = do_polonomial_fit(
            deltaT_vs_amplitud,
            time_walk,
            configs,
            channel,
            "4",  # degree of the polynomial
            "ch%i_deltaT_vs_amplitud_TimeWalk_fitted" % (channel),
            plot_path,
        )

        deltaT_vs_amplitud_corrected = get_deltaT_vs_amplitud_corrected(
            file_path,
            time_walk_fit,
            configs,
            channel,
            "ch%i_deltaT_vs_amplitud_TimeWalk_fitted_corrected" % (channel),
            plot_path,
        )

        sigma, sigma_error = get_time_resolution(
            deltaT_vs_amplitud,
            configs,
            channel,
            "ch%i_Time_resolution" % (channel),
            plot_path,
        )
        sigma_corrected, sigma_error_corrected = get_time_resolution(
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
        # for configs, bv in BV.items():
        for configs in configurations:
            if bias == BV[configs[0]]:
                lv.append(float(low_voltage[configs[0]]))
                sigma.append(sigma_config[configs[0]])
                error.append(error_config[configs[0]])

        draw_bias_scan(
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
        for configs in configurations:
            if voltage == low_voltage[configs[0]]:
                bias.append(float(BV[configs[0]]))
                sigma.append(sigma_config[configs[0]])
                error.append(error_config[configs[0]])

        draw_bias_scan(
            bias,
            sigma,
            error,
            channel,
            "ch%i_bias_voltage_scan_lv_" % (channel) + voltage.replace(".", "p"),
            plot_path,
            use_low_voltage=False,
        )

    print("Ready with scan plots")
