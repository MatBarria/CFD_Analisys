import numpy as np
import ROOT as ROOT
import uproot
import os
from configurations import (
    BV,
    amplitude_region,
    x_cuts,
    y_cuts,
    deltaT_binning,
    file_names,
)

from histo_utilities import (
    create_TH1D,
    create_TH2D,
    create_TGraph,
)


def Get_tracker_cuts(branches):

    amplitudes = branches["amp"].array()
    return (branches["nplanes"].array() > 8) & (amplitudes[:, 7] > 100)


def Get_fiducial_cuts(branches, config, channel):

    ch_amp = "ch" + str(channel)
    x = branches["x_dut"].array()[:, 5]
    y = branches["y_dut"].array()[:, 5]

    fiducial_cuts = (
        (x > x_cuts[config][ch_amp][0])
        & (x < x_cuts[config][ch_amp][1])
        & (y > y_cuts[config][ch_amp][0])
        & (y < y_cuts[config][ch_amp][1])
    )
    return fiducial_cuts


def Get_final_branches(file_path, cofigurations, channel):

    channel_time = np.array([])
    tracker_time = np.array([])
    channel_amplitude = np.array([])

    for config in cofigurations:

        branches = uproot.open(file_path + file_names[config])["pulse"]

        tracker_cuts = Get_tracker_cuts(branches)
        fiducial_cuts = Get_fiducial_cuts(branches, config, channel)
        all_cuts = (tracker_cuts) & (fiducial_cuts)

        amplitudes = branches["amp"].array()
        LP2_50 = branches["LP2_50"].array()
        tracker_time = np.concatenate((tracker_time, LP2_50[:, 7][all_cuts]))
        channel_time = np.concatenate((channel_time, LP2_50[:, channel - 3][all_cuts]))
        channel_amplitude = np.concatenate(
            (channel_amplitude, amplitudes[:, channel][all_cuts])
        )

    return channel_amplitude, channel_time, tracker_time


def Draw_amp(file_path, configurations, channel, plot_name, plot_path):

    amplitude_channel, _, _ = Get_final_branches(file_path, configurations, channel)

    ch_amp = "ch" + str(channel)
    title = ["Amplitude [mV]", "Counts"]

    amplitude_hist = create_TH1D(
        amplitude_channel,
        title="amplitude " + ch_amp,
        axis_title=[title[0], title[1]],
        binning=[100, 0, 500],
    )

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)

    amplitude_hist.Draw("hist e")
    canvas.Draw()
    canvas.SaveAs(plot_path + plot_name + ".png")


# TODO: Finish WIP
# def Get_efficiency_2D_map(file_name, file_path, config, channel, plot_name, plot_path, bias_voltage, draw=False):

#     title = ['x [mm]', 'y [mm]']

#     branches = uproot.open(file_path + file_name)['pulse']
#     amplitudes = branches['amp'].array()
#     LP2_50 = branches['LP2_50'].array()
#     x = branches['x_dut'].array()[:, 5]
#     y = branches['y_dut'].array()[:, 5]

#     ch_amp = "ch" + str(channel)
#     amplitud_min = amplitude_region[config][ch_amp][0]
#     amplitud_max = amplitude_region[config][ch_amp][1]
#     if channel == 3:
#         bin_width_amplitud = 5
#     else:
#         bin_width_amplitud = 15
#     bin_number_amplitud = int(
#         (amplitud_max - amplitud_min) / bin_width_amplitud)

#     bins = [bin_number_amplitud, amplitud_min, amplitud_max, deltaT_binning[config][ch_amp][0],
#             deltaT_binning[config][ch_amp][1], deltaT_binning[config][ch_amp][2]]

#     tracker_cuts = (branches['nplanes'].array() > 8) & (amplitudes[:, 7] > 100)
#     # tracker_cuts = (branches['nplanes'].array() > 8) & (branches['npix'].array() > 3) &(amplitudes[:, 7] > 100)
#     fiducial_cuts = (x > x_cuts[config][ch_amp][0]) & (x < x_cuts[config][ch_amp][1]) & (
#         y > y_cuts[config][ch_amp][0]) & (y < y_cuts[config][ch_amp][1])
#     # amplitude_cuts = np.greater_equal(.5*amplitudes[:,3] , amplitudes[:,2])
#     # amplitude_cuts = .5*amplitudes[:,3] > amplitudes[:,2]

#     all_cuts = (tracker_cuts) & (fiducial_cuts)
#     deltaT = (LP2_50[:, channel - 3] - LP2_50[:, 7])*1e9

#     hist_name = "Amp vs Delta T " + ch_amp
#     deltaT_vs_amplitud = create_TH2D(np.column_stack(
#         (amplitudes[:, channel][all_cuts], deltaT[all_cuts])),
#         hist_name, axis_title=[title[0], title[1], 'counts'],
#         binning=bins)

#     if draw:
#         map_canvas = ROOT.TCanvas("2DMap", "2DMap", 800, 600)
#         deltaT_vs_amplitud.Draw("COLZ")
#         map_canvas.SaveAs(plot_path + plot_name + ".png")
#     return deltaT_vs_amplitud


def Get_deltaT_vs_amplitud_2D_map(
    file_path,
    configurations,
    channel,
    plot_name,
    plot_path,
    draw=False,
):

    title = ["Amplitude [mV]", "Delta T [ns]"]

    channel_amplitude, channel_time, tracker_time = Get_final_branches(
        file_path, configurations, channel
    )

    ch_amp = "ch" + str(channel)
    amplitud_min = amplitude_region[configurations[0]][ch_amp][0]
    amplitud_max = amplitude_region[configurations[0]][ch_amp][1]
    if channel == 3:
        bin_width_amplitud = 5
    else:
        bin_width_amplitud = 15
    bin_number_amplitud = int((amplitud_max - amplitud_min) / bin_width_amplitud)

    bins = [
        bin_number_amplitud,
        amplitud_min,
        amplitud_max,
        deltaT_binning[ch_amp][0],
        deltaT_binning[ch_amp][1],
        deltaT_binning[ch_amp][2],
    ]

    deltaT = (channel_time - tracker_time) * 1e9
    hist_name = "Amp vs Delta T " + ch_amp
    deltaT_vs_amplitud = create_TH2D(
        np.column_stack((channel_amplitude, deltaT)),
        hist_name,
        axis_title=[title[0], title[1], "counts"],
        binning=bins,
    )

    if draw:
        map_canvas = ROOT.TCanvas("2DMap", "2DMap", 800, 600)
        deltaT_vs_amplitud.Draw("COLZ")
        map_canvas.SaveAs(plot_path + plot_name + ".png")
    return deltaT_vs_amplitud


def Get_distribution_centers(hist, channel, plot_name, plot_path, draw=True):

    os.makedirs(plot_path + "fits/", exist_ok=True)
    bins_number = hist.GetXaxis().GetNbins()

    time_walk = ROOT.TGraph()
    for bin in range(1, bins_number + 1):

        tmp_hist = hist.ProjectionY("py", bin, bin)

        rms = tmp_hist.GetRMS()
        mean = tmp_hist.GetMean()

        fit_low = mean - 1.5 * rms
        fit_high = mean + 1.5 * rms

        fit = ROOT.TF1("fit", "gaus", fit_low, fit_high)
        tmp_hist.Fit(fit, "Q", "", fit_low, fit_high)

        tmp_canvas = ROOT.TCanvas("canvas" + str(bin), "canvas" + str(bin))
        tmp_canvas.cd()
        tmp_hist.Draw("hist")
        fit.Draw("same")
        if draw:
            tmp_canvas.SaveAs(
                plot_path
                + "fits/fit_bin_channel_"
                + str(channel)
                + "_"
                + str(bin)
                + ".png"
            )

        mean_y = hist.GetMean(2)
        if fit.GetParameter(1) > mean_y + 0.1 or fit.GetParameter(1) < mean_y - 0.1:
            time_walk.AddPoint(hist.GetXaxis().GetBinCenter(bin), mean_y)
        else:
            time_walk.AddPoint(hist.GetXaxis().GetBinCenter(bin), fit.GetParameter(1))
    draw = False
    if draw:
        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        hist.Draw("COLZ")
        time_walk.SetMarkerStyle(8)
        time_walk.SetMarkerSize(0.8)
        time_walk.SetMarkerColor(ROOT.kBlack)
        time_walk.Draw("P Same")
        canvas.SaveAs(plot_path + plot_name + ".png")
    return time_walk


def Do_polonomial_fit(
    hist,
    graph,
    configurations,
    channel,
    number_of_paremeters,
    plot_name,
    plot_path,
    draw=True,
):

    ch_amp = "ch" + str(channel)
    amplitud_min = amplitude_region[configurations[0]][ch_amp][0]
    amplitud_max = amplitude_region[configurations[0]][ch_amp][1]

    fit = ROOT.TF1("fit", "pol " + number_of_paremeters, amplitud_min, amplitud_max)

    # graph.Fit(fit, "R", "", amplitud_min, amplitud_max);
    graph.Fit(fit, "QR", "", amplitud_min, amplitud_max)
    # for i in range(0, fit.GetNpar()):
    # print("Fit parameter %i: %.5f"%(i, fit.GetParameter(i)))

    graph.SetMarkerStyle(8)
    graph.SetMarkerSize(0.8)
    graph.SetMarkerColor(ROOT.kBlack)
    fit.SetLineColor(ROOT.kRed)

    if draw:
        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        hist.Draw("COLZ")
        graph.Draw("P Same")
        fit.Draw("Same")
        canvas.SaveAs(plot_path + plot_name + ".png")

    return fit


def Get_deltaT_vs_amplitud_Corrected(
    file_path, fit, configurations, channel, plot_name, plot_path
):

    title = ["Amplitude [mV]", "Delta T[ns]"]

    channel_amplitude, channel_time, tracker_time = Get_final_branches(
        file_path, configurations, channel
    )

    ch_amp = "ch" + str(channel)
    amplitud_min = amplitude_region[configurations[0]][ch_amp][0]
    amplitud_max = amplitude_region[configurations[0]][ch_amp][1]
    if channel == 3:
        bin_width_amplitud = 5
    else:
        bin_width_amplitud = 15
    bin_number_amplitud = int((amplitud_max - amplitud_min) / bin_width_amplitud)
    deltaT_range = deltaT_binning[ch_amp][2] - deltaT_binning[ch_amp][1]

    bins = [
        bin_number_amplitud,
        amplitud_min,
        amplitud_max,
        deltaT_binning[ch_amp][0],
        -deltaT_range / 2,
        deltaT_range / 2,
    ]

    deltaT = (channel_time - tracker_time) * 1e9

    # TODO: I know this can be automatized but idk how to get a power of an array...
    number_of_paremeters = fit.GetNpar()
    corrected_deltaT = np.array(deltaT)
    for i in range(0, number_of_paremeters):
        # print(amplitudes[:, channel].array())
        corrected_deltaT -= fit.GetParameter(i) * np.power(
            np.array(channel_amplitude), i
        )

    deltaT_vs_amplitud_corrected = create_TH2D(
        np.column_stack((channel_amplitude, corrected_deltaT)),
        "Corrected Amp vs Delta T " + ch_amp,
        axis_title=[title[0], title[1], "counts"],
        binning=bins,
    )

    # map_canvas = ROOT.TCanvas("2DMap", "2DMap", 800, 600)
    # deltaT_vs_amplitud.Draw("COLZ")
    # map_canvas.SaveAs(plot_path + plot_name + ".png")
    time_walk_corrected = Get_distribution_centers(
        deltaT_vs_amplitud_corrected,
        channel,
        "ch%i_deltaT_vs_amplitud_TimeWalk_corrected" % (channel),
        plot_path,
        False,
    )

    Do_polonomial_fit(
        deltaT_vs_amplitud_corrected,
        time_walk_corrected,
        configurations,
        channel,
        "2",
        plot_name,
        plot_path,
    )

    return deltaT_vs_amplitud_corrected


def Get_time_resolution(hist, configurations, channel, plot_name, plot_path):

    time_difference = hist.ProjectionY("Time Res")
    htitle = "Delta T ch" + str(channel)
    time_difference.SetTitle(htitle)
    time_difference.GetYaxis().SetTitle("Counts")

    rms = time_difference.GetRMS()
    mean = time_difference.GetMean()
    fit_low = mean - 1.5 * rms
    fit_high = mean + 1.5 * rms
    fit = ROOT.TF1("fit_ch" + str(channel), "gaus", fit_low, fit_high)
    time_difference.Fit(fit, "Q", "", fit_low, fit_high)

    # print("For Channel: ", channel)
    # print("Sigma from the fit : %.3f +- %.3f [ns]" %
    # (fit.GetParameter(2), fit.GetParError(2)))

    sigma = fit.GetParameter(2)
    sigma_error = fit.GetParError(2)

    event_text = ROOT.TText(0.55, 0.64, "Events: " + str(time_difference.Integral()))
    event_text.SetNDC()
    # event_text.SetTextAlign(22)
    sigma_text = ROOT.TText(
        0.55, 0.6, "Sigma: %.3f +- %.3f [ns]" % (sigma, sigma_error)
    )
    sigma_text.SetNDC()
    # event_text.SetTextAlign(22)
    channel_text = ROOT.TText(
        0.55, 0.7, "channel: %i BV:" % (channel) + BV[configurations[0]] + "V"
    )
    channel_text.SetNDC()

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    time_difference.Draw("hist e")
    fit.Draw("Same")
    event_text.Draw("Same")
    sigma_text.Draw("Same")
    channel_text.Draw("Same")
    canvas.Draw()
    canvas.SaveAs(plot_path + plot_name + ".png")

    return sigma, sigma_error


def Draw_bias_scan(
    voltages, sigma, error, channel, plot_name, plot_path, use_low_voltage=False
):

    xtitle = "Bias Current [V]" if not use_low_voltage else "Board low voltage [V]"
    htitle = "Bias scan" if not use_low_voltage else "Board low voltage scan"
    if channel == 3:
        htitle += " Dt Corrected"

    graph = create_TGraph(
        voltages,
        sigma,
        ex=np.zeros(len(voltages)),
        ey=error,
        axis_title=[xtitle, "Time resolution [ns]"],
    )
    graph.SetTitle(htitle)

    graph.SetMarkerStyle(8)
    graph.SetMarkerSize(0.8)
    graph.SetMarkerColor(ROOT.kBlack)

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
    graph.Draw("ALP")
    canvas.SaveAs(plot_path + plot_name + ".png")
