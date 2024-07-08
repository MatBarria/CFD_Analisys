import os
from sys import argv
import ROOT
import uproot
import subprocess

# from configurations import configurations as configRuns

CFD_REPO = (
    subprocess.Popen(["git", "rev-parse", "--show-toplevel"], stdout=subprocess.PIPE)
    .communicate()[0]
    .rstrip()
    .decode("utf-8")
)

out_path = CFD_REPO + "/data/"
path = CFD_REPO
os.makedirs(out_path, exist_ok=True)


ROOT.gROOT.SetBatch(True)
configRuns = {
    "131": [715535, 715650],
    "134": [715679, 715689],
    "135": [715690, 715720],
    "136": [715721, 715752],
    "137": [715753, 715783],
    "138": [715784, 715826],
    "139": [715827, 715895],
    "140": [715896, 715978],
    "143": [716003, 716018],
    "147": [716031, 716127],
    "150": [716128, 716229],
    "152": [716280, 716352],
    "153": [716353, 716395],
    "154": [716396, 716436],
    "156": [716437, 716456],
}

if __name__ == "__main__":

    if len(argv) == 1:
        configurations = [str(x) for x in range(131, 160)]
    elif len(argv) == 2:
        configurations = [argv[1]]
    elif len(argv) == 3:
        configurations = [str(x) for x in range(int(argv[1]), int(argv[2]) + 1)]
    else:
        print("Invalid number of arguments")
        exit()

    for config in configurations:

        if config not in configRuns.keys():
            continue

        print("Running for config:", config)
        allRunsList = [
            str(x) for x in range(configRuns[config][0], configRuns[config][1] + 1)
        ]
        runNumList = []

        print("Total runs: ", len(allRunsList))

        for run in allRunsList:

            print("Checking run: ", run)
            path = "/eos/uscms/store/group/cmstestbeam/2024_05_SNSPD_FCFD_ETL//LecroyScope/RecoData/TimingDAQRECO/RecoWithTracks/v14/"
            # Check if the File exists
            if not os.path.isfile(path + "run" + run + "_info.root"):
                print("File does not exist: ", run)
                continue

            # Check the efficiency of the LGAD to remove bad runs

            file = ROOT.TFile(path + "run" + run + "_info.root", "READ")
            tuple = file.Get("pulse")
            if int(config) >= 131:
                try:
                    tuple.Draw(
                        "(amp[0]>100):y_dut[0]:x_dut[0]>>tmp(1,-1.1,2.8,1,3,3.8)",
                        "nplanes>8 && npix>3 && amp[7]>100",
                        "profcolz",
                    )
                    htemp = ROOT.gPad.GetPrimitive("tmp")
                    print(
                        "Eff is : "
                        + str(htemp.GetBinContent(1, 1))
                        + " for run: "
                        + run
                    )
                    if htemp.GetBinContent(1, 1) > -1:
                        runNumList.append(run)
                    else:
                        print("Bad Run: ", run)
                        continue
                except:
                    print("Bad run: " + str(run) + "(we can probably recover this one)")
                    continue

            branches = ["nplanes", "npix", "amp", "LP2_50", "x_dut", "y_dut"]

            filename = "run" + run + "_tmp.root"

            sel_br = uproot.open(path + "run" + run + "_info.root")["pulse"].arrays(
                branches, library="np"
            )
            events_number = len(sel_br["nplanes"])

            cut_expression = sel_br["nplanes"] > 8
            selected_events = {
                name: array[cut_expression] for name, array in sel_br.items()
            }

            with uproot.recreate(filename) as file:
                file["pulse"] = selected_events

        print("There are " + str(len(runNumList)) + " good runs")
        output = out_path + "run{}_{}_skimpy.root".format(
            configRuns[config][0], configRuns[config][1]
        )
        print("Attempting to processes the following runs:", runNumList)
        files = []
        cmd = "hadd -f {}".format(output)
        for runNum in runNumList:
            f = "./run" + runNum + "_tmp.root"
            if os.path.exists(f):
                cmd += " " + f
            else:
                print("WARNING: {} doesn't exist".format(f))
        print(cmd)
        os.system(cmd)
        os.system("rm *_tmp.root")
