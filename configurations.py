from collections import OrderedDict

configurations = [
    ["131"],  # 150V
    ["135"],  # 170V
    ["136"],  # 180V
    ["138"],  # 185V
    ["139"],  # 190V , lv = 1.06V
    ["137"],  # 190V , lv = 1.08V
    ["140"],  # 190V , lv = 1.12V
    ["147"],  # 196V
    ["143", "150"],  # 199V
]

file_names = {
    "131": "run715535_715650_skimpy.root",
    "134": "run715679_715689_skimpy.root",
    "135": "run715690_715720_skimpy.root",
    "136": "run715721_715752_skimpy.root",
    "137": "run715753_715783_skimpy.root",
    "138": "run715784_715826_skimpy.root",
    "139": "run715827_715895_skimpy.root",
    "140": "run715896_715978_skimpy.root",
    "143": "run716003_716018_skimpy.root",
    "147": "run716031_716127_skimpy.root",
    "150": "run716128_716229_skimpy.root",
}


# channels = [3]
channels = [3, 4, 5]

BV = OrderedDict()

# Its there is more than one cofiguration with the same BV write first the one
# with the lower low voltage
BV = {
    "131": "150",
    # "134": "080V"
    "135": "170",
    "136": "180",
    "138": "185",
    "139": "190",
    "137": "190",
    "140": "190",
    "147": "196",
    "143": "199",  # 199V few runs
    "150": "199",
}

low_voltage = OrderedDict()

# Its there is more than one cofiguration with the same low voltage write first the one
# with the lower BV
low_voltage = {
    "139": "1.06",
    "131": "1.08",
    # "134": "080V"
    "135": "1.08",
    "136": "1.08",
    "138": "1.08",
    "137": "1.08",
    "147": "1.08",
    "143": "1.08",  # 199V few runs
    "150": "1.08",
    "140": "1.12",
}

amplitude_region = {
    "131": {
        "ch3": [70, 110],
        "ch4": [145, 230],
        "ch5": [60, 200],
    },
    "134": {
        "ch3": [40, 80],
        "ch4": [150, 230],
        "ch5": [60, 200],
    },
    "135": {
        "ch3": [75, 130],
        "ch4": [180, 300],
        "ch5": [130, 300],
    },
    "136": {
        "ch3": [75, 140],
        "ch4": [200, 340],
        "ch5": [160, 360],
    },
    "138": {
        "ch3": [75, 150],
        "ch4": [220, 350],
        "ch5": [180, 350],
    },
    "137": {
        "ch3": [75, 170],
        "ch4": [230, 380],
        "ch5": [230, 400],
    },
    "139": {
        "ch3": [70, 200],
        "ch4": [230, 410],
        "ch5": [230, 410],
    },
    "140": {
        "ch3": [70, 190],
        "ch4": [240, 400],
        "ch5": [200, 450],
    },
    "143": {
        "ch3": [100, 185],
        "ch4": [270, 450],
        "ch5": [250, 450],
    },
    "147": {
        "ch3": [80, 180],
        "ch4": [250, 400],
        "ch5": [230, 400],
    },
    "150": {
        "ch3": [100, 185],
        "ch4": [230, 380],
        "ch5": [230, 400],
    },
}

x_cuts = {
    # Day 1
    "131": {
        "ch3": [-0.8, 3.0],
        "ch4": [-1.4, 0.4],
        "ch5": [-1.4, 0.4],
    },
    "134": {
        "ch3": [-0.8, 3.0],
        "ch4": [-1.4, 0.4],
        "ch5": [-1.4, 0.4],
    },
    # Day 2
    "135": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.2, -3.5],
        "ch5": [-4.2, -3.6],
    },
    "136": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.2, -3.5],
        "ch5": [-4.2, -3.6],
    },
    "137": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.2, -3.5],
        "ch5": [-4.2, -3.6],
    },
    "138": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.2, -3.5],
        "ch5": [-4.2, -3.6],
    },
    # Day 3
    "139": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.6, -3.4],
        "ch5": [-4.5, -3.5],
    },
    "140": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.6, -3.4],
        "ch5": [-4.5, -3.5],
    },
    # Day 4
    "143": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.6, -3.4],
        "ch5": [-4.5, -3.5],
    },
    "147": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.6, -3.4],
        "ch5": [-4.5, -3.5],
    },
    "150": {
        "ch3": [-3.4, 0.4],
        "ch4": [-4.6, -3.4],
        "ch5": [-4.5, -3.5],
    },
}

# Nominal cuts
y_cuts = {
    # Day 1
    "131": {
        "ch3": [3.0, 3.8],
        "ch4": [5.4, 6.4],
        "ch5": [6.5, 7.5],
    },
    "134": {
        "ch3": [3.0, 3.8],
        "ch4": [5.4, 6.4],
        "ch5": [6.5, 7.5],
    },
    # Day 2
    "135": {
        "ch3": [2.5, 3.3],
        "ch4": [5.1, 6.1],
        "ch5": [6.5, 7.5],
    },
    "136": {
        "ch3": [2.5, 3.3],
        "ch4": [5.1, 6.1],
        "ch5": [6.5, 7.5],
    },
    "137": {
        "ch3": [2.5, 3.3],
        "ch4": [5.1, 6.1],
        "ch5": [6.5, 7.5],
    },
    "138": {
        "ch3": [2.5, 3.3],
        "ch4": [5.1, 6.1],
        "ch5": [6.5, 7.5],
    },
    # Day 3
    "139": {
        "ch3": [2.6, 3.4],
        "ch4": [5.0, 6.3],
        "ch5": [6.3, 7.5],
    },
    "140": {
        "ch3": [2.6, 3.4],
        "ch4": [5.0, 6.3],
        "ch5": [6.5, 7.5],
    },
    # Day 4
    "143": {
        "ch3": [2.6, 3.4],
        "ch4": [5.0, 6.3],
        "ch5": [6.5, 7.5],
    },
    "147": {
        "ch3": [2.6, 3.4],
        "ch4": [5.0, 6.3],
        "ch5": [6.5, 7.5],
    },
    "150": {
        "ch3": [2.6, 3.4],
        "ch4": [5.0, 6.3],
        "ch5": [6.5, 7.5],
    },
}

# Test
# y_cuts = {
# "131": {
# "ch3": [3.35, 3.45],
# "ch4": [5.4, 6.4],
# "ch5": [6.5, 7.5],
# },
# "134": {
# "ch3": [3.0, 3.8],
# "ch4": [5.4, 6.4],
# "ch5": [6.5, 7.5],
# },
# "135": {
# "ch3": [2.85, 2.95],
# "ch4": [5.1, 6.1],
# "ch5": [6.5, 7.5],
# },
# "136": {
# "ch3": [2.85, 2.95],
# "ch4": [5.1, 6.1],
# "ch5": [6.5, 7.5],
# },
# "137": {
# "ch3": [2.85, 2.95],
# "ch4": [5.1, 6.1],
# "ch5": [6.5, 7.5],
# },
# "138": {
# "ch3": [2.85, 2.95],
# "ch4": [5.1, 6.1],
# "ch5": [6.5, 7.5],
# },
# "139": {
# "ch3": [2.95, 3.05],
# "ch4": [5.0, 6.3],
# "ch5": [6.3, 7.5],
# },
# "140": {
# "ch3": [2.95, 3.05],
# "ch4": [5.0, 6.3],
# "ch5": [6.5, 7.5],
# },
# }

deltaT_binning = {
    "ch3": [60, -10.7, -9],
    "ch4": [60, -10.2, -9.7],
    "ch5": [60, -10.2, -9.6],
}
