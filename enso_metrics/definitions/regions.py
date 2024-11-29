# -*- coding:UTF-8 -*-
# ---------------------------------------------------------------------------------------------------------------------#
# Definition of climate regions
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------#
# Import packages
# ---------------------------------------------------#
# basic python package
from copy import deepcopy
# ---------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
# Functions
# ---------------------------------------------------------------------------------------------------------------------#
def regions_param(region=True):
    dict_reference_regions = {
        "global_no_poles": {"long_name": "Global 60S-60N", "short_name": "GL60", "maskland": False, "maskocean": False,
                            "latitude": (-60., 60.), "longitude": (0., 360.)},
        "global": {"long_name": "Global", "short_name": "GLOB", "maskland": False, "maskocean": False,
                   "latitude": (-90., 90.), "longitude": (0., 360.)},
        "equatorial_atlantic": {"long_name": "Equatorial Atlantic", "short_name": "EATL", "maskland": True,
                                "maskocean": False, "latitude": (-5., 2.), "longitude": (-30., 10.)},
        "equatorial_pacific": {"long_name": "Equatorial Pacific", "short_name": "EPAC", "maskland": True,
                               "maskocean": False, "latitude": (-5., 5.), "longitude": (150., 270.)},
        "equatorial_pacific_LonExt": {"long_name": "Equatorial Pacific extended in longitude", "short_name": "EPXL",
                                      "maskland": True, "maskocean": False, "latitude": (-5., 5.),
                                      "longitude": (140., 280.)},
        "equatorial_pacific_LatExt": {"long_name": "Equatorial Pacific extended in latitude", "short_name": "EPYL",
                                      "maskland": True, "maskocean": False, "latitude": (-15., 15.),
                                      "longitude": (150., 270.)},
        "equatorial_pacific_LatExt2": {"long_name": "Equatorial Pacific extended in latitude and longitude",
                                       "short_name": "EPXYL", "maskland": True, "maskocean": False,
                                       "latitude": (-15., 15.), "longitude": (120., 285.)},
        "eastern_equatorial_pacific": {"long_name": "Eastern Equatorial Pacific", "short_name": "EEPAC",
                                       "maskland": True, "maskocean": False, "latitude": (-5., 5.),
                                       "longitude": (205., 280.)},
        "grad_east": {"long_name": "Eastern Tropical Pacific", "short_name": "ETPAC", "maskland": True,
                      "maskocean": False, "latitude": (-8, 8), "longitude": (205., 280.)},
        "grad_equ": {"long_name": "Eastern Near-Equatorial Pacific", "short_name": "ENPAC", "maskland": True,
                     "maskocean": False, "latitude": (-2.5, 2.5), "longitude": (210., 270.)},
        "grad_off": {"long_name": "Eastern Off-Equatorial Pacific", "short_name": "EOPAC", "maskland": True,
                     "maskocean": False, "latitude": (5., 10.), "longitude": (210., 270.)},
        "grad_west": {"long_name": "Western Tropical Pacific", "short_name": "WTPAC", "maskland": True,
                      "maskocean": False, "latitude": (-8, 8), "longitude": (130., 205.)},
        "grad_i": {"long_name": "ITCZ Gradient (N3NO - N3SO)", "short_name": "GRAI", "maskland": True,
                   "maskocean": False, "latitude": None, "longitude": None, "regions": ["itcz_ne", "itcz_se"]},
        "grad_x": {"long_name": "Zonal Gradient (WTPAC - ETPAC)", "short_name": "GRAX", "maskland": True,
                   "maskocean": False, "latitude": None, "longitude": None, "regions": ["grad_west", "grad_east"]},
        "grad_y": {"long_name": "Meridional Gradient (EOPAC - ENPAC)", "short_name": "GRAY", "maskland": True,
                   "maskocean": False, "latitude": None, "longitude": None, "regions": ["grad_off", "grad_equ"]},
        "itcz_ne": {"long_name": "Northeastern Equatorial Pacific (N3 north)", "short_name": "N3NO", "maskland": True,
                    "maskocean": False, "latitude": (5., 15.), "longitude": (210., 270.)},
        "itcz_nw": {"long_name": "Northwestern Equatorial Pacific (N4 north)", "short_name": "N4NO", "maskland": True,
                    "maskocean": False, "latitude": (5., 15.), "longitude": (160., 210.)},
        "itcz_se": {"long_name": "Southeastern Equatorial Pacific (N3 south)", "short_name": "N3SO", "maskland": True,
                    "maskocean": False, "latitude": (-15., -5.), "longitude": (210., 270.)},
        "itcz_sw": {"long_name": "Southwestern Equatorial Pacific (N4 south)", "short_name": "N4SO", "maskland": True,
                    "maskocean": False, "latitude": (-15., -5.), "longitude": (160., 210.)},
        "nino1+2": {"long_name": "Nino 1+2", "short_name": "N1+2", "maskland": True, "maskocean": False,
                    "latitude": (-10., 0.), "longitude": (270., 280.)},
        "nino3": {"long_name": "Nino 3", "short_name": "NIN3", "maskland": True, "maskocean": False,
                  "latitude": (-5., 5.), "longitude": (210., 270.)},
        "nino3_LatExt": {"long_name": "Nino 3 extended in latitude", "short_name": "N3YL", "maskland": True,
                         "maskocean": False, "latitude": (-15., 15.), "longitude": (210., 270.)},
        "nino3.4": {"long_name": "Nino 3.4", "short_name": "N3.4", "maskland": True, "maskocean": False,
                    "latitude": (-5., 5.), "longitude": (190., 240.)},
        "nino4": {"long_name": "Nino 4", "short_name": "NIN4", "maskland": True, "maskocean": False,
                  "latitude": (-5., 5.), "longitude": (160., 210.)},
        "nino4small": {"long_name": "Near-Equatorial Nino 4", "short_name": "WCPAC", "maskland": True,
                       "maskocean": False, "latitude": (-2., 2.), "longitude": (160., 210.)},
        "north_tropical_atlantic": {"long_name": "North Tropical Atlantic (NTA)", "short_name": "NTATL",
                                    "maskland": True, "maskocean": False, "latitude": (0., 30.),
                                    "longitude": (-60., 20.)},
        "southeastern_equatorial_indian": {"long_name": "Southeastern Equatorial Indian", "short_name": "SEEI",
                                           "maskland": True, "maskocean": False, "latitude": (-10., 0.),
                                           "longitude": (90., 110.)},
        "tropic": {"long_name": "Tropic", "short_name": "TROP", "maskland": True, "maskocean": False,
                   "latitude": (-30., 30.), "longitude": (0., 360.)},
        "tropical_atlantic": {"long_name": "Tropical Atlantic", "short_name": "TATL", "maskland": True,
                              "maskocean": False, "latitude": (-30., 30.), "longitude": (-60., 20.)},
        "tropical_indian": {"long_name": "Tropical Indian", "short_name": "TIND", "maskland": True, "maskocean": False,
                            "latitude": (-30., 30.), "longitude": (40., 120.)},
        "tropical_pacific": {"long_name": "Tropical Pacific", "short_name": "TPAC", "maskland": True,
                             "maskocean": False, "latitude": (-30., 30.), "longitude": (120., 280.)},
        "western_equatorial_indian": {"long_name": "Western Equatorial Indian (WEI)", "short_name": "WEIN",
                                      "maskland": True, "maskocean": False, "latitude": (-10., 10.),
                                      "longitude": (50., 70.)},
        "western_equatorial_pacific": {"long_name": "Western Equatorial Pacific", "short_name": "WEPAC",
                                       "maskland": True, "maskocean": False, "latitude": (-5., 5.),
                                       "longitude": (120., 205.)},
        # IPCC climate reference regions defined in Iturbide et al. (2020; https://doi.org/10.5194/essd-12-2959-2020)
        "GIC": {
            "long_name": "POLAR, Land, Greenland / Iceland", "polygon": True, "maskland": False, "maskocean": True,
            "latitude": (62.0, 62.0, 58.0, 58.0, 85.0, 85.0), "longitude": (350.0, 322.0, 318.0, 310.0, 278.0, 350.0)},
        "NWN": {"long_name": "NORTH - AMERICA, Land, N.W.North - America", "maskland": False, "maskocean": True,
                "latitude": (50.0, 50.0, 58.0, 52.5, 72.6, 72.6, 77.6, 81.0), "polygon": True,
                "longitude": (255.0, 230.0, 217.0, 192.0, 192.0, 231.0, 235.0, 255.0)},
        "NEN": {"long_name": "NORTH - AMERICA, Land, N.E.North - America", "maskland": False, "maskocean": True,
                "latitude": (50.0, 58.0, 85.0, 81.0, 50.0), "polygon": True,
                "longitude": (310.0, 310.0, 278.0, 255.0, 255.0)},
        "WNA": {"long_name": "NORTH - AMERICA, Land, W.North - America", "maskland": False, "maskocean": True,
                "latitude": (50.0, 33.8, 33.8, 50.0), "longitude": (230.0, 237.5, 255.0, 255.0), "polygon": True},
        "CNA": {"long_name": "NORTH - AMERICA, Land, C.North - America", "maskland": False, "maskocean": True,
                "latitude": (50.0, 25.0, 33.8, 50.0), "longitude": (270.0, 270.0, 255.0, 255.0), "polygon": True},
        "ENA": {"long_name": "NORTH - AMERICA, Land, E.North - America", "maskland": False, "maskocean": True,
                "latitude": (25.0, 25.0, 50.0, 50.0, 31.0), "polygon": True,
                "longitude": (290.0, 270.0, 270.0, 310.0, 283.0)},
        "NCA": {"long_name": "CENTRAL - AMERICA, Land, N.Central - America", "maskland": False, "maskocean": True,
                "latitude": (25.0, 16.0, 33.8, 33.8), "longitude": (270.0, 255.5, 237.5, 255.0), "polygon": True},
        "SCA": {"long_name": "CENTRAL - AMERICA, Land, S.Central - America", "maskland": False, "maskocean": True,
                "latitude": (12.0, 2.2, 16.0, 25.0), "longitude": (285.0, 276.6, 255.5, 270.0), "polygon": True},
        "CAR_L": {"long_name": "CENTRAL - AMERICA, Land - Ocean, Caribbean", "maskland": False, "maskocean": True,
                  "latitude": (12.0, 25.0, 25.0, 12.0), "longitude": (285.0, 270.0, 290.0, 305.0), "polygon": True},
        "CAR_O": {"long_name": "CENTRAL - AMERICA, Land - Ocean, Caribbean", "maskland": True, "maskocean": False,
                  "latitude": (12.0, 25.0, 25.0, 12.0), "longitude": (285.0, 270.0, 290.0, 305.0), "polygon": True},
        "NWS": {"long_name": "SOUTH - AMERICA, Land, N.W.South - America", "maskland": False, "maskocean": True,
                "latitude": (12.0, 2.2, -10.0, -15.0, -15.0, 12.0), "polygon": True,
                "longitude": (285.0, 276.6, 276.6, 281.0, 288.0, 288.0)},
        "NSA": {"long_name": "SOUTH - AMERICA, Land, N.South - America", "maskland": False, "maskocean": True,
                "latitude": (12.0, -8.0, -8.0, 7.6, 12.0), "polygon": True,
                "longitude": (288.0, 288.0, 310.0, 310.0, 305.0)},
        "NES": {"long_name": "SOUTH - AMERICA, Land, N.E.South - America", "maskland": False, "maskocean": True,
                "latitude": (-20.0, -20.0, 0.0, 0.0), "longitude": (-34.0, -50.0, -50.0, -34.0), "polygon": True},
        "SAM": {"long_name": "SOUTH - AMERICA, Land, South - American - Monsoon", "maskland": False, "maskocean": True,
                "latitude": (-20.0, -15.0, -8.0, -8.0, -20.0), "polygon": True,
                "longitude": (293.6, 288.0, 288.0, 310.0, 310.0)},
        "SWS": {"long_name": "SOUTH - AMERICA, Land, S.W.South - America", "maskland": False, "maskocean": True,
                "latitude": (-15.0, -20.0, -47.0, -47.0, -20.0, -15.0), "polygon": True,
                "longitude": (288.0, 293.6, 288.5, 281.0, 285.4, 281.0)},
        "SES": {"long_name": "SOUTH - AMERICA, Land, S.E.South - America", "maskland": False, "maskocean": True,
                "latitude": (-20.0, -40.0, -40.0, -20.0), "longitude": (326.0, 304.0, 289.8, 293.6), "polygon": True},
        "SSA": {"long_name": "SOUTH - AMERICA, Land, S.South - America", "maskland": False, "maskocean": True,
                "latitude": (-56.0, -47.0, -47.0, -40.0, -40.0, -56.0),
                "longitude": (281.0, 281.0, 288.5, 289.8, 304.0, 304.0), "polygon": True},
        "NEU": {"long_name": "EUROPE, Land, N.Europe", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (48.0, 72.6, 72.6, 61.3), "longitude": (-10.0, -10.0, 40.0, 40.0)},
        "WCE": {"long_name": "EUROPE, Land, West & Central - Europe", "maskland": False, "maskocean": True,
                "latitude": (45.0, 48.0, 61.3, 45.0), "longitude": (-10.0, -10.0, 40.0, 40.0), "polygon": True},
        "EEU": {"long_name": "EUROPE, Land, E.Europe", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (45.0, 45.0, 65.0, 65.0), "longitude": (40.0, 60.0, 60.0, 40.0)},
        "MED_L": {"long_name": "EUROPE - AFRICA, Land - Ocean, Mediterranean", "maskland": False, "maskocean": True,
                  "latitude": (30.0, 45.0, 45.0, 30.0), "longitude": (-10.0, -10.0, 40.0, 40.0), "polygon": True},
        "MED_O": {"long_name": "EUROPE - AFRICA, Land - Ocean, Mediterranean", "maskland": True, "maskocean": False,
                  "latitude": (30.0, 45.0, 45.0, 30.0), "longitude": (-10.0, -10.0, 40.0, 40.0), "polygon": True},
        "SAH": {"long_name": "AFRICA, Land, Sahara", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (14.7, 30.0, 30.0, 14.7), "longitude": (-20.0, -20.0, 33.0, 42.1)},
        "WAF": {"long_name": "AFRICA, Land, Western - Africa", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (7.6, 14.7, 14.7, 0.0), "longitude": (-20.0, -20.0, 15.0, 8.0)},
        "CAF": {"long_name": "AFRICA, Land, Central - Africa", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (-10.0, 0.0, 14.7, 14.7, -10.0), "longitude": (8.0, 8.0, 15.0, 27.0, 27.0)},
        "NEAF": {"long_name": "AFRICA, Land, N.Eastern - Africa", "polygon": True, "maskland": False, "maskocean": True,
                 "latitude": (2.3, 14.7, 14.7, 12.0, 15.0, 7.0, 2.3),
                 "longitude": (27.0, 27.0, 42.1, 43.7, 53.0, 53.0, 46.5)},
        "SEAF": {"long_name": "AFRICA, Land, S.Eastern - Africa", "maskland": False, "maskocean": True,
                 "latitude": (-10.0, -10.0, 2.3, 2.3), "longitude": (27.0, 46.5, 46.5, 27.0), "polygon": True},
        "WSAF": {"long_name": "AFRICA, Land, W.Southern - Africa", "maskland": False, "maskocean": True,
                 "latitude": (-36.0, -36.0, -10.0, -10.0), "longitude": (8.0, 25.0, 25.0, 8.0), "polygon": True},
        "ESAF": {"long_name": "AFRICA, Land, E.Southern - Africa", "maskland": False, "maskocean": True,
                 "latitude": (-10.0, -36.0, -36.0, -10.0), "longitude": (25.0, 25.0, 31.0, 46.5), "polygon": True},
        "MDG": {"long_name": "AFRICA, Land, Madagascar", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (-27.0, -10.0, -10.0, -27.0), "longitude": (36.2, 46.5, 53.0, 53.0)},
        "RAR": {"long_name": "ASIA, Land, Russian - Arctic", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (65.0, 72.6, 82.0, 72.6, 65.0), "longitude": (40.0, 40.0, 94.0, 192.0, 192.0)},
        "WSB": {"long_name": "ASIA, Land, W.Siberia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (45.0, 45.0, 65.0, 65.0), "longitude": (60.0, 90.0, 90.0, 60.0)},
        "ESB": {"long_name": "ASIA, Land, E.Siberia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (45.0, 45.0, 65.0, 65.0), "longitude": (90.0, 130.0, 130.0, 90.0)},
        "RFE": {"long_name": "ASIA, Land, Russian - Far - East", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (45.0, 65.0, 65.0, 59.9, 50.0, 45.0),
                "longitude": (130.0, 130.0, 180.0, 180.0, 157.0, 152.0)},
        "WCA": {"long_name": "ASIA, Land, W.C.Asia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (30.0, 45.0, 45.0, 30.0, 30.0, 23.5, 30.0),
                "longitude": (40.0, 40.0, 75.0, 75.0, 60.0, 60.0, 47.6)},
        "ECA": {"long_name": "ASIA, Land, E.C.Asia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (37.0, 45.0, 45.0, 37.0), "longitude": (75.0, 75.0, 117.0, 108.0)},
        "TIB": {"long_name": "ASIA, Land, Tibetan - Plateau", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (30.0, 37.0, 37.0, 30.0, 26.0), "longitude": (75.0, 75.0, 100.0, 100.0, 88.0)},
        "EAS": {"long_name": "ASIA, Land, E.Asia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (19.5, 37.0, 37.0, 45.0, 45.0, 25.0, 19.5),
                "longitude": (100.0, 100.0, 108.0, 117.0, 152.0, 132.0, 132.0)},
        "ARP": {"long_name": "ASIA, Land, Arabian - Peninsula", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (30.0, 30.0, 23.5, 19.5, 15.0, 12.0), "longitude": (33.0, 47.6, 60.0, 60.0, 53.0, 43.7)},
        "SAS": {"long_name": "ASIA, Land, S.Asia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (23.5, 30.0, 30.0, 26.0, 30.0, 19.5, 19.5, 19.5, 7.0, 7.0, 19.5, 23.5),
                "longitude": (60.0, 60.0, 75.0, 88.0, 100.0, 100.0, 95.0, 87.0, 79.0, 76.0, 70.0, 66.5)},
        "SEA_L": {"long_name": "ASIA, Land - Ocean, S.E.Asia", "polygon": True, "maskland": False, "maskocean": True,
                  "latitude": (-10.0, 19.5, 19.5, 5.0, -10.0), "longitude": (93.0, 93.0, 132.0, 132.0, 155.0)},
        "SEA_O": {"long_name": "ASIA, Land - Ocean, S.E.Asia", "polygon": True, "maskland": True, "maskocean": False,
                  "latitude": (-10.0, 19.5, 19.5, 5.0, -10.0), "longitude": (93.0, 93.0, 132.0, 132.0, 155.0)},
        "NAU": {"long_name": "OCEANIA, Land, N.Australia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (-20.0, -10.0, -10.0, -20.0), "longitude": (110.0, 110.0, 155.0, 155.0)},
        "CAU": {"long_name": "OCEANIA, Land, C.Australia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (-30.0, -20.0, -20.0, -32.9, -30.0), "longitude": (110.0, 110.0, 145.5, 145.5, 140.0)},
        "EAU": {"long_name": "OCEANIA, Land, E.Australia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (-32.9, -20.0, -20.0, -38.0), "longitude": (145.5, 145.5, 155.0, 155.0)},
        "SAU": {"long_name": "OCEANIA, Land, S.Australia", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (-36.0, -30.0, -30.0, -38.0, -50.0), "longitude": (110.0, 110.0, 140.0, 155.0, 155.0)},
        "NZ": {"long_name": "OCEANIA, Land, New - Zealand", "polygon": True, "maskland": False, "maskocean": True,
               "latitude": (-50.0, -30.0, -30.0, -50.0), "longitude": (155.0, 155.0, 180.0, 180.0)},
        "EAN": {"long_name": "POLAR, Land, E.Antarctica", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (-90.0, -83.0, -83.0, -75.0, -75.0, -64.0, -64.0, -90.0),
                "longitude": (-180.0, -180.0, -56.0, -56.0, -25.0, 5.0, 180.0, 180.0)},
        "WAN": {"long_name": "POLAR, Land, W.Antarctica", "polygon": True, "maskland": False, "maskocean": True,
                "latitude": (-83.0, -70.0, -70.0, -62.0, -62.0, -83.0),
                "longitude": (180.0, 180.0, 280.0, 295.0, 304.0, 304.0)},
        "ARO": {"long_name": "ARCTIC, Ocean, Arctic - Ocean", "polygon": True, "maskland": True, "maskocean": False,
                "latitude": (90.0, 73.8, 72.6, 72.6, 77.6, 85.0, 85.0, 72.6, 72.6, 82.0, 73.8, 90.0),
                "longitude": (-180.0, -180.0, -168.0, -129.0, -125.0 , -82.0, -10.0, -10.0, 40.0, 94.0, 180.0, 180.0)},
        "NPO": {"long_name": "PACIFIC, Ocean, N.Pacific - Ocean", "polygon": True, "maskland": True, "maskocean": False,
                "latitude": (7.6, 25.0, 50.0, 59.9, 65.0, 65.0, 52.5, 58.0, 50.0, 33.8, 16.0, 7.6),
                "longitude": (132.0, 132.0, 157.0, 180.0, 180.0, 192.0, 192.0, 217.0, 230.0, 237.5, 255.5, 268.3)},
        "EPO": {"long_name": "PACIFIC, Ocean, Equatorial.Pacific - Ocean", "maskland": True, "maskocean": False,
                "latitude": (-10.0, 5.0, 7.6, 7.6, 2.2, -10.0), "polygon": True,
                "longitude": (155.0, 132.0, 132.0, 268.3, 276.6, 276.6)},
        "SPO": {"long_name": "PACIFIC, Ocean, S.Pacific - Ocean", "polygon": True, "maskland": True, "maskocean": False,
                "latitude": (-30.0, -10.0, -10.0, -20.0, -47.0, -56.0, -56.0, -30.0),
                "longitude": (155.0, 155.0, 276.6, 284.4, 281.0, 281.0, 180.0, 180.0)},
        "NAO": {"long_name": "ATLANTIC, Ocean, N.Atlantic - Ocean", "maskland": True, "maskocean": False,
                "latitude": (7.6, 31.0, 50.0, 58.0, 58.0, 62.0, 62.0, 30.0, 30.0, 7.6), "polygon": True,
                "longitude": (310.0, 283.0, 310.0, 310.0, 312.0, 322.0, 350.0, 350.0, 340.0, 340.0)},
        "EAO": {"long_name": "ATLANTIC, Ocean, Equatorial.Atlantic - Ocean", "maskland": True, "maskocean": False,
                "latitude": (-10.0, 0.0, 0.0, 7.6, 7.6, 0.0, -10.0), "polygon": True,
                "longitude": (-34.0, -34.0, -50.0, -50.0, -20.0, 8.0, 8.0)},
        "SAO": {"long_name": "ATLANTIC, Ocean, S.Atlantic - Ocean", "maskland": True, "maskocean": False,
                "latitude": (-56.0, -40.0, -20.0, -10.0, -10.0, -36.0), "polygon": True,
                "longitude": (-56.0, -56.0, -34.0, -34.0, 8.0, 8.0)},
        "ARS": {"long_name": "INDIAN, Ocean, Arabian - Sea", "polygon": True, "maskland": True, "maskocean": False,
                "latitude": (7.0, 15.0, 19.5, 23.5, 23.5, 19.5, 7.0),
                "longitude": (53.0, 53.0, 60.0, 60.0, 66.5, 70.0, 76.0)},
        "BOB": {"long_name": "INDIAN, Ocean, Bay - of - Bengal", "polygon": True, "maskland": True, "maskocean": False,
                "latitude": (7.0, 19.5, 19.5, 7.0), "longitude": (79.0, 87.0, 93.0, 93.0)},
        "EIO": {"long_name": "INDIAN, Ocean, Equatorial.Indic - Ocean", "maskland": True, "maskocean": False,
                "latitude": (-10.0, 2.3, 7.0, 7.0, -10.0), "polygon": True,
                "longitude": (46.5, 46.5, 53.0, 93.0, 93.0)},
        "SIO": {"long_name": "INDIAN, Ocean, S.Indic - Ocean", "polygon": True, "maskland": True, "maskocean": False,
                "latitude": (-27.0, -27.0, -10.0, -10.0, -36.0, -36.0),
                "longitude": (36.2, 53.0, 53.0, 110.0, 110.0, 31.0)},
        "SOO": {"long_name": "SOUTHERN, Ocean, Southern - Ocean", "polygon": True, "maskland": True, "maskocean": False,
                "latitude":
                    (-56.0, -70.0, -70.0, -62.0, -62.0, -75.0, -75.0, -64.0, -64.0, -50.0, -50, -36.0, -36.0, -56.0),
                "longitude":
                    (-180.0, -180.0, -80.0, -65.0, -56.0, -56.0, -25.0, 5.0, 180.0, 180.0, 155.0, 110.0, 8.0, -56.0)},
        # # IPCC SREX modified by Perry et al. (2020; https://doi.org/10.1007/s00382-019-05006-6)
        # "ALA": {"long_name": "Alaska/N.W. Canada", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (60., 72.6), "longitude": (192., 255.)},
        # "AMZ": {"long_name": "Amazon", "polygon": True, "maskland": False, "maskocean": True,
        #         "latitude": (-20., -1.2, 11.4, 11.4, -20.), "longitude": (293.6, 280.3, 291.2, 310., 310.)},
        # "CAM": {"long_name": "Central America/Mexico", "polygon": True, "maskland": False, "maskocean": True,
        #         "latitude": (11.4, -1.2, 25., 25.), "longitude": (291.2, 280.3, 241.7, 269.7)},
        # "EAF": {"long_name": "East Africa", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (-11.4, 15.), "longitude": (25., 52.)},
        # "EAS": {"long_name": "East Asia", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (20., 50.), "longitude": (100., 145.)},
        # "ERU": {"long_name": "East Russia", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (50., 75.), "longitude": (130., 190.)},
        # "CGI": {"long_name": "Canada Greenland Iceland", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (60., 85.), "longitude": (255., 350.)},
        # "NEB": {"long_name": "North-East Brazil", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (-20., 0.), "longitude": (310., 326.)},
        # "NAU": {"long_name": "North Australia", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (-30., -10.), "longitude": (110., 155.)},
        # "NNA": {"long_name": "Northern North America", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (40., 60.), "longitude": (230., 300.)},
        # "NZE": {"long_name": "New Zealand", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (-50., -30.), "longitude": (165., 180.)},
        # "SAF": {"long_name": "Southern Africa", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (-35., -11.4), "longitude": (0., 52.)},
        # "SAH": {"long_name": "Sahara", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (15., 30.), "longitude": (340., 400.)},
        # "SAS": {"long_name": "South Asia", "polygon": True, "maskland": False, "maskocean": True,
        #         "latitude": (5., 30., 30., 20., 20., 5.), "longitude": (60., 60., 100., 100., 95., 95.)},
        # "SAU": {"long_name": "South Australia", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (-50., -30.), "longitude": (110., 155.)},
        # "SEA": {"long_name": "South East Asia", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (-10., 20.), "longitude": (95., 155.)},
        # "SNA": {"long_name": "Southern North America", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (25., 40.), "longitude": (230., 300.)},
        # "SSA": {"long_name": "Southeastern South America", "polygon": True, "maskland": False, "maskocean": True,
        #         "latitude": (-20., -56.7, -56.7, -50., -20.), "longitude": (320.6, 320.6, 292.7, 287.9, 293.6)},
        # "WAF": {"long_name": "West Africa", "polygon": False, "maskland": False, "maskocean": True,
        #         "latitude": (-11.4, 15.), "longitude": (340., 385.)},
        # "WSA": {"long_name": "West Coast South America", "polygon": True, "maskland": False, "maskocean": True,
        #         "latitude": (-1.2, -20., -50., -56.7, -56.7, 0.5),
        #         "longitude": (280.3, 293.6, 287.9, 292.7, 278., 278.)},
    }
    if region in list(dict_reference_regions.keys()):
        output = deepcopy(dict_reference_regions[region])
    else:
        output = deepcopy(dict_reference_regions)
    return output
# ---------------------------------------------------------------------------------------------------------------------#
