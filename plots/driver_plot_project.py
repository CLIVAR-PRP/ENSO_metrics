# -*- coding:UTF-8 -*-
# ---------------------------------------------------#
# Aim of the program:
#      Create plots for ENSO_metrics
# ---------------------------------------------------#


# ---------------------------------------------------#
# Import the right packages
# ---------------------------------------------------#
from copy import deepcopy
from glob import iglob as GLOBiglob
import json
from os.path import join as OSpath__join
# ENSO_metrics functions
from EnsoCollectionsLib import defCollection
from EnsoMetricPlot import main_plotter


# ---------------------------------------------------#
# Arguments
# ---------------------------------------------------#
metric_collection = "ENSO_tel"
experiment = "historical"  # "piControl" #
list_project = ["cmip5", "cmip6"]
my_project = ["CMIP"]#["CMIP", "select8"]
big_ensemble = False  # True  #
dict_selection = {
    #
    # EnsoSeasonality within 25% of obs &
    #
    # BiasSstLonRmse less than 1°C
    "select1": ["ACCESS1-0", "BCC-ESM1", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "GISS-E2-R", "GISS-E2-R-CC",
                "MIROC4h", "NorESM1-M"],
    # BiasPrLonRmse less than 1 mm/day
    "select2": ["BCC-CSM1-1", "BCC-ESM1", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2",
                "NorESM1-M", "NorESM1-ME"],
    # SeasonalSstLonRmse less than 0.25°C
    "select3": ["BCC-CSM1-1", "BCC-ESM1", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2",
                "NorESM1-M", "NorESM1-ME"],
    # SeasonalPrLonRmse less than Q2 (median)
    "select4": ["CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2", "GISS-E2-R", "GISS-E2-H-CC",
                "NorESM1-M", "NorESM1-ME"],
    # EnsoFbSstTaux within 50% of obs
    "select5": ["CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2", "GISS-E2-R", "GISS-E2-R-CC",
                "MIROC4h", "NorESM1-M", "NorESM1-ME"],
    # EnsoFbSstThf within 50% of obs
    "select6": ["ACCESS1-0", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2", "GISS-E2-R",
                "GISS-E2-R-CC", "MIROC4h", "NorESM1-M", "NorESM1-ME"],
    # EnsoFbSstSwr within 100% of obs
    "select7": ["ACCESS1-0", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CNRM-CM5", "CNRM-CM5-2", "GISS-E2-R",
                "GISS-E2-R-CC", "MIROC4h", "NorESM1-M", "NorESM1-ME"],
    #
    # EnsoFbSstTaux & EnsoFbSstThf within Q1 (approximately)
    #
    # "select8": ["CCSM4", "CESM1-FASTCHEM", "CNRM-CM5", "GFDL-ESM2M", "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-1-H",
    #             "GISS-E2-R", "GISS-E2-R-CC", "MIROC6"],
    "select8": ["CCSM4", "CESM1-FASTCHEM", "GFDL-ESM2M", "GISS-E2-1-G", "GISS-E2-1-G-CC", "GISS-E2-1-H",
                "GISS-E2-R", "GISS-E2-R-CC", "MIROC-ES2L", "MIROC6"],
    #
    # EnsoFbSstTaux within Q2 & EnsodSstOce_2 within Q1
    #
    "select9": ["CESM1-BGC", "CESM2", "CNRM-CM5", "EC-Earth3-Veg", "FGOALS-g2", "GISS-E2-R-CC", "MIROC4h", "NorESM1-ME",
                "SAM0-UNICON"],
    #
    # EnsoFbSstThf within Q2 & EnsodSstOce_2 within Q1
    #
    "select10": ["ACCESS1-0", "CESM1-BGC", "CESM2", "CNRM-CM5", "EC-Earth3-Veg", "FGOALS-g2", "GFDL-CM3", "GISS-E2-H",
                 "GISS-E2-R-CC", "MIROC4h", "NorESM1-ME"],
    #
    # EnsoSeasonality within Q2 & EnsoSstSkew within Q2
    #
    "select11": ["BCC-CSM1-1-M", "CCSM4", "CESM1-CAM5", "CESM1-FASTCHEM", "CESM1-WACCM", "CESM2-WACCM", "CNRM-CM5",
                 "GFDL-ESM2G", "GISS-E2-R", "MIROC4h", "MIROC6", "NorESM1-ME"],
}

path_main = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2018_06_ENSO_metrics/2019_10_report"
path_in = OSpath__join(path_main, "Data_grouped")
# path_out = OSpath__join(path_main, "Plots_v6")
# path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_12_09_AGU/Poster"
path_out = "/Users/yannplanton/Documents/Yann/Fac/2016_2018_postdoc_LOCEAN/2019_10_ENSO_evaluation/v02"

expe = "hist" if experiment == "historical" else "pi"
pattern = "cmip?_" + experiment + "_" + metric_collection + "_v2019????_"


# ---------------------------------------------------#
# Main
# ---------------------------------------------------#
# read json file
data_json = dict()
for proj in list_project:
    lpath = OSpath__join(path_in, proj + "/" + experiment + "/" + metric_collection)
    lname = proj + "_" + experiment + "_" + metric_collection + "_v2019????.json"
    filename_js = list(GLOBiglob(OSpath__join(lpath, lname)))[0]
    with open(filename_js) as ff:
        data = json.load(ff)
    ff.close()
    data_json[proj] = data["RESULTS"]["model"]
    del data, ff, filename_js, lpath, lname
alt_json = dict()
for proj in list_project:
    lpath = deepcopy(path_in)  # OSpath__join(path_in, proj + "/" + experiment + "/" + metric_collection)
    lname = proj + "_" + experiment + "_" + metric_collection + "_v2019*.json"
    filename_js = list(GLOBiglob(OSpath__join(lpath, lname)))[0]
    with open(filename_js) as ff:
        data = json.load(ff)
    ff.close()
    alt_json[proj] = data["RESULTS"]["model"]
    del data, ff, filename_js, lpath, lname
# loop on metrics
list_metrics = sorted(defCollection(metric_collection)['metrics_list'].keys(), key=lambda v: v.upper())
if metric_collection == "ENSO_perf":
    # list_metrics = [
    #     'BiasPrLatRmse', 'BiasPrLonRmse', 'BiasSstLatRmse', 'BiasSstLonRmse', 'BiasTauxLatRmse', 'BiasTauxLonRmse',
    #     'EnsoAmpl', 'EnsoDuration', 'EnsoPrTsRmse', 'EnsoSeasonality', 'EnsoSstLonRmse', 'EnsoSstSkew', 'EnsoSstTsRmse',
    #     'EnsoTauxTsRmse', 'NinaSstDur_1', 'NinaSstDur_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstTsRmse_1',
    #     'NinaSstTsRmse_2', 'NinoSstDiversity_1', 'NinoSstDiversity_2', 'NinoSstDur_1', 'NinoSstDur_2',
    #     'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstTsRmse_1', 'NinoSstTsRmse_2', 'SeasonalPrLatRmse',
    #     'SeasonalPrLonRmse', 'SeasonalSstLatRmse', 'SeasonalSstLonRmse', 'SeasonalTauxLatRmse', 'SeasonalTauxLonRmse']
    # list_metrics = [
    #     'BiasPrLatRmse', 'BiasPrLonRmse', 'BiasSstLatRmse', 'BiasSstLonRmse', 'BiasTauxLatRmse', 'BiasTauxLonRmse',
    #     'EnsoAmpl', 'EnsoDuration', 'EnsoPrTsRmse', 'EnsoSeasonality', 'EnsoSstLonRmse', 'EnsoSstSkew', 'EnsoSstTsRmse',
    #     'EnsoTauxTsRmse', 'SeasonalPrLatRmse',
    #     'SeasonalPrLonRmse', 'SeasonalSstLatRmse', 'SeasonalSstLonRmse', 'SeasonalTauxLatRmse', 'SeasonalTauxLonRmse']
    list_metrics = ['BiasPrLatRmse']
elif metric_collection == "ENSO_proc":
    # list_metrics = ['EnsoAmpl', 'EnsodSstOce_1', 'EnsodSstOce_2', 'EnsoFbSshSst', 'EnsoFbSstLhf', 'EnsoFbSstLwr',
    #                 'EnsoFbSstShf', 'EnsoFbSstSwr', 'EnsoFbSstTaux', 'EnsoFbSstThf', 'EnsoFbTauxSsh']
    # list_metrics = ['EnsoAmpl', 'EnsodSstOce_1', 'EnsodSstOce_2', 'EnsoFbSstLhf', 'EnsoFbSstLwr',
    #                 'EnsoFbSstShf', 'EnsoFbSstSwr', 'EnsoFbSstTaux', 'EnsoFbSstThf']
    # list_metrics = ['EnsodSstOce_1', 'EnsodSstOce_2', 'EnsoFbSshSst', 'EnsoFbSstLhf', 'EnsoFbSstLwr',
    #                 'EnsoFbSstShf', 'EnsoFbSstSwr', 'EnsoFbSstTaux', 'EnsoFbSstThf', 'EnsoFbTauxSsh']
    list_metrics = ['EnsoFbTauxSsh']  # ['EnsoFbSstTaux']  # ['EnsoFbSstThf']  #
elif metric_collection == "ENSO_tel":
    # list_metrics = [
    #     'EnsoAmpl', 'EnsoPrMap', 'EnsoSlpMap', 'EnsoSstMap', 'NinaPrMap_1', 'NinaPrMap_2', 'NinaSlpMap_1',
    #     'NinaSlpMap_2', 'NinaSstLonRmse_1', 'NinaSstLonRmse_2', 'NinaSstMap_1', 'NinaSstMap_2', 'NinoPrMap_1',
    #     'NinoPrMap_2', 'NinoSlpMap_1', 'NinoSlpMap_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstMap_1',
    #     'NinoSstMap_2']
    list_metrics = ['EnsoPrMap']  # ['EnsoPrMap', 'EnsoSlpMap', 'EnsoSstMap']
    # for selection2:
    # list_metrics = [
    #     'NinaSstMap_1', 'NinaSstMap_2', 'NinoPrMap_1',
    #     'NinoPrMap_2', 'NinoSstLonRmse_1', 'NinoSstLonRmse_2', 'NinoSstMap_1',
    #     'NinoSstMap_2']
# loop on metrics
for met in list_metrics:
    print met
    my_data = {"diagnostic_values": {}, "metric_values": {}, "filenames": {}, "models": {}}
    if big_ensemble is False:
        # one ensemble for CMIP5 and one for CMIP6
        for proj in list_project:
            lpath = OSpath__join(path_in, proj + "/" + experiment + "/" + metric_collection)
            list_models = sorted(data_json[proj].keys(), key=lambda v: v.upper())
            # # one model has wrong wind direction
            # if "BCC-ESM1" in list_models:
            #     if "Taux" in met:
            #         list_models.remove("BCC-ESM1")
            # several models do not work...
            if "BNU-ESM" in list_models:
                list_models.remove("BNU-ESM")
            if "HadCM3" in list_models:
                list_models.remove("HadCM3")
            if metric_collection == "ENSO_tel":
                if "GFDL-CM4" in list_models:
                    list_models.remove("GFDL-CM4")
            dict1 = dict()
            tab1, tab2 = list(), list()
            for mod in list_models:
                # NetCDF file names
                try: tab1.append(list(GLOBiglob(OSpath__join(lpath, pattern + mod + "_" + met + "_???.nc")))[0])
                # try: tab1.append(list(GLOBiglob(OSpath__join(lpath, pattern + mod + "_" + met + "_yyp.nc")))[0])
                except: pass
                else:
                    # models
                    tab2.append(mod)
                    # diagnostic values
                    if metric_collection == "ENSO_tel" and "Map" in met:
                        dict1[mod] = None
                    else:
                        dict1[mod] = data_json[proj][mod]["value"][met]["diagnostic"][mod]["value"]
            my_data["diagnostic_values"][proj] = dict1
            my_data["filenames"][proj] = tab1
            my_data["models"][proj] = tab2
            if proj == list_project[0]:
                if metric_collection == "ENSO_tel" and "Map" in met:
                    tab = None
                else:
                    tab = dict((key1, data_json[proj][list_models[0]]["value"][met]["diagnostic"][key1]["value"])
                               for key1 in data_json[proj][list_models[0]]["value"][met]["diagnostic"].keys()
                               if key1 != list_models[0])
                my_data["diagnostic_values"]["obs"] = tab
            # metric values
            if metric_collection == "ENSO_tel" and "Map" in met:
                dicttmp1 = dict()
                for key1 in data_json[proj][list_models[0]]["value"][met + "Corr"]["metric"].keys():
                    dicttmp2 = dict()
                    for mod in my_data["models"][proj]:
                        dicttmp2[mod] = [data_json[proj][mod]["value"][met+suffix]["metric"][key1]["value"]
                                         if key1 in data_json[proj][mod]["value"][met+suffix]["metric"].keys() else None
                                         for suffix in ["Corr", "Rmse"]]
                    dicttmp1[key1] = dicttmp2
                    del dicttmp2
                my_data["metric_values"][proj] = dicttmp1
                del dicttmp1
            else:
                # my_data["metric_values"][proj] =\
                #     dict((key1, dict((mod, data_json[proj][mod]["value"][met]["metric"][key1]["value"]
                #                       if key1 in data_json[proj][mod]["value"][met]["metric"].keys() else None)
                #                      for mod in my_data["models"][proj]))
                #          for key1 in data_json[proj][list_models[0]]["value"][met]["metric"].keys())
                my_data["metric_values"][proj] = \
                    dict((key1, dict((mod, alt_json[proj][mod]["value"][met]["metric"][key1]["value"]
                                      if key1 in alt_json[proj][mod]["value"][met]["metric"].keys() else None)
                                     for mod in my_data["models"][proj]))
                         for key1 in alt_json[proj][list_models[0]]["value"][met]["metric"].keys())
            del lpath
            print str().ljust(5) + proj + " " + str(len(tab1)).zfill(2) + " file(s) " + \
                  str(len(list_models)).zfill(2) + " name(s)"
            # print tab1
            # print list_models
    else:
        # two project in one big ensemble with CMIP5 and CMIP6
        for grp in my_project:
            dict1 = dict()
            tab1, tab2 = list(), list()
            for proj in list_project:
                lpath = OSpath__join(path_in, proj + "/" + experiment + "/" + metric_collection)
                list_models = sorted(data_json[proj].keys(), key=lambda v: v.upper())
                if grp in dict_selection.keys():
                    my_mod = dict_selection[grp]
                    list_models = [mod for mod in list_models if mod in my_mod]
                    del my_mod
                # # one model has wrong wind direction
                # if "BCC-ESM1" in list_models and "Taux" in met:
                #     list_models.remove("BCC-ESM1")
                # several models do not work...
                if "BNU-ESM" in list_models:
                    list_models.remove("BNU-ESM")
                if "HadCM3" in list_models:
                    list_models.remove("HadCM3")
                tmp = list()
                for mod in list_models:
                    # NetCDF file names
                    try: tab1.append(list(GLOBiglob(OSpath__join(lpath, pattern + mod + "_" + met + "_???.nc")))[0])
                    except: pass
                    else:
                        # models
                        tab2.append(mod)
                        tmp.append(mod)
                        # diagnostic values
                        if metric_collection == "ENSO_tel" and "Map" in met:
                            dict1[mod] = None
                        else:
                            dict1[mod] = data_json[proj][mod]["value"][met]["diagnostic"][mod]["value"]
                if grp == my_project[0] and proj == list_project[0]:
                    if metric_collection == "ENSO_tel" and "Map" in met:
                        tab = dict()
                    else:
                        tab = dict((key1, data_json[proj][list_models[0]]["value"][met]["diagnostic"][key1]["value"])
                                   for key1 in data_json[proj][list_models[0]]["value"][met]["diagnostic"].keys()
                                   if key1 != list_models[0])
                    my_data["diagnostic_values"]["obs"] = tab
                # metric values
                if len(list_models) > 0:
                    if metric_collection == "ENSO_tel" and "Map" in met:
                        for suffix in ["Corr", "Rmse"]:
                            for key1 in data_json[proj][list_models[0]]["value"][met+suffix]["metric"].keys():
                                try:
                                    dict2
                                except:
                                    dict2 = {key1: dict(
                                        (mod, [data_json[proj][mod]["value"][met+suffix]["metric"][key1]["value"]])
                                        for mod in tmp)}
                                else:
                                    try:
                                        dict2[key1]
                                    except:
                                        dict2[key1] = dict(
                                            (mod, [data_json[proj][mod]["value"][met+suffix]["metric"][key1]["value"]]
                                             if key1 in data_json[proj][mod]["value"][met+suffix]["metric"].keys()
                                             else [None])
                                            for mod in tmp)
                                    else:
                                        for mod in tmp:
                                            try:
                                                dict2[key1][mod]
                                            except:
                                                dict2[key1][mod] = [
                                                    data_json[proj][mod]["value"][met+suffix]["metric"][key1]["value"]
                                                if key1 in data_json[proj][mod]["value"][met+suffix]["metric"].keys() else None]
                                            else:
                                                dict2[key1][mod] += \
                                                    [data_json[proj][mod]["value"][met+suffix]["metric"][key1]["value"]
                                                     if key1 in data_json[proj][mod]["value"][met+suffix]["metric"].keys() else None]
                    else:
                        # for key1 in data_json[proj][list_models[0]]["value"][met]["metric"].keys():
                        #     try:
                        #         dict2
                        #     except:
                        #         dict2 = {key1: dict(
                        #             (mod, data_json[proj][mod]["value"][met]["metric"][key1]["value"]
                        #              if key1 in data_json[proj][mod]["value"][met]["metric"].keys() else None)
                        #             for mod in tmp)}
                        #     else:
                        #         try:
                        #             dict2[key1]
                        #         except:
                        #             dict2[key1] = dict(
                        #                 (mod, data_json[proj][mod]["value"][met]["metric"][key1]["value"]
                        #                  if key1 in data_json[proj][mod]["value"][met]["metric"].keys() else None)
                        #                 for mod in tmp)
                        #         else:
                        #             for mod in tmp:
                        #                 dict2[key1][mod] = data_json[proj][mod]["value"][met]["metric"][key1]["value"]\
                        #                     if key1 in data_json[proj][mod]["value"][met]["metric"].keys() else None
                        for key1 in alt_json[proj][list_models[0]]["value"][met]["metric"].keys():
                            try:
                                dict2
                            except:
                                dict2 = {key1: dict(
                                    (mod, alt_json[proj][mod]["value"][met]["metric"][key1]["value"]
                                     if key1 in alt_json[proj][mod]["value"][met]["metric"].keys() else None)
                                    for mod in tmp)}
                            else:
                                try:
                                    dict2[key1]
                                except:
                                    dict2[key1] = dict(
                                        (mod, alt_json[proj][mod]["value"][met]["metric"][key1]["value"]
                                         if key1 in alt_json[proj][mod]["value"][met]["metric"].keys() else None)
                                        for mod in tmp)
                                else:
                                    for mod in tmp:
                                        dict2[key1][mod] =\
                                            alt_json[proj][mod]["value"][met]["metric"][key1]["value"]\
                                                if key1 in alt_json[proj][mod]["value"][met]["metric"].keys() else None
                del lpath, tmp
            my_data["diagnostic_values"][grp] = dict1
            my_data["metric_values"][grp] = dict2
            my_data["filenames"][grp] = tab1
            my_data["models"][grp] = tab2
            print str().ljust(5) + grp + " " + str(len(tab1)).zfill(2) + " file(s) " + \
                  str(len(list_models)).zfill(2) + " name(s)"
    if metric_collection == "ENSO_tel" and "Map" in met:
        units = ""
    else:
        units = data_json[proj][data_json[proj].keys()[0]]["metadata"]["metrics"][met]["diagnostic"]["units"]
    my_data["diagnostic_units"] = units  # .replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    if metric_collection == "ENSO_tel" and "Map" in met:
        units = [data_json[proj][data_json[proj].keys()[0]]["metadata"]["metrics"][met+suffix]["metric"]["units"]
                 for suffix in ["Corr", "Rmse"]]
        # units = [uni.replace("C", "$^\circ$C").replace("long", "$^\circ$long") for uni in units]
    else:
        units = data_json[proj][data_json[proj].keys()[0]]["metadata"]["metrics"][met]["metric"]["units"]
        # units = units.replace("C", "$^\circ$C").replace("long", "$^\circ$long")
    my_data["metric_units"] = units
    # main_plotter(metric_collection, met, my_project, experiment, my_data["filenames"], my_data["diagnostic_values"],
    #              my_data["diagnostic_units"], my_data["metric_values"], my_data["metric_units"], path_png=path_out,
    #              models2=my_data["models"], shading=True)
    main_plotter(metric_collection, met, list_project, experiment, my_data["filenames"], my_data["diagnostic_values"],
                 my_data["diagnostic_units"], my_data["metric_values"], my_data["metric_units"], path_png=path_out,
                 models2=my_data["models"], shading=True)
    del dict1, list_models, my_data, tab1, tab2, units
    try: del dict2
    except: pass

