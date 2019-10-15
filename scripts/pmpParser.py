from pcmdi_metrics.driver.pmp_parser import PMPParser
import copy
import re


def ReadOptions():

    P = PMPParser() # Includes all default options

    #P.use("--mip")
    #P.use("--exp")
    
    P.add_argument("--mip",
                   type=str,
                   default="cmip5",
                   help="A WCRP MIP project such as CMIP3 and CMIP5")
    P.add_argument("--exp",
                   type=str,
                   default="historical",
                   help="An experiment such as AMIP, historical or pi-contorl")
    P.use("--modpath")
    P.add_argument("--modpath_lf",
                   type=str,
                   dest='modpath_lf',
                   help="Directory path to model land fraction field")
    P.add_argument("--modnames",
                   type=str,
                   nargs='+',
                   default=None,
                   help="List of models")
    P.add_argument("-r", "--realization",
                   type=str,
                   default="r1i1p1",
                   help="Consider all accessible realizations as idividual\n"
                        "- r1i1p1: default, consider only 'r1i1p1' member\n"
                        "          Or, specify realization, e.g, r3i1p1'\n"
                        "- *: consider all available realizations")
    P.use("--reference_data_path")
    P.add_argument("--reference_data_lf_path",
                   type=str,
                   dest='reference_data_lf_path',
                   help="Data path to land fraction of reference dataset")
    P.add_argument("--metricsCollection",
                   type=str,
                   dest='metricsCollection',
                   help="Metrics Collection e.g. MC1, ENSO_perf, or ENSO_tel")
    P.add_argument("--json_name",
                   type=str,
                   dest='json_name',
                   help="File name for output JSON")
    P.add_argument("--netcdf_name",
                   type=str,
                   dest='netcdf_name',
                   help="File name for output NetCDF")
    P.use("--results_dir")

    # Switches
    P.add_argument("-d", "--debug", nargs='?',
                   const=True, default=False,
                   type=bool,
                   help="Option for debug: True / False (defualt)")
    P.add_argument("--nc_out", nargs='?',
                   const=True, default=True,
                   type=bool,
                   help="Option for generate netCDF file output: True (default) / False")
    
    param = P.get_parameter()

    return param


def sort_human(input_list):
    tmp_list = copy.copy(input_list)
    convert = lambda text: float(text) if text.isdigit() else text
    alphanum = lambda key: [convert(c) for c in re.split('([-+]?[0-9]*\.?[0-9]*)', key)]
    tmp_list.sort(key=alphanum)
    return tmp_list
