from pcmdi_metrics.pcmdi.pmp_parser import PMPParser

def ReadOptions():

    P = PMPParser() # Includes all default options

    P.use("--results_dir")
    P.use("--reference_data_path")
    P.use("--modpath")
    
    P.add_argument("--modpath_lf",
                   type=str,
                   dest='modpath_lf',
                   help="Directory path to model land fraction field")
    P.add_argument("--reference_data_lf_path",
                   type=str,
                   dest='reference_data_lf_path',
                   help="Data path to land fraction of reference dataset")
    P.add_argument("--modnames",
                   type=list,
                   default=None,
                   help="List of models")
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
    P.add_argument("-d", "--debug",
                   type=bool,
                   default=False,
                   help="Option for debug: False (defualt) or True")
    
    param = P.get_parameter()

    return param
