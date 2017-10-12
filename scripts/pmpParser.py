from pcmdi_metrics.pcmdi.pmp_parser import PMPParser

def ReadOptions():

    P = PMPParser() # Includes all default options
    
    P.add_argument("-mp", "--modpath",
                   type=str,
                   dest='modpath',
                   help="Directory path to model monthly field")
    P.add_argument("-op", "--obspath",
                   type=str,
                   dest='obspath',
                   help="Directory path to obs monthly field")
    P.add_argument("-mns", "--modnames",
                   type=str,
                   dest='modnames',
                   help="List of models to apply")
    P.add_argument("-varobs", "--variableobs",
                   type=str,
                   dest='variableobs',
                   help="Variable name in observation")
    P.add_argument("-outpj", "--outpathjsons",
                   type=str,
                   dest='outpathjson',
                   help="Output path for json")
    P.add_argument("-outnj", "--outnamejson",
                   type=str,
                   dest='outnamejson',
                   help="Output file name for json")
    P.add_argument("-outpd", "--outpathdata",
                   type=str,
                   dest='outpathdata',
                   help="Output path for data")
    P.add_argument("-mx", "--metrics",
                   type=str,
                   dest='metrics',
                   help="List of metrics")
    
    P.add_argument("--metricsCollection",
                   type=str,
                   dest='metricsCollection',
                   help="Metrics Collection e.g. MC1, MC2, etc")

    P.add_argument("--sstName",
                   type=str,
                   dest='sstName',
                   help="Variable name for SST in the model")
    P.add_argument("--tauxName",
                   type=str,
                   dest='tauxName',
                   help="Variable name for taux in the model")
    
    P.add_argument("--sstNameObs",
                   type=str,
                   dest='sstName',
                   help="Variable name for SST in the observation")
    P.add_argument("--tauxNameObs",
                   type=str,
                   dest='tauxName',
                   help="Variable name for taux in the observation")
    
    P.add_argument("--sstObsPath",
                   type=str,
                   dest='sstObsPath',
                   help="Directory path to obs monthly SST field")
    P.add_argument("--tauuObsPath",
                   type=str,
                   dest='tauuObsPath',
                   help="Directory path to obs monthly tauu field")
    
    param = P.get_parameter()

    return param
