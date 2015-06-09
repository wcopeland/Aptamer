__author__ = 'Wilbert'
__version__ = 0.2

from libApt import *
import matplotlib.pyplot as mplot
import random


""" Create system model """
systemOfEquations = """
$dyeExt -> dye; kt * dyeExt
dye -> dyeExt; kt * dye
dye -> ; dil * dye

 -> rna; syn
rna -> ; (dil + degRna) * rna
rna + dye -> bound; kf * rna * dye

bound -> rna + dye; kr * bound
bound -> ; dil * bound
bound -> dye; degBound * bound
"""

""" Specify default parameters """
defaultValues = {'syn':5e-3,        'degRna':1.2e-4,    'degBound':1.2e-4,
                 'kf':1e-1,         'kr':1e-8,          'ci':4000,
                 'dil':4.167e-3,    'kt':1.0,
                 'rna': 0.0,        'bound':0.0,        'dye':0.0,
                 'dyeExt':0.0}

''' LOAD EXPERIMENTAL DATA '''
expt = Experiment()
expt.ImportData("June6\\yaml.txt", default_antimony=systemOfEquations, default_model_values=defaultValues)

''' PLOT ALL EXPERIMENTAL DATA '''
# 5uM Malachite Green & ~0.35 dilution rate
mplot.figure()
for k,v in expt.samples.items():
    if v.metadata["dye_external"] == '5uM' and v.metadata["dilution_rate"] > 0.31:
        label_data = '{0} ({1} hr-1, {2} MG)'.format(v.metadata["strain"], v.metadata["dilution_rate"], v.metadata["dye_external"])
        mplot.plot(v.measurement[:,0], v.measurement[:,1], 'o', label=label_data)
mplot.xlabel("Time (seconds)")
mplot.ylabel("MGA5S Fluorescence (RFUs)")
mplot.legend(loc='upper left')
mplot.show()

# 5uM Malachite Green & ~0.25 dilution rate
mplot.figure()
for k,v in expt.samples.items():
    if v.metadata["dye_external"] == '5uM' and v.metadata["dilution_rate"] < 0.31:
        label_data = '{0} ({1} hr-1, {2} MG)'.format(v.metadata["strain"], v.metadata["dilution_rate"], v.metadata["dye_external"])
        mplot.plot(v.measurement[:,0], v.measurement[:,1], 'o', label=label_data)
mplot.xlabel("Time (seconds)")
mplot.ylabel("MGA5S Fluorescence (RFUs)")
mplot.legend(loc='upper left')
mplot.show()


''' SPECIFY SHARED PARAMETERS BETWEEN SAMPLE MODELS '''
linkedParamList = [ ['!dye_external=5uM.ci'],
                    ['!dye_external=3uM.ci'],
                    ['!dye_external=2uM.ci'],
                    ['!strain=J23100.syn'],
                    ['!strain=J23101.syn'],
                    ['!strain=J23104.syn'],
                    ['!strain=J23107.syn'],
                    ['!strain=J23110.syn'],
                    ['!strain=J23111.syn'],
                    ['*.kf'],
                    ['*.kr'],
                    ['*.kt']]
expt.CreateParameterMap(linked_parameters=linkedParamList)


''' CREATE OPTIMIZATION VECTOR '''
random.seed()
expt.optimization_vector = [round(random.uniform(0,10),1) for x in range(expt.GetNumberUniqueParametersInMap())]

print('\n\n')
expt.UpdateOptVecParameters()

for key in expt.parameter_map:
    print('{0} = {1}'.format(key, expt.GetOptVecParameter(key)))
