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
                 'kf':5e-1,         'kr':5e-8,          'ci':3800,
                 'dil':4.167e-3,    'kt':3.5e-3,
                 'rna': 0.0,        'bound':0.0,        'dye':0.0,
                 'dyeExt':0.0}

''' LOAD EXPERIMENTAL DATA '''
expt = Experiment()
expt.ImportData("June6\\yaml.txt", default_antimony=systemOfEquations, default_model_values=defaultValues)
for key in expt.samples:
    expt.samples[key].model.LockParameters(["degRna", "degBound"])

''' PLOT ALL EXPERIMENTAL DATA '''

# 5uM Malachite Green & ~0.35 dilution rate
mplot.figure()
for k,v in expt.samples.items():
    if v.metadata["dye_external"] == '5uM' and v.metadata["dilution_rate"] > 0.31:
        label_data = '{0} ({1} hr-1, {2} MG)'.format(v.metadata["strain"], v.metadata["dilution_rate"], v.metadata["dye_external"])
        mplot.plot(v.measurement[:,0], v.measurement[:,1], 'o', label=label_data)
        fl_val = v.measurement[:,1][v.measurement[:,0] > 150][0]
        print('{0} {1}: {2}'.format(v.metadata['strain'], v.metadata['dilution_rate'], fl_val))
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

# 5uM Malachite Green & ~0.25 dilution rate
mplot.figure()
for k,v in expt.samples.items():
    if v.metadata["strain"] == 'J23101' and v.metadata["dilution_rate"] > 0.31:
        label_data = '{0} ({1} hr-1, {2} MG)'.format(v.metadata["strain"], v.metadata["dilution_rate"], v.metadata["dye_external"])
        mplot.plot(v.measurement[:,0], v.measurement[:,1], 'o', label=label_data)
mplot.xlabel("Time (seconds)")
mplot.ylabel("MGA5S Fluorescence (RFUs)")
mplot.legend(loc='upper left')
mplot.show()


''' SPECIFY SHARED PARAMETERS BETWEEN SAMPLE MODELS '''
linkedParamList = [ ['!dye_external:5uM.ci'],
                    ['!dye_external:3uM.ci'],
                    ['!dye_external:2uM.ci'],
                    ['!strain:J23100.syn'],
                    ['!strain:J23101.syn'],
                    ['!strain:J23104.syn'],
                    ['!strain:J23107.syn'],
                    ['!strain:J23110.syn'],
                    ['!strain:J23111.syn'],
                    ['*.kf'],
                    ['*.kr'],
                    ['*.kt']]

relations = ["*.kf 1 <",
             "*.kt 0.1 <",
             "*.kr 1e-7 * *.kf",
             "!strain:J23100.syn 1.23 * !strain:J23101.syn",
             "!strain:J23104.syn 1.99 * !strain:J23101.syn",
             "!strain:J23107.syn 0.3 * !strain:J23101.syn",
             "!strain:J23110.syn 0.63 * !strain:J23101.syn",
             "!strain:J23111.syn 2.36 * !strain:J23101.syn", ]

random_member_min_max = { '*.kf':(-9, 1),
                          '*.kr':(-9, 1),
                          '*.kt':(-7, 1),
                          '*.syn':(-5, 0),
                          '*.ci':(1,4)}

expt.CreateParameterMap(linked_parameters=linkedParamList)
expt.BuildRelations(relations)


"""
# Manually fitting Sample4 (J23111)
expt.samples['Sample4'].model.SetValue('syn', 4.8e-3)
expt.samples['Sample4'].PlotObserved()

#Sweep kt
kt = expt.samples['Sample4'].model.GetValue('kt')
for x in [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20]:
    expt.samples['Sample4'].model.SetValue('kt', (x*kt))
    expt.samples['Sample4'].PlotSimulation()
expt.samples['Sample4'].model.SetValue('kt', kt)

#Sweep kf
mplot.figure()
expt.samples['Sample4'].PlotObserved()
kf = expt.samples['Sample4'].model.GetValue('kf')
for x in [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20]:
    expt.samples['Sample4'].model.SetValue('kf', (x*kf))
    expt.samples['Sample4'].PlotSimulation()
    print('Ration kf/kt = {0}'.format((x*kf/kt)))
expt.samples['Sample4'].model.SetValue('kf', kf)
"""

