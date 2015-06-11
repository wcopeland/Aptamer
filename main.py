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

"""
syn: 1e-6, 1e-1
kf: 1e-8, 1e3
kr:
kt:
ci: 1e1, 1e4

"""

''' LOAD EXPERIMENTAL DATA '''
expt = Experiment()
expt.ImportData("June6\\yaml.txt", default_antimony=systemOfEquations, default_model_values=defaultValues)

''' PLOT ALL EXPERIMENTAL DATA '''

'''# 5uM Malachite Green & ~0.35 dilution rate
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

# 5uM Malachite Green & ~0.25 dilution rate
mplot.figure()
for k,v in expt.samples.items():
    if v.metadata["strain"] == 'J23101' and v.metadata["dilution_rate"] > 0.31:
        label_data = '{0} ({1} hr-1, {2} MG)'.format(v.metadata["strain"], v.metadata["dilution_rate"], v.metadata["dye_external"])
        mplot.plot(v.measurement[:,0], v.measurement[:,1], 'o', label=label_data)
mplot.xlabel("Time (seconds)")
mplot.ylabel("MGA5S Fluorescence (RFUs)")
mplot.legend(loc='upper left')
mplot.show()'''


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





# Differential evolution
GENERATION_COUNT = 0
MAX_GENERATIONS = 200
FITNESS_THRESHOLD = 1e-6

NI = 1      # Number of islands
MF = 0.45   # Migration frequency
NM = 1      # Number of migrants
SP = 3      # Selection policy = Randomly choose one of the top 3 for migration.
RP = 3      # Replacement policy = Randomly choose one of the top 3 for migration.
MT = range(NI)[1:] + [0]    # Migration topology: Ciricular, unidirectional



from datetime import datetime
print('\n\nStarting DE search.\n\n')
clock = datetime.now()

# Initialize Diff Evo routine
de = DEOpt(experiment=expt)

# Create islands
[de.CreateIsland(6) for i in range(NI)]
for island in de.Islands:
    assert len(island) > 3


while True:
    GENERATION_COUNT += 1

    for island in de.Islands:
        samples = [random.sample(island, k=3) for i in range(len(island))]
        trial_values = [de.CreateTrialMember(island[i], samples[i]) for i in range(len(island))]

        print("\n")
        for i in range(len(island)):
            if trial_values[i].Fitness < island[i].Fitness:
                island[i] = trial_values[i]
            print("Member\tFitness: {0}\tVector: {1}".format(island[i].Fitness, island[i].Vector))

        de.SortIslandsByFitness()

    if GENERATION_COUNT % 5 == 0:
        print('\n')
        top_members = [x[0] for x in de.Islands]
        for member in top_members:
            print('Fitness: {0} Vector: {1}'.format(round(member.Fitness,3), member.Vector))

    if GENERATION_COUNT >= MAX_GENERATIONS or min([x[0].Fitness for x in de.Islands]) < FITNESS_THRESHOLD:
        top_members = [x[0] for x in de.Islands]
        for member in top_members:
            print('Fitness: {0} Vector: {1}'.format(round(member.Fitness,3), member.Vector))
        break

print('Done optimizing.')
print('Optimization time: {0}'.format(datetime.now()-clock))
