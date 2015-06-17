__author__ = 'Wilbert'
__version__ = 0.2

from datetime import datetime
from libApt import *
import matplotlib.pyplot as mplot
import random



""" CREATE A SYSTEMS MODEL """
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

""" SPECIFY DEFAULT VALUES FOR PARAMETERS AND SPECIES """
defaultValues = {'syn':5e-3,        'degRna':1.2e-4,    'degBound':1.2e-4,
                 'kf':5e-1,         'kr':5e-8,          'ci':3800,
                 'dil':4.167e-3,    'kt':3.5e-3,
                 'rna': 0.0,        'bound':0.0,        'dye':0.0,
                 'dyeExt':0.0}

""" LOAD EXPERIMENTAL DATA """
expt = Experiment()
expt.ImportData("June6\\yaml.txt", default_antimony=systemOfEquations, default_model_values=defaultValues)
for key in expt.samples:
    expt.samples[key].model.LockParameters(["degRna", "degBound"])


""" SPECIFY PARAMETER CONNECTIONS AND CONSTRAINTS """
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
             "*.kt 1 <",
             "*.kr 1.17e-7 * *.kf",]
             #"!strain:J23100.syn 1.4 * !strain:J23101.syn",
             #"!strain:J23104.syn 2.2 * !strain:J23101.syn",
             #"!strain:J23107.syn 0.3 * !strain:J23101.syn",
             #"!strain:J23110.syn 0.69 * !strain:J23101.syn",
             #"!strain:J23111.syn 2.6 * !strain:J23101.syn", ]

random_member_min_max = { '*.kf':(-9, 1),
                          '*.kr':(-9, 1),
                          '*.kt':(-7, 1),
                          '*.syn':(-5, 0),
                          '*.ci':(2,5)}

expt.CreateParameterMap(linked_parameters=linkedParamList)
expt.BuildRelations(relations)


# Differential evolution
GENERATION_COUNT = 0
MAX_GENERATIONS = 1000
FITNESS_THRESHOLD = 1e-4

NI = 2      # Number of islands
NP = 30      # Number of members per island
MF = 0.45   # Migration frequency
NM = 1      # Number of migrants
SP = 3      # Selection policy = Randomly choose one of the top 3 for migration.
RP = 3      # Replacement policy = Randomly choose one of the top 3 for migration.
MT = range(NI)[1:] + [0]    # Migration topology: Ciricular, unidirectional





random.seed()
print('\n\nStarting DE search.\n\n')
clock = datetime.now()

# Initialize Diff Evo routine
de = DEOpt(experiment=expt)

# Create islands
[de.CreateIsland(NP, random_member_min_max) for i in range(NI)]
for island in de.Islands:
    assert len(island) > 3


while True:
    GENERATION_COUNT += 1

    for island in de.Islands:
        samples = [random.sample(island, k=3) for i in range(len(island))]
        trial_values = [de.CreateTrialMember(island[i], samples[i]) for i in range(len(island))]

        for i in range(len(island)):
            if trial_values[i].Fitness < island[i].Fitness:
                island[i] = trial_values[i]

        de.SortIslandsByFitness()

    if GENERATION_COUNT % 10 == 0:
        print('\n')
        #top_members = [x[0] for x in de.Islands]
        for island in de.Islands:
            for member in island:
                print('Fitness: {0}\tVector: {1}'.format(round(member.Fitness,3), member.Vector))

    if GENERATION_COUNT >= MAX_GENERATIONS or min([x[0].Fitness for x in de.Islands]) < FITNESS_THRESHOLD:
        top_members = [x[0] for x in de.Islands]
        for member in top_members:
            print('Fitness: {0} Vector: {1}'.format(round(member.Fitness,3), member.Vector))
        break

print('Done optimizing.')
print('Optimization time: {0}'.format(datetime.now()-clock))

de.PlotResults(de.Islands[0][0])
report = de.GenerateReport()



