__author__ = 'Wilbert'

import math
import matplotlib.pyplot as mplot
import numpy as np
import random
import re
import scipy.stats
from scipy.special import expit
import tellurium as te
import yaml


class Model(object):
    def __init__(self, antimony=None, default_values=None, lockedParameters=None):
        self.antimony = antimony
        self.lockedParameters = [] if lockedParameters == None else lockedParameters
        self.handle = None

        # Create RoadRunner model from antimony and default values.
        if antimony != None:
            if default_values is not None:
                for k,v in default_values.items():
                    antimony = antimony + "\n{0}={1}".format(k,v)
            self.handle = te.loada(antimony)
        return

    def GetParameters(self):
        return {x: self.handle.__getattr__(x) for x in self.handle.getGlobalParameterIds()}

    def GetFloatingSpecies(self):
        return {x: self.handle.__getattr__(x) for x in self.handle.getFloatingSpeciesIds()}

    def GetBoundarySpecies(self):
        return {x: self.handle.__getattr__(x) for x in self.handle.getBoundarySpeciesIds()}

    def GetAll(self):
        list_of_all = self.handle.getGlobalParameterIds() + self.handle.getFloatingSpeciesIds() + self.handle.getBoundarySpeciesIds()
        return {x: self.handle.__getattr__(x) for x in list_of_all}

    def GetValue(self, key):
        try:
            return self.handle.__getattr__(key)
        except:
            return "Value Not Found: {0}".format(key)

    def SetValue(self, key, value):
        try:
            self.handle.__setattr__(key, value)
        except:
            return "Value Not Found: {0}".format(key)

    def LockParameters(self, parameterNames):
        if type(parameterNames) == str:
            parameterNames = [parameterNames, ]

        for parameterName in parameterNames:
            if parameterName in self.lockedParameters:
                print("Warning: {0} has already been fixed.".format(parameterName))
            else:
                self.lockedParameters.append(parameterName)
        return

    def UnlockParameters(self, parameterNames):
        if type(parameterNames) == str:
            parameterNames = [parameterNames, ]

        for parameterName in parameterNames:
            if parameterName in self.lockedParameters:
                self.lockedParameters.remove(parameterName)
            else:
                print("Warning: {0} was not fixed.")
        return

    def GetOptimizableParameters(self):
        result = [x for x in self.GetParameters() if x not in self.lockedParameters]
        return result


class Sample(object):
    ''' '''
    __uid = 0

    def __init__(self, uid=None, model=None, measurement=None):
        Sample.__uid += 1
        self.id = 'Sample{0}'.format(Sample.__uid) if uid == None else uid
        self.description = "No description provided."
        self.metadata = {}
        self.model = model
        self.measurement = measurement

    def Plot(self, new_figure=False, color=None):
        if new_figure:
            mplot.figure()
        c = self.PlotObserved(color=color)
        self.PlotSimulation(color=c)
        return

    def PlotObserved(self, new_figure=False, color=None):
        if new_figure:
            mplot.figure()
        if color is None:
            line = mplot.plot(self.measurement[:,0], self.measurement[:,1], 'o')
            color = line[0].get_color()
        else:
            line = mplot.plot(self.measurement[:,0], self.measurement[:,1], 'o', color)
        mplot.show()
        return color

    def PlotSimulation(self, new_figure=False, color=None):
        if new_figure:
            mplot.figure()
        ci = self.model.GetValue("ci")

        #TODO Hard-coded for now. Think of a better solution for handling species concentrations.
        self.model.SetValue('rna', 0.0)
        self.model.SetValue('bound', 0.0)
        self.model.SetValue('dye', 0.0)
        self.model.SetValue('dyeExt',0.0)
        self.model.handle.steadyState()
        self.model.SetValue('dyeExt', 5.0)
        # !End to do
        self.model.handle.selections = ['time', 'bound']
        MAX_TIME = 180
        simulation = [x for x in self.model.handle.simulate(0,MAX_TIME,MAX_TIME*10-1, stiff=True)]
        result = np.asarray([(x[0],ci*x[1]) for x in simulation])
        if color is None:
            line = mplot.plot(result[:,0], result[:,1])
            color = line[0].get_color()
        else:
            mplot.plot(result[:,0], result[:,1], color=color)
        mplot.show()
        return color


class Experiment(object):
    def __init__(self):
        self.samples = {}
        self.parameter_map = {}
        self.relations = None
        self.optimization_vector = None
        return

    def ImportData(self, yaml_filename, default_antimony=None, default_model_values=None):
        file = open(yaml_filename, 'r')
        data = file.read()
        yaml_data = yaml.load(data)

        for m in yaml_data:
            file = open(yaml_data[m]["filename"], 'r')
            data = file.read()
            data = [[float(y) for y in x.split('\t')] for x in data.split('\n')]

            # Get time points
            raw_time = [x[0] for x in data]

            # Refactor time based on the moment samples were exposed to dye.
            exposure_index = yaml_data[m]["exposure_after"]
            t_exposure = raw_time[exposure_index] - float(yaml_data[m]["time_offset"])
            raw_time = [round(x - t_exposure,1) for x in raw_time]
            read_progression = {int(x):i for i,x in enumerate(yaml_data[m]["read_progression"].split(','))}

            # Get the average of negative control for blank subtraction.
            if yaml_data[m]["samples"].get("NAC") is not None:
                avg_blank_data = []
                meta = yaml_data[m]["samples"].get("NAC")
                sample_columns = [int(x) for x in meta["columns"].split(',')]
                for row in range(len(data)):
                    avg_blank_data.append(sum([data[row][x+1] for x in sample_columns])/len(sample_columns))

            # Calculate and store more precise data points, time points, and sample.
            for name,meta in yaml_data[m]["samples"].items():
                if name != "NAC":
                    sample_data = []
                    sample_columns = [int(x) for x in meta["columns"].split(',')]
                    for row in range(len(data)):
                        for col in sample_columns:
                            # Precise time at which data point was measured.
                            time_point = raw_time[row] + read_progression[col]
                            # Data point. col+2 for time and temp offset. col-1 for zero-based indexing.
                            #  Net = col+1
                            data_point = data[row][col+1] - avg_blank_data[row]
                            sample_data.append((time_point, data_point))

                    # Add t=0 point
                    t0_index = next((i for i in range(len(sample_data)) if sample_data[i][0] > 0))
                    t0_data_estimate = round(sum([x[1] for x in sample_data[t0_index-len(sample_columns):t0_index]])/len(sample_columns),1)
                    sample_data.insert(t0_index, (0.,t0_data_estimate))

                    # Create and store sample.
                    new_model = Model(antimony=default_antimony, default_values=default_model_values)
                    sample = Sample(model=new_model)
                    sample.metadata["strain"] = re.match('(\w+)(\.\d+)*', name).group(1) #name
                    if(meta.get("dilution_rate") is not None):
                        sample.metadata["dilution_rate"] = float(meta["dilution_rate"])
                        # <CAREFUL!> Dilution Rate is converted from hr^-1 to min^-1 in the following line.
                        sample.model.SetValue("dil", float(sample.metadata["dilution_rate"])/60.)
                        sample.model.LockParameters("dil")
                    if(meta.get("dye_external") is not None):
                        sample.metadata["dye_external"] = meta["dye_external"]
                        sample.model.SetValue("dyeExt", float(re.match('\d+',sample.metadata["dye_external"]).group()))
                        sample.model.LockParameters("dyeExt")
                    sample.measurement = np.asarray(sample_data)
                    self.samples[sample.id] = sample
        return

    def AddSamples(self, samples):
        if type(samples) != list:
            samples = [samples]

        for sample in samples:
            self.samples.update({sample.id: sample})

    def GetOptVecParameter(self, key):
        assert len(self.optimization_vector) == (max([x for x in self.parameter_map.values()]) + 1), "ERROR. Parameter Map and Optimization Vector have incompatible dimensions."
        sampleName, paramName = key.split('.')
        return self.samples[sampleName].model.GetValue(paramName)

    def SetOptVecParameter(self, key):
        assert len(self.optimization_vector) == (max([x for x in self.parameter_map.values()]) + 1), "ERROR. Parameter Map and Optimization Vector have incompatible dimensions."
        val = self.optimization_vector[self.parameter_map[key]]
        sampleName, paramName = key.split('.')
        self.samples[sampleName].model.SetValue(paramName, val)
        return

    def UpdateOptVecParameters(self):
        assert len(self.optimization_vector) == (max([x for x in self.parameter_map.values()]) + 1), "ERROR. Parameter Map and Optimization Vector have incompatible dimensions."
        for key in self.parameter_map:
            self.SetOptVecParameter(key)
        return

    def CreateParameterMap(self, linked_parameters=None):
        self.parameter_map = {}
        for sample_name in self.samples:
            for param_name in self.samples[sample_name].model.GetOptimizableParameters():
                param_id = "{0}.{1}".format(sample_name, param_name)
                param_index = len(self.parameter_map)
                self.parameter_map.update({param_id:param_index})

        # Account for linked parameters
        if(linked_parameters != None):
            # Expand parameter selection based on query terms.
            for row in range(len(linked_parameters)):
                for col in range(len(linked_parameters[row])):
                    linked_parameters[row] = self.FindValuesByMeta(linked_parameters[row][col])

        # Update indices of linked parameters and track deleted indices.
        removedIndices = []
        for row in linked_parameters:
            if len(row) == 0:
                continue
            shared_index = self.parameter_map[row[0]]
            for param_name in row[1:]:
                removedIndices.append(self.parameter_map[param_name])
                self.parameter_map[param_name] = shared_index

        # Remove indices spaces.
        for key in self.parameter_map:
            self.parameter_map[key] = self.parameter_map[key] - len([x for x in removedIndices if self.parameter_map[key] > x])
        return


    def FindValuesByMeta(self, query):
        assert query[0] in ['!', '*'], "Not a proper query. {0}".format(query)
        results = []

        meta_query, parameter_handle = query.split('.')
        if query[0] == '*':
            #parameter_handle = query.split('.')
            for sample_name in self.samples:
                results.append("{0}.{1}".format(self.samples[sample_name].id, parameter_handle))
            pass
        elif query[0] == '!':
            #meta_query, parameter_handle = query.split('.')
            meta_name, meta_value = meta_query[1:].split(':')
            for sample_name in self.samples:
                if self.samples[sample_name].metadata[meta_name] == meta_value:
                    results.append('{0}.{1}'.format(self.samples[sample_name].id, parameter_handle))
        return results

    def FindIndicesByValue(self, query):
        result = []
        for key in self.FindValuesByMeta(query):
            result.append(self.parameter_map.get(key))
        return result

    def BuildRelations(self, text):
        self.relations = []
        for line in text:
            target, value, operator = line.split(' ')[:3]
            indices = self.FindIndicesByValue(target)
            if len(set(indices)) > 1:
                print("Warning: Query points to multiple indices! {0}".format(x))
            target = int(indices[0])
            if operator in ['*', '+', '-', '/']:
                reference = line.split(' ')[3]
                value = float(value)
                indices = self.FindIndicesByValue(reference)
                if len(set(indices)) > 1:
                    print("Warning: Query points to multiple indices! {0}".format(x))
                reference = int(indices[0])
                self.relations.append((target, operator, value, reference))
            elif operator in ['~']:
                value = tuple([float(x) for x in value.split(',')])
                self.relations.append((target, operator, value))
            elif operator in ['>', '<']:
                value = float(value)
                self.relations.append((target, operator, value))
        return

    def GetNumberUniqueParametersInMap(self):
        return len(set(self.parameter_map.values()))


class DEOpt(object):
    def __init__(self, experiment, min_max=None):
        self.Islands = []
        self.experiment = experiment
        vector_length = experiment.GetNumberUniqueParametersInMap()
        self.min_max = min_max if min_max is not None else np.array([[0., 1.e2] for x in range(vector_length)])
        self.rules = None
        return

    def CreateIsland(self, population_size, min_max_rules=None):
        island = []
        for i in range(population_size):
            island.append(self.CreateRandomMember(min_max_rules))
        self.Islands.append(island)
        return

    def LogRandom(self, lower, upper):
        return 10 ** random.uniform(lower, upper)


    def CreateRandomMember(self, min_max_rules=None):
        vector_length = self.experiment.GetNumberUniqueParametersInMap()
        success = False

        while not success:
            try:
                vector = [random.uniform(0,1e3) for x in range(vector_length)]

                if min_max_rules is not None:
                    for entry in min_max_rules:
                        opt_vec_indices = set(self.experiment.FindIndicesByValue(entry))
                        for opt_vec_index in opt_vec_indices:
                            vector[opt_vec_index] = self.LogRandom(*min_max_rules[entry])

                fitness = self.GetFitness(vector)
                success = True
            except:
                pass
        print("RandomMember\tFitness: {0}\tVector: {1}".format(fitness, vector))
        return Member(vector, fitness)

    def CreateTrialMember(self, original, samples, CR=0.6, F=0.8):
        o = original.Vector
        a = samples[0].Vector
        b = samples[1].Vector
        c = samples[2].Vector

        new_vector = []
        for i in range(len(o)):
            if random.random() <= CR:
                v = round(a[i] + F * (b[i] - c[i]), 7)
                if v>0:
                    new_vector.append(v)
                else:
                    new_vector.append(o[i]/2.)
            else:
                new_vector.append(o[i])
        try:
            new_fitness = self.GetFitness(new_vector)
            return Member(new_vector, new_fitness)
        except:
            return original

    def EnforceRelations(self, vector):
        for relation in self.experiment.relations:
            target_index = relation[0]
            operation = relation[1]
            if operation == '*':
                value, reference_index = relation[2:]
                vector[target_index] = vector[reference_index] * value
            elif operation == '/':
                value, reference_index = relation[2:]
                vector[target_index] = vector[reference_index] / value
            elif operation == '+':
                value, reference_index = relation[2:]
                vector[target_index] = vector[reference_index] + value
            elif operation == '-':
                value, reference_index = relation[2:]
                vector[target_index] = vector[reference_index] - value
            elif operation == '~':
                min_val, max_val = relation[2]
                if vector[target_index] < min_val:
                    vector[target_index] = min_val
                elif vector[target_index] > max_val:
                    vector[target_index] = max_val
            elif operation == '<':
                value = relation[2]
                if vector[target_index] > value:
                    vector[target_index] = value
            elif operation == '>':
                value = relation[2]
                if vector[target_index] < value:
                    vector[target_index] = value
        return

    def GetFitness(self, vector):
        # Update model based on vector.
        self.EnforceRelations(vector)
        self.experiment.optimization_vector = vector
        self.experiment.UpdateOptVecParameters()
        # Calculate error.
        error = 0.
        for key,sample in self.experiment.samples.items():
            #TODO Clean up code that "smartly" assigns species concentrations.
            MAX_TIME = 300

            # Simulate steady-state in absence of dye
            sample.model.SetValue('rna', 0.0)
            sample.model.SetValue('bound', 0.0)
            sample.model.SetValue('dye', 0.0)
            sample.model.SetValue('dyeExt', 0.0)
            sample.model.handle.steadyState()
            # Simulate with 0.1 second time steps.
            sample.model.SetValue('dyeExt', 5.0)
            sample.model.handle.selections = ['time', 'bound']
            sim_result = sample.model.handle.simulate(0,MAX_TIME, MAX_TIME*10-1, stiff=True)
            # Convert concentration to fluorescence
            ci = sample.model.GetValue("ci")
            sim_result = {round(x[0],1): ci*x[1] for x in sim_result}
            # Collate time, observed data, and simulated data.
            data = []
            for index in range(len(sample.measurement)):
                time_point = sample.measurement[index][0]
                if 0 <= time_point < MAX_TIME:
                    #log_obs_point = math.log(sample.measurement[index][1]) if sample.measurement[index][1] > 0 else 0
                    obs_point = sample.measurement[index][1]
                    # Attempt to get the closest simulation point, then log transform that value.
                    sim_point = sim_result.get(time_point)
                    if sim_point is None and time_point > 0:
                        count = 1
                        while sim_point is None:
                            sim_point = sim_result.get(time_point - 0.1 * count)
                            count += 1
                    #log_sim_point = math.log(sim_point) if sim_point > 0 else 0
                    data.append((time_point, obs_point, sim_point))
            data = np.asarray(data)
            diff_squared = abs(data[:,1] - data[:,2])**2

            #TODO Identify the best weighting scheme.
            # SCHEME 1:  NO WEIGHTING
            weights = np.ones(data.shape)
            # SCHEME 2: WEIGHTING WITH SIGMOID
            weight_lower_threshold = 0.1
            steepness=15.
            hshift=10.
            weights = expit(-data[:,0]/steepness + hshift)
            weights[weights<weight_lower_threshold] = weight_lower_threshold
            # SCHEME 3: WEIGHTING WITH EXPONENTIAL
            #pdf_scalar = 4./300.
            #weights = scipy.stats.expon.pdf(pdf_scalar*data[:,0])

            # Sum error
            error += sum(diff_squared * weights)
        return error

    def SortIslandsByFitness(self):
        for island in self.Islands:
            island = sorted(island, key=lambda o: o.Fitness)
        return

    def PlotResults(self, member):
        self.experiment.optimization_vector = member.Vector
        self.experiment.UpdateOptVecParameters()
        mplot.figure()
        for key,sample in self.experiment.samples.items():
            color = sample.PlotObserved()
            sample.PlotSimulation(color=color)
        return

    def GenerateReport(self):
        report = dict()
        # Save Parameter Map
        report.update({"parameter_map": self.experiment.parameter_map})
        # Save Members
        report.update({"islands": []})
        for island_index in range(len(self.Islands)):
            population = []
            for member in self.Islands[island_index]:
                population.append((member.Fitness, member.Vector))
            report["islands"].append(population)
        return report

class Member(object):
    def __init__(self, vector, fitness):
        self.Vector = vector
        self.Fitness = fitness
        return