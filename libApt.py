__author__ = 'Wilbert'

import math
import numpy as np
import random
import re
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


class Experiment(object):
    def __init__(self):
        self.samples = {}
        self.parameter_map = {}
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
                    # Wild card query
                    if '*.' in linked_parameters[row][col]:
                        suffix = linked_parameters[row][col].split('*.')[-1]
                        for sample_name in self.samples:
                            linked_parameters[row].append('{0}.{1}'.format(self.samples[sample_name].id, suffix))
                        del linked_parameters[row][0]
                    # Attribute query
                    if linked_parameters[row][col][0] == '!':
                        attr_info, suffix = linked_parameters[row][col].split('.')
                        attr_name, attr_value = attr_info[1:].split('=')
                        for sample_name in self.samples:
                            if self.samples[sample_name].metadata[attr_name] == attr_value:
                                linked_parameters[row].append('{0}.{1}'.format(self.samples[sample_name].id, suffix))
                        del linked_parameters[row][0]

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

    def GetNumberUniqueParametersInMap(self):
        return len(set(self.parameter_map.values()))


class DEOpt(object):
    def __init__(self, experiment, min_max=None):
        self.Islands = []
        self.experiment = experiment
        vector_length = experiment.GetNumberUniqueParametersInMap()
        self.min_max = min_max if min_max is not None else np.array([[0., 5.e1] for x in range(vector_length)])
        return

    def CreateIsland(self, population_size):
        island = []
        for i in range(population_size):
            island.append(self.CreateRandomMember())
        self.Islands.append(island)
        return

    def CreateRandomMember(self):
        vector_length = self.experiment.GetNumberUniqueParametersInMap()
        vector = [round(random.uniform(*self.min_max[x,:]),6) for x in range(vector_length)]
        fitness = self.GetFitness(vector)
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
        new_fitness = self.GetFitness(new_vector)
        return Member(new_vector, new_fitness)

    def GetFitness(self, vector):
        # Update model based on vector.
        self.experiment.optimization_vector = vector
        self.experiment.UpdateOptVecParameters()
        # Calculate error.
        error = 0.
        for key,sample in self.experiment.samples.items():
            # Simulate with 0.1 second time steps.
            sim_result = sample.model.handle.simulate(0,300,2999, ['time', 'bound'])
            # Convert concentration to fluorescence
            ci = sample.model.GetValue("ci")
            sim_result = {round(x[0],1): ci*x[1] for x in sim_result}
            # Collate time, observed data, and simulated data.
            data = []
            for index in range(len(sample.measurement)):
                time_point = sample.measurement[index][0]
                if 0 <= time_point < 300:
                    log_obs_point = math.log10(sample.measurement[index][1]) if sample.measurement[index][1] > 0 else 0
                    # Attempt to get the closest simulation point, then log transform that value.
                    sim_point = sim_result.get(time_point)
                    if sim_point is None and time_point > 0:
                        count = 1
                        while sim_point is None:
                            sim_point = sim_result.get(time_point - 0.1 * count)
                            count += 1
                    log_sim_point = math.log10(sim_point) if sim_point > 0 else 0
                    data.append((time_point, log_obs_point, log_sim_point))
            data = np.asarray(data)
            error += sum(abs(data[:,1] - data[:,2])**2)
        return error

    def SortIslandsByFitness(self):
        for island in self.Islands:
            island = sorted(island, key=lambda o: o.Fitness)
        return

class Member(object):
    def __init__(self, vector, fitness):
        self.Vector = vector
        self.Fitness = fitness
        return