import optimizer
import numpy as np
import os
import subprocess
from matplotlib import pyplot as plt

num_agents = 10
num_iterations = 50


lb = [0.000, 0.00, 0.0, 0.0, 1.0 / 3.8 / 0.9, 0, 400, 
        400, 400, 400, 400, 400, 400, 400, 400, 400, 
        400, 400, 400, 400, 400, 400, 400, 400, 400,
        400, 400]

ub = [0.003, 0.01, 0.2, 1.0, 1.0 / 3.8 / 0.3, 1, 1000,
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
        1000, 1000]

optimizer.Logger.setLevel('DEBUG')

# params-good.txt
# 0.00200000009499,0.00300000002608,0.15000000596,0.25,0.0328407224959,522,730,730,522,626,626,522,522,626,626,626,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522
# 0.00200000009499,0.00300000002608,0.15000000596,0.25,0.0328407224959,1,522,730,730,522,626,626,522,522,626,626,626,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522

def read_csv(filename):
    with open(filename, 'r') as f:
        return list(map(float, f.read().split(',')))

default_params = read_csv('params-good.txt')[:27]

def write_csv(filename, x):
    with open(filename, 'w') as f:
        f.write(','.join(map(str, x)))
    
# $ cat objectives.txt
# efficiency 0.109467
# fake_rate 0.968843
# fakes 31748
# matching 1021
# reconstructed 32769
# simulated 9327

def get_objectives(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        objectives = {}
        for line in lines:
            key, value = line.split()
            if value == '-nan':
                value = 1.0
            objectives[key] = float(value)
        return objectives
    

# $ ./serial --maxEvents 1 --objective --paramsFromFile
def pixeltrack(x):
    write_csv('params.txt', x)
    # write_csv(f'optimization/params_{_this_iteration}.txt', x)
    # print(x)
    try:
        subprocess.check_call(['./serial', '--maxEvents', '100', '--numberOfThreads', '16','--objective', '--paramsFromFile'])
        objective = get_objectives('objectives.txt')
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        objective = {"efficiency":0, "fake_rate":1}
    # copy ojectives.txt to objectives_{iteration}.txt
    # os.rename('objectives.txt', f'optimization/objectives_{_this_iteration}.txt')
    # print(objective['efficiency'], objective['fake_rate'])
    return [1 - objective['efficiency'], objective['fake_rate']]

optimizer.Randomizer.rng = np.random.default_rng(42)

optimizer.FileManager.working_dir = "optimization/"
optimizer.FileManager.loading_enabled = False
optimizer.FileManager.saving_enabled = True
optimizer.FileManager.saving_csv_enabled = True
optimizer.FileManager.headers_enabled = True

objective = optimizer.ElementWiseObjective([pixeltrack], 2, objective_names=['1-efficiency', 'fake_rate'])

param_names = [
    "CAThetaCutBarrel",
    "CAThetaCutForward",
    "dcaCutInnerTriplet",
    "dcaCutOuterTriplet",
    "hardCurvCut",
    "doZ0Cut",
    "phiCut_0", "phiCut_1", "phiCut_2", "phiCut_3", "phiCut_4", "phiCut_5", "phiCut_6", "phiCut_7",
    "phiCut_8", "phiCut_9", "phiCut_10", "phiCut_11", "phiCut_12", "phiCut_13", "phiCut_14", "phiCut_15",
    "phiCut_16", "phiCut_17", "phiCut_18", "phiCut_19", "phiCut_20"
]

pso = optimizer.MOPSO(objective=objective, lower_bounds=lb, upper_bounds=ub,
                      num_particles=num_agents,
                      inertia_weight=0.6, cognitive_coefficient=1, social_coefficient=2,
                      initial_particles_position='gaussian', default_point=default_params,
                      topology='higher_weighted_crowding_distance',
                      param_names=param_names)

# run the optimization algorithm
pso.optimize(num_iterations)
