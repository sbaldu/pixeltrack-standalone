import optimizer
import numpy as np
import os
import subprocess
from matplotlib import pyplot as plt

num_agents = 20
num_iterations = 200

_this_iteration = 0

lb = [0.0000001, 0.0000001, 0.0000001, 0.0000001, 1.0 / 3.8 / 0.9, 5.0, 400, 
        400, 400, 400, 400, 400, 400, 400, 400, 400, 
        400, 400, 400, 400, 400, 400, 400, 400, 400]

ub = [0.006, 0.03, 0.2, 1.0, 1.0 / 3.8 / 0.3, 20.0, 1000,
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 
        1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]

optimizer.Logger.setLevel('INFO')

# params-good.txt
# 0.00200000009499,0.00300000002608,0.15000000596,0.25,0.0328407224959,522,730,730,522,626,626,522,522,626,626,626,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522,522

def read_csv(filename):
    with open(filename, 'r') as f:
        return list(map(float, f.read().split(',')))

default_params = read_csv('params-good.txt')[:25]

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
    global _this_iteration
    print(f'Iteration {_this_iteration}')
    write_csv('params.txt', x)
    # write_csv(f'optimization/params_{_this_iteration}.txt', x)
    # print(x)
    subprocess.check_call(['./serial', '--maxEvents', '1', '--objective', '--paramsFromFile'])
    objective = get_objectives('objectives.txt')
    # copy ojectives.txt to objectives_{iteration}.txt
    # os.rename('objectives.txt', f'optimization/objectives_{_this_iteration}.txt')
    # print(objective['efficiency'], objective['fake_rate'])
    _this_iteration +=1
    return [1 - objective['efficiency'], objective['fake_rate']]

optimizer.Randomizer.rng = np.random.default_rng(42)

optimizer.FileManager.working_dir = "optimization/"
optimizer.FileManager.loading_enabled = False
optimizer.FileManager.saving_enabled = True

objective = optimizer.ElementWiseObjective([pixeltrack], 2)

pso = optimizer.MOPSO(objective=objective, lower_bounds=lb, upper_bounds=ub,
                      num_particles=num_agents,
                      inertia_weight=0.6, cognitive_coefficient=1, social_coefficient=2, initial_particles_position='random', default_point=default_params, incremental_pareto=True)

# run the optimization algorithm
pso.optimize(num_iterations)

pareto_front = pso.pareto_front
n_pareto_points = len(pareto_front)
pareto_x = [1-particle.fitness[0] for particle in pareto_front]
pareto_y = [particle.fitness[1] for particle in pareto_front]

plt.scatter(pareto_x, pareto_y, s=5)

plt.xlabel('Efficiency')
plt.ylabel('Fake Rate')

plt.savefig('optimization/pf.png')