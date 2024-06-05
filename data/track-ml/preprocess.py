import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calculate_index(volumes: pd.Series, layers: pd.Series) -> list:
    ids = []

    for v, l in zip(volumes, layers):
        id = 0
        if v == 7:
            id += 3
            id += 8 - l / 2
        elif v == 9:
            id += 10
            id += l / 2
        else:
            id += l / 2 - 1
        ids.append(int(id))

    return ids

def calculate_module_id(global_index: pd.Series, module_id: pd.Series) -> list:
    max_module_id = [224, 448, 728, 1092, 108, 108,
                     108, 108, 108, 108, 108, 108,
                     108, 108, 108, 108, 108, 108]
    start_id = []
    for i in range(len(max_module_id)):
        start_id.append(sum(max_module_id[:i]))
        
    detId = []
    for g, m in zip(global_index, module_id):
        detId.append(m + start_id[g])
    return detId

print('Reading event')

df_hits = pd.read_csv('/eos/user/s/srossiti/track-ml/train_100_events/event000001010-hits.csv')
df_truth = pd.read_csv('/eos/user/s/srossiti/track-ml/train_100_events/event000001010-truth.csv')
df_particles = pd.read_csv('/eos/user/s/srossiti/track-ml/train_100_events/event000001010-particles.csv')
# add an empty line for particle 0
df_particles.loc[len(df_particles)] = 0

print('Preprocessing')

df = pd.merge(df_hits, df_truth, on='hit_id')
df = pd.merge(df, df_particles, on='particle_id')


to_drop = df[df['volume_id'] > 9].index
df = df.drop(to_drop)
df['global_index'] = calculate_index(df['volume_id'], df['layer_id'])
df['new_module_id'] = calculate_module_id(df['global_index'], df['module_id'])
# Convert the data to cm
df['x'] = df['x'] / 10
df['y'] = df['y'] / 10
df['z'] = df['z'] / 10
df['vx'] = df['vx'] / 10
df['vy'] = df['vy'] / 10
df['vz'] = df['vz'] / 10

# constexpr short phi2short(float x) {
#   constexpr float p2i = ((int)(std::numeric_limits<short>::max()) + 1) / M_PI;
#   return std::round(x * p2i);
# }

# def phi2short(x):
#     p2i = ((int)(np.iinfo(np.int16).max) + 1) / np.pi
#     return np.round(x * p2i)

df['phi'] = np.arctan2(df['y'], df['x'])
# print(df[df['global_index']==0][['new_module_id','z', 'phi']].sort_values(by='new_module_id').to_string())


# count the number of unique module_id grouped by phi in global_index 0
df_barrel_0 = df[df['global_index'].isin([0])]
print(df[df['global_index'].isin([0])].groupby('new_module_id')['phi'].mean().to_string())
exit()

#select those with z between 20 and 21
df_barrel_0 = df_barrel_0[(df_barrel_0['z'] > 0) & (df_barrel_0['z'] < 41)]


fig = plt.figure()
plt.scatter(df_barrel_0['x'], df_barrel_0['y'], c=df_barrel_0['new_module_id'], cmap='rainbow')
# show colorbar
plt.savefig('scatter2d.png')


exit()

# Compute pT and dR
df['pt'] = (df['px'] ** 2 + df['py'] ** 2) ** 0.5
df['dr'] = (df['vx'] ** 2 + df['vy'] ** 2) ** 0.5
# Compute nhits
df['nhits'] = df['particle_id'].map(df['particle_id'].value_counts())
# Compute the number of unique global_indexes per particle
df['nlayers'] = df.groupby('particle_id')['global_index'].transform('nunique')
df = df.drop(['hit_id', 'volume_id', 'layer_id', 'module_id','px','py','pz','vx','vy', 'tx', 'ty', 'tz', 'tpx', 'tpy', 'tpz', 'weight', 'q'], axis=1)

print('Saving')

df.to_csv('trackml_1000.csv', index=False)
