import pandas as pd


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


print('Reading event')

df_hits = pd.read_csv('event_data/event000001000-hits.csv')
df_truth = pd.read_csv('event_data/event000001000-truth.csv')
df_particles = pd.read_csv('event_data/event000001000-particles.csv')
# add an empty line for particle 0
df_particles.loc[len(df_particles)] = 0

print('Preprocessing')

df = pd.merge(df_hits, df_truth, on='hit_id')
df = pd.merge(df, df_particles, on='particle_id')


to_drop = df[df['volume_id'] > 9].index
df = df.drop(to_drop)
df['global_index'] = calculate_index(df['volume_id'], df['layer_id'])
# Convert the data to cm
df['x'] = df['x'] / 10
df['y'] = df['y'] / 10
df['z'] = df['z'] / 10
df['vx'] = df['vx'] / 10
df['vy'] = df['vy'] / 10
df['vz'] = df['vz'] / 10
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
