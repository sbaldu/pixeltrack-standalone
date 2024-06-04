
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


df_hits = pd.read_csv('event000001000-hits.csv')
# data = data[data['volume_id'] <= 9]
to_drop = df_hits[df_hits['volume_id'] > 9].index
df_hits = df_hits.drop(to_drop) 
df_hits['global_index'] = calculate_index(df_hits['volume_id'], df_hits['layer_id'])
# print(df_hits.to_markdown())
df_hits = df_hits.drop(['volume_id', 'layer_id', 'module_id'], axis=1)
# df_hits.to_csv('hits_1000.csv', index=False)

# 'event000001000-particles.csv'
# particle_id,vx,vy,vz,px,py,pz,q,nhits
# 4503668346847232,-0.00928816,0.00986098,-0.0778789,-0.0552689,0.323272,-0.203492,-1,8
# 4503737066323968,-0.00928816,0.00986098,-0.0778789,-0.948125,0.470892,2.01006,1,11
# 4503805785800704,-0.00928816,0.00986098,-0.0778789,-0.886484,0.105749,0.683881,-1,0

df_particles = pd.read_csv('event000001000-particles.csv')
# compute pT as sqrt(px^2 + py^2)
df_particles['pt'] = (df_particles['px'] ** 2 + df_particles['py'] ** 2) ** 0.5


df_truth = pd.read_csv('event000001000-truth.csv')
df_truth =df_truth.drop(to_drop)
# # assign the correct pt based on the particle id
# df_truth['pt'] = df_truth['particle_id'].map(df_particles.set_index('particle_id')['pt'])
# # if particle_id is 0 set pt to 1
# df_truth['pt'] = df_truth['pt'].fillna(0)
# # df_truth.to_csv('truth_1000.csv', index=False)

# Merge all datasets based on hit_id
df = df_hits.merge(df_truth, on='hit_id')

# Merge all datasets based on particle_id
df = df.merge(df_particles, on='particle_id')

# Fill particle_id, vx, vy, vz, px, py, pz, q, nhits, pt with 0 if NaN
df = df.fillna(0)

# print the first 5 rows of the dataframe
print(df.head())

