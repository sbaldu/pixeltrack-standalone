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

folder_path ="/eos/user/s/srossiti/track-ml/train_100_events"
output_folder = "/eos/user/s/srossiti/track-ml/preprocessed_train_100_events"
# file names are like: event000001049-cells.csv event000001049-hits.csv event000001049-particles.csv event000001049-truth.csv


print("Reading data")
sumdf = pd.DataFrame()
for event_id in range(10):
    print(f"Reading event {event_id}")
    df_hits = pd.read_csv(f'{folder_path}/event000001{str(event_id).zfill(3)}-hits.csv')
    df_particles = pd.read_csv(f'{folder_path}/event000001{str(event_id).zfill(3)}-particles.csv')
    df_truth = pd.read_csv(f'{folder_path}/event000001{str(event_id).zfill(3)}-truth.csv')
    print("Preprocessing")
    df_particles.loc[len(df_particles)] = 0

    df = pd.merge(df_hits, df_truth, on='hit_id')
    df = pd.merge(df, df_particles, on='particle_id')
    # Dropping the hits not in the barrel or endcaps
    to_drop = df[df['volume_id'] > 9].index
    df = df.drop(to_drop)
    # Calculate the layer index
    df['global_index'] = calculate_index(df['volume_id'], df['layer_id'])
    df['det_id'] = calculate_module_id(df['global_index'], df['module_id'])
    # Convert the data to cm
    df['x'] = df['x'] / 10
    df['y'] = df['y'] / 10
    df['z'] = df['z'] / 10
    df['phi'] = np.arctan2(df['y'], df['x'])

    df_barrel_0 = df[df['global_index'].isin([0])]
    df_barrel_0 = df_barrel_0[(df_barrel_0['z'] > 10) & (df_barrel_0['z'] < 10.5)]
    # add df_barrel_0[['x', 'y', 'z', 'det_id']] to sumdf
    sumdf = pd.concat([sumdf, df_barrel_0[['x', 'y', 'z','phi', 'det_id']]])

# print mean std min max and delta phi for each det_id
report_df = sumdf.groupby('det_id').agg({'phi': ['mean', 'std', 'min', 'max']}).reset_index()
report_df['delta_phi'] = report_df['phi']['max'] - report_df['phi']['min']
print(report_df)


fig = plt.figure()
plt.scatter(sumdf['x'], sumdf['y'], c=sumdf['det_id'], cmap='tab20', s=0.1)
# show colorbar
plt.savefig('scatter2d.png')
