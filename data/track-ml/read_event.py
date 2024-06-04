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


folder_path ="/eos/user/s/srossiti/track-ml/train_100_events"
output_folder = "/eos/user/s/srossiti/track-ml/preprocessed_train_100_events"
# file names are like: event000001049-cells.csv event000001049-hits.csv event000001049-particles.csv event000001049-truth.csv


print("Reading data")
for event_id in range(100):
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
    # Dropping the columns that are not needed
    df = df.drop(['hit_id', 'volume_id', 'layer_id', 'module_id','px','py','pz','vx','vy', 'tx', 'ty', 'tz', 'tpx', 'tpy', 'tpz', 'weight', 'q'], axis=1)
    print('Saving')
    df.to_csv(f'{output_folder}/trackml_{str(event_id).zfill(3)}.csv', index=False)

print("Done")