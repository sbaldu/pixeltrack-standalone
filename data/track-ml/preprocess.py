
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


# data = pd.read_csv('event000001000-hits.csv')
# data = data[data['volume_id'] <= 9]
# data['global_index'] = calculate_index(data['volume_id'], data['layer_id'])
# print(data.to_markdown())
# data = data.drop(['hit_id', 'volume_id', 'layer_id', 'module_id'], axis=1)
# data.to_csv('hits_1000.csv', index=False)

data = pd.read_csv('event000001000-truth.csv')
data = data.drop(['hit_id'], axis=1)
data.to_csv('truth_1000.csv', index=False)
