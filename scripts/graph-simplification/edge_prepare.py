import sys
import pandas as pd
import numpy as np
import torch
from sklearn.model_selection import train_test_split
from torch.utils.data import TensorDataset

edge_name_list = ["ID","Trans_ID1","Trans_ID2","HET","HOM","AD","SN","Mini","Homopolymer"] + ["Graphlet" + str(i) for i in range(0,68)] + ["training_y","testing_y","316_y","318_y"]

edge_input = sys.argv[1]
edge_output_statistics = sys.argv[2]
edge_output_dataset = sys.argv[3]

val_ratio = 0.05
test_ratio = 0

#Generate the training, validation, testing set
data = pd.read_csv(edge_input, sep='\t', header=None, names=edge_name_list, na_values=['NONE','nan'])
train_id_df, val_id_df = train_test_split(data[data['training_y'].notna()]["ID"], test_size=(val_ratio + test_ratio), random_state=42)
data['train_mask'] = data.index.isin(train_id_df.index)
data['val_mask'] = data.index.isin(val_id_df.index)
data['test_mask'] = data.index.isin(data[data['testing_y'].notna()].index)


#Calculate the statistics of features
fo = open(edge_output_statistics,"w")
feature_dict = {}
for feature in ["HET","HOM","AD","SN"]:
    feature_stat1 = data[feature].min()
    feature_stat2 = data[feature].max()
    feature_stat3 = data[feature].median()
    feature_dict[feature] = [feature_stat1, feature_stat2, feature_stat3]
    fo.write(feature + '\t' + str(feature_stat1) + '\t' + str(feature_stat2) + '\t' + str(feature_stat3) + '\n')

for feature in ["Graphlet" + str(i) for i in range(0,68)]:
    feature_stat1 = data[feature].quantile(0.25)
    feature_stat2 = data[feature].quantile(0.75)
    feature_stat3 = data[feature].median()
    feature_dict[feature] = [feature_stat1, feature_stat2, feature_stat3]
    fo.write(feature + '\t' + str(feature_stat1) + '\t' + str(feature_stat2) + '\t' + str(feature_stat3) + '\n')
fo.close()


#Fill the missing value with median
for feature in ["Graphlet" + str(i) for i in range(0,68)]:
    data[feature] = data[feature].fillna(feature_dict[feature][2])


#Normalization
for feature in ["HET","HOM","AD","SN"] + ["Graphlet" + str(i) for i in range(0,68)]:
    if feature_dict[feature][1] - feature_dict[feature][0] == 0:
        continue
    else:
        data[feature] = (data[feature] - feature_dict[feature][0])/(feature_dict[feature][1] - feature_dict[feature][0])


#Create Dataset
X = torch.tensor(data[["HET","HOM","AD","SN","Mini","Homopolymer"] + ["Graphlet" + str(i) for i in range(0,68)]].to_numpy(), dtype=torch.float32)
y = torch.tensor(data[['training_y', 'testing_y']].stack().to_list(), dtype=torch.float32)
train_mask = torch.tensor(data['train_mask'], dtype=torch.bool).squeeze()
val_mask = torch.tensor(data['val_mask'], dtype=torch.bool).squeeze()
test_mask = torch.tensor(data['test_mask'], dtype=torch.bool).squeeze()

tensor_dataset = TensorDataset(X, y, train_mask, val_mask, test_mask)
torch.save(tensor_dataset, edge_output_dataset)

