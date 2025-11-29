import sys
import torch
from torch import nn
import numpy as np
import pandas as pd
from torch.utils.data import TensorDataset, DataLoader
import torch.utils.data as data
from sklearn.metrics import roc_curve


def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

class ResidualBlock(nn.Module):
    def __init__(self, hidden_size):
        super(ResidualBlock, self).__init__()
        self.fc1 = nn.Linear(hidden_size, hidden_size)
        self.bn1 = nn.LayerNorm(hidden_size)
        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.bn2 = nn.LayerNorm(hidden_size)
    def forward(self, x):
        identity = x
        x = torch.relu(self.bn1(self.fc1(x)))
        x = self.bn2(self.fc2(x))
        return torch.relu(x + identity)

class MLP(nn.Module):
    def __init__(self, input_size, hidden_size):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.bn1 = nn.LayerNorm(hidden_size)
        self.blocks1 = nn.Sequential(*[ResidualBlock(hidden_size) for i in range(4)])
        self.fc2 = nn.Linear(hidden_size, 1)
    def forward(self, x):
        x = self.bn1(self.fc1(x))
        x = self.blocks1(x)
        x = self.fc2(x)
        return x

####Load the data
input_file = sys.argv[1]
edge_name_list = ["ID","Trans_ID1","Trans_ID2","HET","HOM","AD","SN","Mini","Homopolymer"] + ["Graphlet" + str(i) for i in range(0,68)] + ["training_y"]
data = pd.read_csv(input_file, sep='\t', header=None, names=edge_name_list, na_values=['NONE','nan'])

####Fill the missing value
statistics_file = sys.argv[2]
feature_dict = {}
for line in openfile(statistics_file):
    feature = line.split()[0]
    feature_dict[feature] = [ float(statistic) for statistic in line.strip().split()[1:] ]
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
loader = DataLoader(TensorDataset(X), batch_size=524288, shuffle=False)

model_file = sys.argv[3]
model = torch.load(model_file, map_location=torch.device('cpu'))
model.eval()
threshold = float(sys.argv[4]) 

outputs = []
for inputs in loader:
    outputs += torch.sigmoid(model(inputs[0])).detach().squeeze().tolist()
labels = ["TP" if output >= threshold else "FP" for output in outputs]

for i in range(len(data)):
    ID = data["ID"][i]
    label = labels[i]
    print(ID + "\t" + label)



