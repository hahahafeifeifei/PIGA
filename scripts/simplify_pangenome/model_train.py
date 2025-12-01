import sys
import torch
from torch.utils.data import TensorDataset, DataLoader
from sklearn.metrics import roc_curve,auc
import torch.utils.data as data
from torch import nn
import numpy as np

dataset = torch.load(sys.argv[1])
torch.manual_seed(42)

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

X_train_tensor = dataset.tensors[0][dataset.tensors[2]]
X_val_tensor = dataset.tensors[0][dataset.tensors[3]]

y_train_tensor = dataset.tensors[1][dataset.tensors[2]]
y_val_tensor = dataset.tensors[1][dataset.tensors[3]]
del dataset

train_dataset = TensorDataset(X_train_tensor, y_train_tensor)
train_loader = DataLoader(train_dataset, batch_size=524288, shuffle=True, num_workers=8)
val_dataset = TensorDataset(X_val_tensor, y_val_tensor)
val_loader = DataLoader(val_dataset, batch_size=524288, shuffle=False, num_workers=8)

model = MLP(input_size=X_train_tensor.shape[1], hidden_size=256)
criterion = nn.BCEWithLogitsLoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)

# Training loop
num_epochs = 100
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

train_losses = []
val_losses = []
val_accuracies = []


for epoch in range(num_epochs):
    model.train()
    train_loss = 0.0
    batch = 0
    for inputs, labels in train_loader:
        inputs = inputs.to(device)
        labels = labels.to(device)

        optimizer.zero_grad()
        outputs = model(inputs)
        loss = criterion(outputs.squeeze(), labels)
        loss.backward()
        optimizer.step()
        train_loss += loss.item()
    train_losses.append(train_loss / len(train_loader))

    # Evaluation on test set
    model.eval()
    with torch.no_grad():
        val_loss = 0
        correct = 0
        tp = 0
        fp = 0
        fn = 0
        tn = 0
        total = 0
        for inputs, labels in val_loader:
            inputs = inputs.to(device)
            labels = labels.to(device)
            outputs = model(inputs)
            loss = criterion(outputs.squeeze(), labels)
            val_loss += loss.item()

            predicted = torch.round(torch.sigmoid(outputs)).squeeze()
            total += labels.size(0)
            tp += ((predicted == labels) & (labels == 0)).sum().item()
            fp += ((predicted != labels) & (labels == 1)).sum().item()
            fn += ((predicted != labels) & (labels == 0)).sum().item()
            tn += ((predicted == labels) & (labels == 1)).sum().item()
            correct += (predicted == labels).sum().item()

        val_losses.append(val_loss / len(val_loader))
        val_ter = tp/(tp + fn) * 100 
        val_tvr = tn/(fp + tn) * 100
        val_accuracy = correct / total * 100
        val_accuracies.append(val_accuracy)

        print(f"Epoch [{epoch+1}/{num_epochs}], Train loss: {train_losses[-1]:.4f}, Val loss: {val_losses[-1]:.4f}, Val accuracy: {val_accuracies[-1]:.2f}%, Val TER: {val_ter:.2f}%, Val TVR: {val_tvr:.2f}%")

torch.save(model, sys.argv[2])

#AUROC threshold
model.eval()
loader = DataLoader(val_dataset, batch_size=524288, shuffle=False, num_workers=8)
outputs = np.empty(0)
for inputs in loader:
    outputs = np.concatenate((outputs, torch.sigmoid(model(inputs[0])).detach().squeeze().numpy()))
labels = y_val_tensor.numpy()
fpr, tpr, thresholds = roc_curve(labels, outputs)
auroc = auc(fpr,tpr)
print(f"Val AUROC: {auroc:.4f}")

fo = open(sys.argv[3], "w")
for i in range(len(tpr)):
    if tpr[i] > 0.95:
        fo.write(str(thresholds[i]) + '\n')
        break
fo.close()
