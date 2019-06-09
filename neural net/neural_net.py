import torch
import torch.nn as nn
import pandas as pd

dataset = pd.read_csv('./salaries.csv')

x_temp = dataset.iloc[:, :-1].values
y_temp = dataset.iloc[:, 1:].values

X_train = torch.FloatTensor(x_temp)
Y_train = torch.FloatTensor(y_temp)
class Model(nn.Module):

    def __init__(self):
        super().__init__()
        self.features = nn.Sequential(
                nn.Linear(1, 128),
                nn.ReLU(),
                nn.Linear(128, 10),
                nn.ReLU(),
                nn.Linear(10, 1)
                )
    def forward(self, x):
        return self.features(x)


num_epochs = 10000
learning_rate = 0.01

model = Model()

criterion = nn.MSELoss()
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)

### TRAINING
for epoch in range(num_epochs):

    y_pred = model(X_train)

    loss = criterion(y_pred, Y_train)

    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    print(loss.item())


