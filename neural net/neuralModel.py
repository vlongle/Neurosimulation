import torch.nn as nn
import torch.nn.functional as F
import torch
from dataGetter import get_data
import numpy as np
from sklearn.metrics import mean_squared_error
import sys


class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()

        self.fc1 = nn.Linear(6,30)

        self.fc2 = nn.Linear(30, 4)


    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x


num_epoch = int(sys.argv[1])

net = Net()

criterion = nn.MSELoss()

optimizer = torch.optim.SGD(net.parameters(), lr=0.001)

train_loader = get_data()

train_file = 'training_data.txt'
### CHECK SOME DATA
test_X = np.loadtxt(train_file, usecols=range(0,6), max_rows=4)
test_Y = np.loadtxt(train_file, usecols=range(6, 10), max_rows=4)

net.eval()
predict_Y = net(torch.FloatTensor(test_X)).detach()
print('test_Y \n', test_Y)
print('predict \n', predict_Y)

mse = mean_squared_error(test_Y, predict_Y)
print('MSE loss:', mse)


net.train()
for epoch in range(num_epoch):
    for i_batch, data in enumerate(train_loader, 0):
        X, Y = data

        optimizer.zero_grad()
        out = net(X)

        loss = criterion(out, Y)
        loss.backward()
        optimizer.step()

    if epoch % 50 == 0:
        print ('Epoch [%d/%d] Loss: %.4f'
                   %(epoch+1, num_epoch, loss.item()))

### POST-TRAINING
net.eval()
predict_Y = net(torch.FloatTensor(test_X)).detach()
print('test_Y \n', test_Y)
print('predict \n', predict_Y)



mse = mean_squared_error(test_Y, predict_Y)
print('MSE loss:', mse)
