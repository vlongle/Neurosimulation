import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from torch.utils.data import Dataset, DataLoader, TensorDataset
import torch


class SpikeDataset(Dataset):
  def __init__(self, file_path):
      # first 6 columns are inputs
    self.data_x = np.loadtxt(file_path,
      usecols=range(0,6), dtype=np.float32)

    # last 4 columns are labels
    self.data_y = np.loadtxt(file_path,
            usecols=range(6,10), dtype=np.float64)

  def __getitem__(self, index):
    features = torch.FloatTensor(self.data_x[index])
    labels = torch.FloatTensor(self.data_y[index])
    return (features, labels)

  def __len__(self):
    return len(self.data_x)


train_file = 'training_data.txt'
#data_x = np.loadtxt(train_file, usecols=range(0,6))
#data_y = np.loadtxt(train_file, usecols=range(6,10))
#
#print(type(data_x))
#train_dataset = TensorDataset(data_x, data_y)

def get_data():
    train_dataset = SpikeDataset(train_file)
    train_loader = DataLoader(train_dataset, batch_size=16,
        shuffle=True)

    return train_loader

#for i_batch, data in enumerate(train_loader,0):
#    X, Y = data
#    print(X.size())
#    print(Y.size())
#

#class IrisDataset(Dataset):
#    def __init__(self, dataset,)
#dataset = pd.read_csv('iris.csv')
#
#dataset.loc[dataset.species=='Iris-setosa', 'species'] = 0
#dataset.loc[dataset.species=='Iris-versicolor', 'species'] = 1
#dataset.loc[dataset.species=='Iris-virginica', 'species'] = 2
#
#
#
#num_epochs = 1000
#all_losses = []
#
#for epoch in range(num_epochs):
#    for :
#        pass
#
#        all_losses.append(loss.data)
#
#
#
#
#all_losses = np.array(all_losses, dtype = np.float)
#plt.plot(all_losses)
#plt.show()
#
#
#
