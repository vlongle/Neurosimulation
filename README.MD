#NEUROSIMULATION/DEEP LEARNING
Given a small population of biological neurons, we have a stochastic simulation to model how that population would react to an external current drive and interact among themselves. Some characteristics of the population like the synapse connection strength between neurons would largely influence the its dynamics. In this project, we aim to use deep learning and leverage the power of big data to infer those characteristics.

## Progress
1. _Summary Statistics_
First, we vary the amount of external current delivered to the population and record the summary statistics: the mean firing rate. Using those data, we infer the synapse strength by a deep feed-forward neural network with a mean-absolute-error of about 0.03. 
2. __Time Series__

