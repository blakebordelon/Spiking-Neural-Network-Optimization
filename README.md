# Supervised-Learning-Rules-in-Spiking-Neural-Networks 
spikesMatchOptInputFixed.c contains an implementation of the spike-matching algorithm described here: https://arxiv.org/abs/1810.03199
It generates a reference and naive weight matrix for Leaky Integrate-and-Fire (LIF) neural networks and changes the weights of the naive matrix to eliminate spikes unpaired after application of a dynamic string matching algorithm. 

spikesMatchOptInputSpikeData.c allows the user to pass in spike data from an external source. In the manuscript referenced above, this file was used to discover a LIF network that reproduced critical PIF neural activity.

vanRossum.m is a Matlab script that calculates the adapted van-Rossum distance measure used in the text. 
