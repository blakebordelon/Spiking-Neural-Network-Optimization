
spikes = csvread('spikes.csv');

figure(1)
imagesc(spikes(:,1:5000))
title('PIF Spike Trains')
xlabel('Time-Step')
ylabel('Neuron ID')