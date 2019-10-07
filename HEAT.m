% heat map 
clear;

refSpikes = csvread('refSpikes.txt');
controlSpikes = csvread('controlSpikes.txt');
newOptSpikes = csvread('newOptSpikes.txt');
oldOptSpikes = csvread('oldOptSpikes.txt');
%csvSpikes = csvread('spikes.csv');

refSpikes = -1 * refSpikes(:,1:5000);
controlSpikes = -1 * controlSpikes(:,1:5000);
newOptSpikes = -1 * newOptSpikes(:,1:5000);
oldOptSpikes = -1 * oldOptSpikes(:,1:5000);

size(refSpikes)

heat = [1,2,2,1;5,3,1,0;2,2,2,2];
imagesc(heat);

myMap = [0,0,0;1,1,1];


set(gca,'Ydir','Normal')

subplot(2,2,1)
imagesc(refSpikes);
colormap(myMap)
title('\fontsize{20}Reference');
xlabel('\fontsize{20}Time Steps');
ylabel('\fontsize{20}Neuron ID');


subplot(2,2,2)
imagesc(oldOptSpikes);
colormap(myMap)
title('\fontsize{20}Naive');
xlabel('\fontsize{20}Time Steps');
ylabel('\fontsize{20}Neuron ID');


subplot(2,2,3)
imagesc(newOptSpikes);
colormap(myMap)
title('\fontsize{20}Optimized');
xlabel('\fontsize{20}Time Steps');
ylabel('\fontsize{20}Neuron ID');

subplot(2,2,4)
imagesc(controlSpikes);
title('\fontsize{20}Control');
xlabel('\fontsize{20}Time Steps');
ylabel('\fontsize{20}Neuron ID');



% subplot(2,2,1)
% %h1 = HeatMap(refSpikes, 'ColorMap', redbluecmap);
% 
% %h1 = HeatMap(refSpikes, 'ColorMap', redbluecmap);
% %addTitle(h1, 'A. Reference');
% %addXLabel(h1, 'Time');
% %addYLabel(h1, 'Neuron ID');
% 
% 
% 
% 
% subplot(2,2,3)
% h2 = HeatMap(newOptSpikes,'ColorMap', redbluecmap);
% addTitle(h2, 'C. Post-Optimization');
% addXLabel(h2, 'Time');
% addYLabel(h2, 'Neuron ID');
% 
% subplot(2,2,2)
% h3 = HeatMap(oldOptSpikes, 'ColorMap', redbluecmap);
% addTitle(h3, 'B. Pre-Optimization');
% addXLabel(h3, 'Time');
% addYLabel(h3, 'Neuron ID');
% 
% subplot(2,2,4)
% h4 = HeatMap(controlSpikes, 'ColorMap',redbluecmap);
% addTitle(h4, 'D. Control');
% addXLabel(h4, 'Time');
% addYLabel(h4, 'Neuron ID');
% %HeatMap(csvSpikes, 'ColorMap',redbluecmap);
