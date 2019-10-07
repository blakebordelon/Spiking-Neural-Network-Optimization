spikes = [];

for i=1:1:30
    
    string1 = ['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\PSPM Data Files\sparse\sparsity = .5\trial ' num2str(i)]
    
    cd(string1) 

    controlSpikes = csvread('controlSpikes.txt');
    refSpikes = csvread('refSpikes.txt');
    optSpikes = csvread('newOptSpikes.txt');
    naiveSpikes = csvread('oldOptSpikes.txt');
    
    
    string2 = ['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\PSPM Data Files Complete\sparse\trial ' num2str(i)]
    
    mkdir(string2);
    cd(string2);
    csvwrite('controlSpikes.csv', round(controlSpikes));
    csvwrite('refSpikes.csv', round(refSpikes));
    csvwrite('newOptSpikes.csv', optSpikes);
    csvwrite('naiveSpikes.csv', naiveSpikes);
    
end
