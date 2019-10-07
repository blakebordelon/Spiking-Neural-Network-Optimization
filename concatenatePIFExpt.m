
allReference = [];
allNaive = [];
allOpt = [];
allControl = [];


for i=1:1:5
    currentstring = ['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\pif\trial ' num2str(i)];
    cd(currentstring);
    ref = csvread('refSpikes.txt');
    opt = csvread('newOptSpikes.txt');
    naive = csvread('oldOptSpikes.txt');
    control = csvread('controlSpikes.txt');
    
    allReference = [allReference, ref];
    allNaive = [allNaive, naive];
    allOpt = [allOpt, opt];
    allControl = [allControl, control];
    
end

