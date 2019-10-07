function varout = convertF(numNeurons)

F = csvread('fluorescence2.csv');
F=F';


new = [];
rows = numNeurons;
buckets = size(F,2)/rows;
alltraces = [];
for i=1:1:rows
    temp = [];
    for j=1:1:buckets
        temp = [temp F(1, ((i-1)*buckets + j))];
    end
    temp = temp';
    new = [new,temp];
end



passV.dt = .025;
passV.T = buckets;
passV.fast_thr = 1;
passV.est_a = 1;
passV.est_b = 1;
passV.est_gam = 1;
passV.est_sig = 1;
passV.est_lam = 1;
passV.Ncells = 1;


alltraces = [];
allparams = [];


desiredMeanRate = 1;
divisor = passV.dt * size(new,2)* size(new,1);
epsilon = .01;
nt = .1;
maxError = 1;
actualFiringRate = 100;




for i = 1:1:rows
    temp = fast_oopsi(new(:,i), passV); 
    alltraces = [alltraces, temp];

    %csvwrite('oopsi.csv', alltraces);
end

csvwrite('oopsi.csv', alltraces);
store = heaviside(alltraces - 3*repmat(std(alltraces')',1,size(alltraces,2)))



%{
iteration = 0;
while abs(actualFiringRate - desiredMeanRate) > maxError
    iteration = iteration + 1
    sum = 0;
    store = [];
    for i=1:1:size(alltraces,1)
        smallstore = [];
        for j=1:1:size(alltraces,2)
            if(alltraces(i, j) < nt)
                smallstore = [smallstore 0];
            
            else
                smallstore = [smallstore 1];
                sum = sum + 1;
            end
        end
        smallstore = smallstore';
        store = [store, smallstore];
    end
    actualFiringRate = sum / divisor
    if(actualFiringRate > desiredMeanRate)
        nt = nt + epsilon;
    else
        nt = nt - epsilon;
    end
end
        

%}

varout = store;