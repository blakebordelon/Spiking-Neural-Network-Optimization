
for c=1:1:29

    
    if c ~= 100
        count = c
        string1 = ['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5\trial ' num2str(count)]
        cd(string1)
    
        vanrossum;
        isiCompSmallWorld;
    
    end
    
    
end

c = 1;
while c < 30
    count = c
    if c ~=100
    
    
    
    string1 = ['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5\trial ' num2str(count)]
    
    count = c
    
    cd(string1) 
    
    %vanrossum;
    %isiCompSmallWorld;
    
    movefile('KSTestISI.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
    movefile('localVanRossum.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
    movefile('meanSpikes.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
    movefile('totalVanRossum.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
    movefile('varSpikes.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
    
    
    
    cd('C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')

    count = c
    isiString = ['KSTestISI' num2str(count) '.csv'];
    
    movefile('KSTestISI.csv', isiString)
    
    count = c
    
    locString = ['localVanRossum' num2str(count) '.csv'];
    movefile('localVanRossum.csv', locString)
    count = c
    movefile('meanSpikes.csv', ['meanSpikes' num2str(count) '.csv'])
    count = c
    movefile('totalVanRossum.csv', ['totalVanRossum' num2str(count) '.csv'])
    count = c
    movefile('varSpikes.csv', ['varSpikes' num2str(count) '.csv'])

    end
    c = c+1
    
    
end

% 
% cd(['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\gaussian\alpha, beta = 1.5\trial 2']) 
% 
% vanrossum;
% isiCompSmallWorld;
% 
% movefile('KSTestISI.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\gaussian\alpha, beta = 1.5')
% 
% cd(['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\gaussian\alpha, beta = 1.5']) 
% 
% movefile('KSTestISI.csv', 'KSTestISI2.csv');
