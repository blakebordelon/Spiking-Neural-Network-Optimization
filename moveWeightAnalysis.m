
for c=1:1:30

    
    if c ~= 21 
        count = c;
        %string1 = ['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparisty = .5\trial ' num2str(count)];
        string1 = ['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\alpha = .5\trial ' num2str(count)];

        cd(string1)
    
        weightanalysis;
    
    end
end

c = 1;
while c < 31
    count = c
    if c ~=21 
    
    
    
    string1 = ['C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\alpha = .5\trial ' num2str(count)]
    
    count = c
    
    cd(string1) 
    
    %vanrossum;
    %isiCompSmallWorld;
    
%    movefile('KSTestISI.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
%    movefile('localVanRossum.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
%    movefile('meanSpikes.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
%    movefile('totalVanRossum.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')
%    movefile('varSpikes.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\sparse\sparsity = .5')

    movefile('weightdiffs.csv','C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\alpha = .5')

    
    
    cd('C:\Users\Blake\OneDrive (2)\Documents\Wash U\Research\Final revision results\differentWeightDists\alpha = .5')

    count = c
    isiString = ['weightdiffs' num2str(count) '.csv'];
    
    movefile('weightdiffs.csv', isiString)
    
%     count = c
%     
%     locString = ['localVanRossum' num2str(count) '.csv'];
%     movefile('localVanRossum.csv', locString)
%     count = c
%     movefile('meanSpikes.csv', ['meanSpikes' num2str(count) '.csv'])
%     count = c
%     movefile('totalVanRossum.csv', ['totalVanRossum' num2str(count) '.csv'])
%     count = c
%     movefile('varSpikes.csv', ['varSpikes' num2str(count) '.csv'])

    end
    c = c+1
    
    
end