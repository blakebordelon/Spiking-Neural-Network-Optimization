% analysis of weights

refweights = csvread('refweights.txt');
optweights = csvread('optweights.txt');
naiveweights = csvread('naiveweights.txt');
controlweights = csvread('controlweights.txt');

refweights = refweights(:,1:size(refweights,1));
optweights = optweights(:,1:size(optweights,1));
naiveweights = naiveweights(:,1:size(naiveweights,1));
controlweights = controlweights(:,1:size(controlweights,1));


%size(refweights)


optdiff = optweights - refweights;
naivediff = naiveweights - refweights;
controldiff = controlweights - refweights;

% 
% for i=1:1:size(optdiff,1)
%     for j=1:1:size(optdiff,2)
%         if abs(optdiff(i,j)) > 1.5e-3
%             optdiff(i,j) = 1.5e-3;
%         end
%         if abs(naivediff(i,j)) > 1.5e-3
%             naivediff(i,j) = 1.5e-3;
%         end
%         if abs(controldiff(i,j)) > 1.5e-3
%             controldiff(i,j) = 1.5e-3;
%         end
% 
%         
%         
%     end
% end



%optdiff = optdiff.^2;
%naivediff = naivediff.^2;
%controldiff = controldiff.^2;


sumoptdiff = sum(sum(abs(optdiff)));
sumnaivediff = sum(sum(abs(naivediff)));
sumcontroldiff = sum(sum(abs(controldiff)));

diffs = [sumnaivediff, sumoptdiff, sumcontroldiff]

naiveFlat = [];
optFlat = [];
controlFlat = [];

for i=1:1:size(optdiff,1)
    for j=1:1:size(optdiff,2)
        naiveFlat = [naiveFlat, naivediff(i,j)];
        optFlat = [optFlat, optdiff(i,j)];
        controlFlat = [controlFlat, controldiff(i,j)];
        
    end
end

meanNaiveFlatBeforeAbs = mean(naiveFlat)
meanOptFlatBeforeAbs = mean(optFlat)
meanControlFlatBeforeAbs = mean(controlFlat)

meanNaiveFlatAfterAbs = mean(abs(naiveFlat))
meanOptFlatAfterAbs = mean(abs(optFlat))
meanControlFlatAfterAbs = mean(abs(controlFlat))


stdN = std(abs(naiveFlat))
stdO = std(abs(optFlat))
stdC = std(abs(controlFlat))




csvwrite('weightDiffs.csv', diffs);


% %imagesc(heat);
% 
% 
set(gca,'Ydir','Normal')

subplot(2,2,1)
%imagesc(refweights);
%title('Reference Weights');
imagesc(naivediff.^2);
title('$\epsilon_{ij} = (N_{ij} - R_{ij})^2$', 'Interpreter', 'LaTeX', 'Fontsize' ,15);
colorbar

subplot(2,2,2)
%imagesc(naiveweights);
%title('Naive weights');
imagesc(optdiff.^2);
title('$\epsilon_{ij} = (O_{ij} - R_{ij})^2$', 'Interpreter', 'LaTeX', 'Fontsize' ,15);
colorbar

subplot(2,2,3)
%imagesc(optweights);
%title('Optimized Weights');
imagesc(controldiff.^2);
title('$\epsilon_{ij} = (C_{ij} - R_{ij})^2$', 'Interpreter', 'LaTeX', 'Fontsize' ,15);
colorbar
%subplot(2,2,4)
%imagesc(controlweights);
%title('Control weights');
% 
% 
% 
% 
% % refmean = mean(mean(refweights))
% % optmean = mean(mean(optweights))
% % naivemean = mean(mean(naiveweights))
% % controlmean = mean(mean(controlweights))
% 
% % opt
% 
% % refgraph = digraph(refweights);
% 
% %D1 = indegree(refgraph)
% %D2 = indegree(refgraph)
% %D3 = outdegree(refgraph)
% 
% %count = 0
% %for i = 1:1:size(refweights,1)
%    % if D1(i) ~= 400
%   %      count = count + 1;
%  %   end
% %end

