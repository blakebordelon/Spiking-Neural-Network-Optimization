clear
numNeurons = 400;
timeSteps = 5001;
%input name of pre-optimization spike train file
spikes1 = csvread('oldOptSpikes.txt');
size(spikes1)
%input name of post-optimization spike train file
spikes2 = csvread('newOptSpikes.txt');
%input name of reference spike train
spikes3pre = csvread('refSpikes.txt');
spikesControl = csvread('controlSpikes.txt');


%Comment out the three lines beginning with "weights" if you aren't
%analyzing a network structure

%input name of pre-optimization network file
%weightsStart = csvread('');
%input name of post-optimization network file
%weightsOut = csvread('weightsOutOptToSmallWorld.txt');
%input name of reference network file
%weightsSmallWorld = csvread('weightsSmallWorld.csv');



meanRef = mean(sum(spikes3pre,2))
meanNaive = mean(sum(spikes1,2))
meanOpt = mean(sum(spikes2,2))
meanControl = mean(sum(spikesControl,2))


meanvals = [meanRef,meanNaive, meanOpt,meanControl];

csvwrite('meanSpikes.csv', meanvals);


refVar = var(sum(spikes3pre,2));
naiveVar = var(sum(spikes1,2));
optVar = var(sum(spikes2,2));
controlVar = var(sum(spikesControl,2));

varVals = [refVar, naiveVar, optVar, controlVar];

csvwrite('varSpikes.csv', varVals);



refRates = [];
preOptRates = [];
controlRates = [];
postOptRates = [];

for i=1:1:numNeurons
    refRates = [refRates, sum(spikes3pre(i, :))];
    preOptRates = [preOptRates, sum(spikes1(i, :))];
    postOptRates = [postOptRates, sum(spikes2(i,:))];
    controlRates = [controlRates, sum(spikesControl(i, :))];
end


spikes3 = spikes3pre(1:numNeurons,1:timeSteps);

intervalsOld = [];

oldCVs = [];

for i=1:1:size(spikes1,1)
   isi = 0;
   intervalsrow = [];
   for j =1:1:size(spikes1,2)

        if(spikes1(i,j)==0)
            isi=isi+1;
        else
        %    if(isi<=100)
                intervalsrow = [intervalsrow,isi];
         %   end
            isi = 0;
        end
   end
    if mean(intervalsrow) ~= 0
        currCV = std(intervalsrow) / mean(intervalsrow);
        oldCVs = [oldCVs currCV];
    end
    
    intervalsOld = [intervalsOld,intervalsrow];
end

CVOLD = mean(oldCVs);

intervalsNew = [];

newCVs = [];

for i=1:1:size(spikes2,1)
   isi = 0;
   intervalsrow = [];
   for j =1:1:size(spikes2,2)

        if(spikes2(i,j)==0)
            isi=isi+1;
        else
         %   if(isi<=100)
                intervalsrow = [intervalsrow,isi];
          %  end
            isi = 0;
        end
   end
    if mean(intervalsrow) ~= 0
        cvcurrent = std(intervalsrow) / mean(intervalsrow);
        newCVs = [newCVs cvcurrent];
    end
    
    intervalsNew = [intervalsNew,intervalsrow];
end



CVNEW = mean(newCVs);

%intervalsNew = intervalsNew - 2*ones(size(intervalsNew,1),size(intervalsNew,2));

intervalsRef = [];

refCVs = [];

for i=1:1:size(spikes3,1)
   isi = 0;
   intervalsrow = [];
   for j =1:1:size(spikes3,2)

        if(spikes3(i,j)==0)
            isi=isi+1;
        else
          %  if(isi<=100)
                intervalsrow = [intervalsrow,isi];
           % end
            isi = 0;
        end
   end
    if mean(intervalsrow) ~= 0
        cvcurrent = std(intervalsrow) / mean(intervalsrow);
        refCVs = [refCVs cvcurrent];
    end
    
    intervalsRef = [intervalsRef,intervalsrow];
end


CVREF = mean(refCVs);

intervalsControl = [];

controlCVs = [];

for i=1:1:size(spikesControl,1)
   isi = 0;
   intervalsrow = [];
   for j =1:1:size(spikesControl,2)

        if(spikesControl(i,j)==0)
            isi=isi+1;
        else
          %  if(isi<=100)
                intervalsrow = [intervalsrow,isi];
           % end
            isi = 0;
        end
   end
    if mean(intervalsrow) ~= 0
        cvcurr = std(intervalsrow) / mean(intervalsrow);
        controlCVs = [controlCVs cvcurr];
    end
    intervalsControl = [intervalsControl,intervalsrow];
end

CVCONTROL = mean(controlCVs);

cvs = [CVREF CVOLD CVNEW CVCONTROL];

csvwrite('CVs.csv', cvs);

[oldtest pOld] = kstest2(refCVs, oldCVs);
[newtest pNew] = kstest2(refCVs, newCVs);
[controltest pControl] = kstest2(refCVs, newCVs);

pOld
pNew
pControl

ps = [pOld pNew pControl];
csvwrite('cvdistkstest.csv', ps);


binSize = 20;
Var1 = zeros(size(spikes1,1),round((size(spikes1,2)-1)/binSize));
Var2 = zeros(size(spikes1,1),round((size(spikes1,2)-1)/binSize));
Var3 = zeros(size(spikes1,1),round((size(spikes1,2)-1)/binSize));
Var4 = zeros(size(spikes1,1),round((size(spikes1,2)-1)/binSize));

Mean1 = zeros(size(spikes1,1),round((size(spikes1,2)-1)/binSize)); 
Mean2 = zeros(size(spikes1,1),round((size(spikes1,2)-1)/binSize));
Mean3 = zeros(size(spikes1,1),round((size(spikes1,2)-1)/binSize));
Mean4 = zeros(size(spikes1,1),round((size(spikes1,2)-1)/binSize));

for i=1:1:size(spikes1,1)
    for j=1:1:(size(spikes1,2)-1)/binSize
        Var1(i,j) = var(spikes1(i,(((j-1)*binSize)+1):(j*binSize)));
        Mean1(i,j) = mean(spikes1(i,((j-1)*binSize)+1:j*binSize));
    end
end

for i=1:1:size(spikes2,1)
    for j=1:1:(size(spikes2,2)-1)/binSize
        Var2(i,j) = var(spikes2(i,((j-1)*binSize)+1:j*binSize));
        Mean2(i,j) = mean(spikes2(i,((j-1)*binSize)+1:j*binSize));
    end
end

for i=1:1:size(spikes3,1)
    for j=1:1:(size(spikes3,2)-1)/binSize
        Var3(i,j) = var(spikes3(i,((j-1)*binSize)+1:j*binSize));
        Mean3(i,j) = mean(spikes3(i,((j-1)*binSize)+1:j*binSize));
    end 
end

for i=1:1:size(spikesControl,1)
    for j=1:1:(size(spikesControl,2)-1)/binSize
        Var4(i,j) = var(spikesControl(i,((j-1)*binSize)+1:j*binSize));
        Mean4(i,j) = mean(spikesControl(i,((j-1)*binSize)+1:j*binSize));
    end 
end


%{
CV1 = zeros(size(Var1));
CV2 = zeros(size(Var2));
CV3 = zeros(size(Var3));
CV4 = zeros(size(Var4));



for i=1:1:size(Var1)
    
    if Mean1 ~= 0
        CV1(i) = Var1(i)/Mean1;
    else
        if Var1(i) == 0
            CV1(i) = 1;
        end
       
    end
       
    if Mean2 ~= 0
        CV2(i) = Var2(i)/Mean2;
    else
        if Var2(i) == 0
            CV2(i) = 1;
        end
       
       end
    if Mean3 ~= 0
        CV3(i) = Var3(i)/Mean3;
    else
        if Var3(i) == 0
            CV3(i) = 1;
        end
       
    end
    if Mean4 ~= 0
        CV4(i) = Var4(i)/Mean4;
    else
        if Var4(i) == 0
            CV4(i) = 1;
        end
       
    end
end

CV1Size = size(CV1)

%}


CV1 = Var1./Mean1;
CV2 = Var2./Mean2;
CV3 = Var3./Mean3;
CV4 = Var4./Mean4;



Fano1 = CV1.^-1;
Fano2 = CV2.^-1;
Fano3 = CV3.^-1;
Fano4 = CV4.^-1;

CV1mean = nanmean(CV1,1);
CV2mean = nanmean(CV2,1);
CV3mean = nanmean(CV3,1);
CV4mean = nanmean(CV4,1);

spikeCountPre = sum(sum(spikes1,1))
spikeCountPost = sum(sum(spikes2,1))
spikeCountRef = sum(sum(spikes3,1))

spikeCountControl = sum(sum(spikesControl,1))

%mean weights for both networks
%smallMean = mean(mean(weightsSmallWorld))
%optMean = mean(mean(weightsOut))
%startMean = mean(mean(weightsStart));

postError = 0;
preError = 0;

%outNew = histogram(intervalsNew, 100);
%outOld = histogram(intervalsOld, 100);
%outRef = histogram(intervalsRef, 100);


%outNew.Values




subplot(2,2,3);
hist(intervalsNew,100)
title('Optimized')
ylabel('Frequency')
xlabel('InterSpike Interval (time steps)')
ylim([0 2000])
%xlim([0 200])


subplot(2,2,2); 
hist(intervalsOld,100)
title('Naive')
ylabel('Frequency')
xlabel('InterSpike Interval (time steps)')
ylim([0 2000])
%xlim([0 200])

subplot(2,2,1)
histogram(intervalsRef,100)
title('Reference')
ylabel('Frequency')
xlabel('InterSpike Interval (time steps)')
ylim([0 2000])
%xlim([0 200])

subplot(2,2,4)
histogram(intervalsControl,100)
title('Control')
ylabel('Frequency')
xlabel('InterSpike Interval (time steps)')
ylim([0 2000])
%xlim([0 200])



%figure(1)
%hist(intervalsNew,100)
%title('D. Post-Optimization')
%ylabel('Frequency')
%xlabel('InterSpike Interval (time steps)')
%ylim([0 2000])
%xlim([0 100])


figure(2)
hist(intervalsOld,75)
title('ISI Distribution Pre-Optimization')
ylabel('Frequency')
xlabel('InterSpike Interval (time steps)')
ylim([0 2000])
xlim([0 100])


figure(3)
histogram(intervalsRef,100)
title('ISI Distribution for Reference Network')
ylabel('Frequency')
xlabel('InterSpike Interval (time steps)')
ylim([0 2000])
xlim([0 100])


figure(4)
histogram(intervalsControl,100,'BinLimits',[1,100])
title('ISI Distribution for Control Network')
ylabel('Frequency')
xlabel('InterSpike Interval (time steps)')
ylim([0 2000])
xlim([0 100])

binLim = 200;



figure(5)
newHist = histogram(intervalsNew,200,'BinLimits',[1,binLim],'facealpha',0.5,'facecolor','b')
hold on
oldHist = histogram(intervalsOld,200,'BinLimits',[1,binLim],'facealpha',0.5,'facecolor','r')
refHist = histogram(intervalsRef,200,'BinLimits',[1,binLim],'facealpha',0.5,'facecolor','y')
controlHist = histogram(intervalsControl,200,'BinLimits',[1,binLim],'facealpha',0.5,'facecolor','y')

%histogram(intervalsControl, 100, 'facealpha',0.5, 'facecolor', 'p');
%legend('Post-Optimization ISI distribution','Pre-Optimization ISI distribution','Reference ISI distribution')
ylabel('Frequency')
xlabel('InterSpike Interval')
ylim([0,5000])

%controlHist = histogram(intervalsControl, edges);

oldVals = oldHist.Values;
newVals = newHist.Values;
refVals = refHist.Values;
controlVals = controlHist.Values;

oldError = 0;
newError = 0;

controlError = 0;




for i=1:1:binLim
    
    oldError = oldError + (oldVals(i)-refVals(i))^2 /(refVals(i));
    newError = newError + (newVals(i) - refVals(i))^2 /(refVals(i));
    controlError = controlError + (controlVals(i) - refVals(i))^2 / (refVals(i));
    
end

CV1Error = 0;
CV2Error = 0;
CV4Error = 0;
for i=1:1:size(CV1mean)
    
    CV1Error = CV1Error + (CV1mean(i) - CV3mean(i))^2 / (CV3mean(i));
    CV2Error = CV2Error + (CV2mean(i) - CV3mean(i))^2 / (CV3mean(i));
    CV4Error = CV4Error + (CV4mean(i) - CV3mean(i))^2 / (CV3mean(i));
    
    
end

KLISIOld = 0;
KLISINew = 0;
KLISICont = 0;

sumOld = sum(oldVals);
sumNew = sum(newVals);
sumRef = sum(refVals);
sumControl = sum(controlVals);

%{
for i=1:1:binLim
    
    if oldVals(i)  ~=0 && refVals(i) ~= 0
        KLISIOld = KLISIOld + (oldVals(i)/sumOld)*( log2(oldVals(i)) - log2(refVals(i)) - log2(sumOld)+log2(sumRef));
    else
        if refVals(i) == 0 && oldVals(i) ~=0
        end
    end

    if newVals(i)  ~=0 && refVals(i) ~= 0
        KLISINew = KLISINew + (newVals(i)/sumNew)*(log2(newVals(i)) - log2(refVals(i)) - log2(sumNew)+log2(sumRef));
    
    end    
    
    if controlVals(i)  ~=0 && refVals(i) ~= 0
        KLISICont = KLISICont + (controlVals(i)/sumControl)*(log2(controlVals(i)) - log2(refVals(i)) - log2(sumControl)+log2(sumRef));
    
    end    
    
end

KLISIOld
KLISINew
KLISICont

CV1Error = 0;
CV2Error = 0;
CV4Error = 0;
for i=1:1:size(CV1mean)
    
    CV1Error = CV1Error + (CV1mean(i) - CV3mean(i))^2 / (CV1mean(i) + CV3mean(i));
    CV2Error = CV2Error + (CV2mean(i) - CV3mean(i))^2 / (CV2mean(i) + CV3mean(i));
    CV4Error = CV4Error + (CV4mean(i) - CV3mean(i))^2 / (CV4mean(i) + CV3mean(i));
    
    
end


%}

figure(6)
plot(CV1mean)
hold on
plot(CV2mean)
plot(CV3mean)
plot(CV4mean)
hold off
xlabel('Time bin (each bin corresponds to 20 ms)')
ylabel('CV')
legend('Before optimization','After optimization','Reference','Control')

figure(7)
histogram(sum(spikes1,2),50,'facealpha',0.5,'facecolor','b')
hold on
histogram(sum(spikes2,2),50,'facealpha',0.5,'facecolor','r')
histogram(sum(spikes3,2),50,'facealpha',0.5,'facecolor','y')
legend('Before Optimization',50,'After Optimization','Reference')
ylim([0,200])
hold off

%{
figure(7)
histogram(intervalsControl,100)
title('ISI Distribution for Control Network')
ylabel('Frequency')
xlabel('InterSpike Interval')
%}

timeScale = 1/(timeSteps * .025)

%refRates = refRates * timeScale;
%preOptRates = preOptRates * timeScale;
%postOptRates = postOptRates * timeScale;
%controlRates = controlRates * timeScale;


figure(8)
subplot(2,2,1)
histogram(refRates,100 ,'BinLimits',[1,100],'facealpha',0.5,'facecolor','b')
xlabel('spike number')
title('Reference')
hold on
subplot(2,2,2)
histogram(preOptRates,100,'BinLimits',[1,100],'facealpha',0.5,'facecolor','r')
xlabel('spike number')
title('Naive')
subplot(2,2,3)
histogram(postOptRates,100,'BinLimits',[1,100],'facealpha',0.5,'facecolor','y')
xlabel('spike number')
title('Optimized')
subplot(2,2,4)
histogram(controlRates,100,'BinLimits',[1,100],'facealpha',0.5,'facecolor','g')
xlabel('spike number');
title('Control')
%legend('reference','naive','optimized','control')
hold off

% normalize the error with respect to the total number of spikes in the
% network

numRefSpikes = sum(sum(spikes3pre));

oldError = oldError / numRefSpikes
newError = newError / numRefSpikes
controlError = controlError / numRefSpikes


CV1Error
CV2Error
CV4Error


%[oldKS pOld] = kstest2(oldVals, refVals);
%[newKS pNew] = kstest2(newVals, refVals);
%[conKS pCon] = kstest2(controlVals, refVals);

[oldKS pOld] = kstest2(intervalsOld, intervalsRef);
[newKS pNew] = kstest2(intervalsNew, intervalsRef);
[conKS pCon] = kstest2(intervalsControl, intervalsRef);



[oldCVKS pCVOld] = kstest2(CV1mean, CV3mean);
[newCVKS pCVNew] = kstest2(CV2mean, CV3mean);
[conCVKS pCVCon] = kstest2(CV4mean, CV3mean);

log(pOld)
log(pNew)
log(pCon)


pCVOld
pCVNew
pCVCon


ksISIpVals = [pOld pNew pCon];
ksCVpVals = [pCVOld pCVNew pCVCon];

csvwrite('KSTestISI.csv', ksISIpVals);
csvwrite('KSTestCV.csv', ksCVpVals);

bothErrors = [oldError newError controlError];

CVErrors = [CV1Error CV2Error CV4Error];
%bothErrors = [oldError newError controlError];

realCVOld = std(oldVals)/mean(oldVals);
realCVNew = std(newVals)/mean(newVals);
realCVRef = std(refVals)/mean(refVals);
realCVControl = std(controlVals)/mean(controlVals);

allCVReal = [realCVOld, realCVNew, realCVRef, realCVControl];

csvwrite('isiDistErrorsAttempt2.csv', bothErrors);
csvwrite('cvDistErrorsAttempt2.csv',CVErrors);
csvwrite('realCV.csv', allCVReal);


figure(9)
hist(refCVs, 100);
hold on
hist(oldCVs, 100);
hist(newCVs, 100);
hist(controlCVs, 100);

%intervalsRef
