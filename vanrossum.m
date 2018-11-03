%Van-Rossum Distance

ref = csvread('refSpikes.txt');
opt = csvread('newOptSpikes.txt');
naive = csvread('oldOptSpikes.txt');
control = csvread('controlSpikes.txt');

%refRobust = csvread('refRobust.txt');
%optRobust = csvread('optRobust.txt');
%naiveRobust = csvread('naiveRobust.txt');
%controlRobust = csvread('controlRobust.txt');


x = linspace(-20,20);
filter = exp(- .01 * x.^2);
opterror = 0;
naiveerror = 0;
controlerror = 0;
refSig = 0;
optSig = 0;

totalRef = 0;
totalOpt = 0;
totalNaive = 0;
totalControl = 0;


for i=1:1:size(ref,1)
    
    refTrain = ref(i,:);
    optTrain = opt(i,:);
    naiveTrain = naive(i,:);
    controlTrain = control(i,:);
    
    %refRobustTrain = refRobust(i,:);
    %optRobustTrain = optRobust(i,:);
    %naiveRobustTrain = naiveRobust(i,:);
    %controlRobustTrain = controlRobust(i,:);
    
    
    
    refSig  = conv(filter, refTrain);
    optSig = conv(filter, optTrain);
    naiveSig = conv(filter, naiveTrain);
    controlSig = conv(filter, controlTrain);
    
    %refRobustSig  = conv(filter, refTrain);
    %optRobustSig = conv(filter, optTrain);
    %naiveRobustSig = conv(filter, naiveTrain);
    %controlRobustSig = conv(filter, controlTrain);

    
    totalRef = totalRef + refSig;
    totalOpt = totalOpt + optSig;
    totalNaive = totalNaive + naiveSig;
    totalControl = totalControl + controlSig;
    
    opterror = opterror + sum((refSig - optSig).^2);
    naiveerror = naiveerror + sum((refSig - naiveSig).^2);
    controlerror = controlerror + sum((refSig - controlSig).^2);

end




opterror
naiveerror
controlerror

totalNaiveError = sum((totalRef - totalNaive).^2);
totalOptError = sum((totalRef - totalOpt).^2);
totalControlError = sum((totalRef - totalControl).^2);


localerrors = [naiveerror, opterror,controlerror]
totalerrors = [totalNaiveError, totalOptError,totalControlError]

csvwrite('localVanRossum.csv', localerrors);
csvwrite('totalVanRossum.csv', totalerrors);






figure(1)
subplot(2,2,1)
plot(totalRef)
title('\fontsize{20}Reference Network Activity')
xlabel('\fontsize{20}Time Step')
ylabel('$A(\mathcal{R},t)$','Interpreter','LaTeX','Fontsize',20);
ylim([0,800])
xlim([0,5000])

subplot(2,2,3)
plot(totalOpt)
title('\fontsize{20}Optimized Network Activity')
xlabel('\fontsize{20}Time Step')
ylabel('$A(\mathcal{O},t)$','Interpreter','LaTeX','Fontsize',20);
ylim([0,800])
xlim([0,5000])


subplot(2,2,2)
plot(totalNaive)
title('\fontsize{20}Naive Network Activity')
xlabel('\fontsize{20}Time Step')
ylabel('$A(\mathcal{N},t)$','Interpreter','LaTeX','Fontsize',20);
ylim([0,800])
xlim([0,5000])


subplot(2,2,4)
plot(totalControl)
title('\fontsize{20}Control Network Activity')
xlabel('\fontsize{20}Time Step')
ylabel('$A(\mathcal{C},t)$','Interpreter','LaTeX','Fontsize',20);
ylim([0,800])
xlim([0,5000])



figure(2)
subplot(2,2,1)
plot(refSig)
title('\fontsize{20}Reference Single Neuron Activity')
xlabel('\fontsize{20}Time Step')
ylabel('$a(r_1(t),t)$','Interpreter','LaTeX','Fontsize',20);
ylim([0,2])
xlim([0,5000])

subplot(2,2,3)
plot(optSig)
title('\fontsize{20}Optimized Single Neuron Activity')
xlabel('\fontsize{20}Time Step')
ylabel('$a(o_1(t),t)$','Interpreter','LaTeX','Fontsize',20);
ylim([0,2])
xlim([0,5000])


subplot(2,2,2)
plot(naiveSig)
title('\fontsize{20}Naive Single Neuron Activity')
xlabel('\fontsize{20}Time Step')
ylabel('$a(n_1(t),t)$','Interpreter','LaTeX','Fontsize',20);
ylim([0,2])
xlim([0,5000])


subplot(2,2,4)
plot(controlSig)
title('\fontsize{20}Control Single Neuron Activity')
xlabel('\fontsize{20}Time Step')
ylabel('$a(c_1(t),t)$','Interpreter','LaTeX','Fontsize',20);
ylim([0,2])
xlim([0,5000])

