function varout = vRossum(path, plot)

path
ref = csvread(strcat(path, '\refSpikes.csv'));
opt = csvread(path + '\newOptSpikes.csv');
naive = csvread(path + '\naiveSpikes.csv');
control = csvread(path + '\controlSpikes.csv');

x = linspace(-20,20);
filter = exp(- .005 * x.^2);
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
    
    
    
    refSig  = conv(filter, refTrain);
    optSig = conv(filter, optTrain);
    naiveSig = conv(filter, naiveTrain);
    controlSig = conv(filter, controlTrain);

    
    totalRef = totalRef + refSig;
    totalOpt = totalOpt + optSig;
    totalNaive = totalNaive + naiveSig;
    totalControl = totalControl + controlSig;
    
    opterror = opterror + mean((refSig - optSig).^2);
    naiveerror = naiveerror + mean((refSig - naiveSig).^2);
    controlerror = controlerror + mean((refSig - controlSig).^2);

end




opterror
naiveerror
controlerror

totalNaiveError = 1/size(ref,1) * mean((totalRef - totalNaive).^2);
totalOptError = 1/size(ref,1) * mean((totalRef - totalOpt).^2);
totalControlError = 1/size(ref,1) * mean((totalRef - totalControl).^2);


localerrors = [naiveerror, opterror,controlerror];
totalerrors = [totalNaiveError, totalOptError,totalControlError];

csvwrite(path + '\localVanRossum.csv', localerrors);
csvwrite(path + '\totalVanRossum.csv', totalerrors);


if plot == 1
    
end

varout = 0
