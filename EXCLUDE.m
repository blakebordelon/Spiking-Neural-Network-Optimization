function [burstMax, burstMin] = EXCLUDE(burst, varargin)
iVarArg = 1;
KS = 1; setmin = 10;
num = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'setmin',      setmin = varargin{iVarArg+1}; iVarArg = iVarArg + 1
        case 'num',       num = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'flag',   flag = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(AV_analysis) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end

dKS = 1; xmin = 1;
while KS > min(num/sqrt(length(burst(burst>xmin))),0.1) & dKS > 0.0005
    setmin
    [alpha xmin] = tplfit(burst, setmin)
    xmax = max(burst);
    N = length(burst);
    % n = length(burst);
    k = 0;
    z = burst;
    z = z(z>=xmin);            n    = length(z);
    cdf = cumsum(hist(z,xmin:xmax)./n);
    
    s = unique(burst(xmin<=burst & burst<=xmax));
    % s = unique(X(xmin<=X & X<=100*landa));    s = unique(X(xmin<=X ));
    smin = min(s);      %smax = max(s);  s = smin:smax;
    A = 1/ sum(s.^-alpha );
    fit = cumsum ( A*(xmin:xmax).^-alpha );
    KS_old  = KS;
    KS = max(abs(cdf - fit));
    dKS = abs(KS_old - KS);
    
%     plot(cdf,'k');
    hold on
%     plot(fit,'r')
    2/sqrt(length(burst(burst>xmin)))
    burst = burst(burst < max(burst));
    burstMax = max(burst);
end
burstMin = xmin;

