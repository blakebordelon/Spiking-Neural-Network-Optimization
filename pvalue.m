function [P_value alpha ks] = pvalue(burst)

[alpha xmin] = tplfit(burst,40) 
% alpha = plmle(burst, 'xmin', xmin)
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

KS = max(abs(cdf - fit))

ks = [];        j = 1;

h = figure; semilogx(xmin:xmax,fit,'*'); hold on; semilogx(xmin:xmax,cdf,'+r'); xlabel('Size','Fontsize',16);   ylabel('Prob(size < S)','Fontsize',16);
title('CDF,field view # ','Fontsize',17);
A = lognfit(z); mu = A(1); sig = A(2);
% plot( ( erfc( (log(xmin)-mu)/sig/sqrt(2) ) - erfc( (log(s+1)-mu)/sig/sqrt(2) ) ) ...
%     ...
%                                                  / ( erfc( (log(xmin)-mu)/sig/sqrt(2) ) - erfc( (log(xmax)-mu)/sig/sqrt(2) ) ),'g' );
%                                              
% plot( 1 - erfc( (log(s)-mu)/sqrt(2)/sig) / erfc( (log(xmin)-mu)/sqrt(2)/sig) ,'m');                                            
                                             
legend('Power law CDF','Experimental CDF','Location','Southeast'); legend boxoff;
% text(18, 0.6, ['\alpha',' = ',num2str(alpha) , '\pm', num2str(SD)] ,'Fontsize',17 );
% pubgraph(h,14,2,'w')

cdfplot(z); grid off

while j<200 
    j
%% ########################## Inverse Method  #############################################
   
syn_data = floor( (xmin-1/2)*(1-rand(1,N)).^(1/(1-alpha)) + 1/2 ) ; 
syn_data = floor(heaviside(xmax-syn_data)) .* syn_data;
% syn_data(syn_data<xmin)=[];  

%% ########################## Accept-Reject Method  #############################################
% beta = alpha - 1;
% umax = xmin^-beta;      u = umax * rand(1,N);       r = floor(u.^-(beta^-1) );
% 
% syn_data = r .* floor(heaviside( (P(r,xmin,alpha) .* Q(xmin,xmin,beta))./(P(xmin,xmin,alpha) .* Q(r,xmin,beta)   ) - rand(1,N) ));
% syn_data = floor(heaviside(xmax-syn_data)) .* syn_data;

%% ############################################################################################
    X = syn_data;
    % [junk xmin]=plfit(burst);
    n = length(syn_data(xmin<=syn_data & syn_data<=xmax)) ;

%     s = unique(X(xmin<=X & X<=xmax));
    s = xmin:xmax;
    LL = @(x) x*sum( log( X(xmin<=X & X<=xmax) ) ) - n*log( 1/sum(s.^-x ) ) ;
    [a,fval] = fminsearch(LL , 2.3);
    
    if abs(a-alpha)<=0.1

        X = X(X>=xmin);            n    = length(X);
        cdf = cumsum(hist(X,xmin:xmax)./n);

        s = unique(X(xmin<=X & X<=xmax));
        % s = unique(X(xmin<=X & X<=100*landa));    s = unique(X(xmin<=X ));
        smin = min(s);      smax = max(s);
%         A = 1/ sum(s.^-alpha );
%         fit = cumsum ( A*(xmin:xmax).^-alpha );
        % ks = max(abs(cdf - fit));
        ks = [ks , max(abs(cdf - fit)) ];
        j = j + 1;
    end
end

P_value = sum( sign( ks(ks>=KS) ) )/200;

display(['P_value = ',num2str(P_value),])