function [alpha xmin ks Loglike] = tplfit(burst,limit, x_start, vargarin)

KS = []; alpha = [];  Loglike = [];

switch nargin
    case 2
        range = 1:limit;
    case 3        
            range = x_start;        
end
        

for x0 = range

% [alpha xmin]=plfit(burst);  
X = burst;
% xmax = size(data,1);
% xmin = 5;
xmax = max(burst);
n = length(burst(x0<=burst & burst<=xmax)) ;
s = unique(X(x0<=X & X<=xmax));
smin = min(s);      
smax = max(s);
LL = @(x) x*sum( log( burst(x0<=burst & burst<=xmax) ) ) - n*log( 1/sum(s.^-x ) ) ;

[a,fval] = fminsearch(LL , 2.3 );

Loglike = [Loglike , -fval];


% X = burst;
% 
% for i = 1:1000
%     
%     X = Shuffle(X);
%     I = randi(length(X),[1 length(X)]);    
%     AV = X(I);
%     [sth xmin]=plfit(AV);
% %     xmax = max(AV);
%     n = length(AV(xmin<=AV & AV<=xmax)) ;
%     s = unique(AV(xmin<=AV & AV<=xmax));
%     smin = min(s);      smax = max(s);
%     LL = @(x) x*sum( log( AV(xmin<=AV & AV<=xmax) ) ) - n*log( 1/sum(s.^-x ) ) ;
% %     options = optimset('PlotFcns',@optimplotfval);
%     [a2,fval] = fminsearch(LL , 2.8);
%     alpha = [alpha , a2];
%     clear a2 I
% end
% 
% SD = std(alpha);
% 
% 
% 
% 
%% ########################## Hypothesis Testing ###############################
% alpha = a;  
% [sth xmin] = plfit(burst);
xmax = max(burst);

N = length(burst);
% n = length(burst);
k = 0;

z = burst;
z = z(z>=x0);            n    = length(z);
cdf = cumsum(hist(z,x0:xmax)./n);

s = unique(X(x0<=X & X<=xmax));
% s = unique(X(xmin<=X & X<=100*landa));    s = unique(X(xmin<=X ));
smin = min(s);      smax = max(s);
A = 1/ sum(s.^-a );
alpha = [alpha , a];
fit{x0} = cumsum ( A*(x0:xmax).^-a );
KS = [KS , max(abs(cdf - fit{x0})) ] ;
% 
% jj = 1;
% ks = [];
% 
% while jj < 5000 
%     
% %% ########################## Inverse Method  #############################################
%    
% syn_data = floor( (xmin-1/2)*(1-rand(1,N)).^(1/(1-alpha)) + 1/2 ) ; 
% syn_data = floor(heaviside(xmax-syn_data)) .* syn_data;
%     
% 
% %% ########################## Accept-Reject Method  #############################################
% % beta = alpha - 1;
% % umax = xmin^-beta;      u = umax * rand(1,N);       r = floor(u.^-(beta^-1) );
% % 
% % syn_data = r .* floor(heaviside( (P(r,xmin,alpha) .* Q(xmin,xmin,beta))./(P(xmin,xmin,alpha) .* Q(r,xmin,beta)   ) - rand(1,N) ));
% % syn_data = floor(heaviside(xmax-syn_data)) .* syn_data;
% 
% %% ############################################################################################
%     X = syn_data;
%     % [junk xmin]=plfit(burst);
%     n = length(syn_data(xmin<=syn_data & syn_data<=xmax)) ;
% 
%     s = unique(X(xmin<=X & X<=xmax));
%     LL = @(x) x*sum( log( X(xmin<=X & X<=xmax) ) ) - n*log( 1/sum(s.^-x ) ) ;
%     [a,fval] = fminsearch(LL , 2.3);
%     
%     if abs(a-alpha)<=0.05
% 
%         X = X(X>=xmin);            n    = length(X);
%         cdf = cumsum(hist(X,xmin:xmax)./n);
% 
%         s = unique(X(xmin<=X & X<=xmax));
%         % s = unique(X(xmin<=X & X<=100*landa));    s = unique(X(xmin<=X ));
%         smin = min(s);      smax = max(s);
%         A = 1/ sum(s.^-alpha );
%         fit = cumsum ( A*(xmin:xmax).^-alpha );
%         % ks = max(abs(cdf - fit));
%         ks = [ks , max(abs(cdf - fit)) ];
%         jj = jj + 1;
%     end
% end
% 
% P_value = sum( sign( ks(ks>=KS) ) )/5000;

end

switch nargin
    case 2
        xmin = find(KS==min(KS)) ;  alpha = alpha( xmin );  Loglike = Loglike(xmin);
    case 3        
            xmin = x_start;        
end


ks = min(KS);

% alpha = a
% Loglike = -fval
% KS = KS
% p = P_value
% ks = KS
% sd = SD






