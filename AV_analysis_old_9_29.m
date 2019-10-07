function [Result, burst, Tm_on_off]=AV_analysis(data, varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0, 'defaultaxesfontsize', 18)
set(0, 'defaulttextfontsize', 15)
set(0, 'defaultaxeslinewidth', 2)
color = [0.82 .2 0.2];
iVarArg = 1; flag = 1; perc = 0.35; frmrt = 1;
while iVarArg <= length(varargin)
    argOkay = true;
    switch varargin{iVarArg},
        case 'color',       color = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'chunk',       chunk = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'flag',   flag = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'perc',   perc = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
        case 'frmrt',  frmrt = varargin{iVarArg+1}; iVarArg = iVarArg + 1;
            
        otherwise,
            argOkay = false;
    end
    if ~argOkay
        disp(['(AV_analysis) Ignoring invalid argument #' num2str(iVarArg+1)]);
    end
    iVarArg = iVarArg + 1;
end
%
[n m]=size(data); %get neurons# and frame#, n*m can represent the total time
% data = sign(data);
network=nansum(data);  % define the network
if network(1) == 0;
    network(1) = 1;
end
Sep = find(network~=0);
dSep = diff(Sep);
LSep = find(dSep > 1000);
Remove_loc(1,:) = Sep(LSep)+5;
Remove_loc(2,:) = Sep(LSep+1);
Remove = [];
for i = 1:length(LSep)
    Remove = [Remove, Remove_loc(1,i) : Remove_loc(2,i)];
end
clear i
network(Remove) = [];
m = length(network);
zdata = network;
z2data=zdata;
% figure;
% plot(network)

sortN = sort(network);
threshold2 = sortN(round(m*perc));
zdata=heaviside(zdata-threshold2*ones(size(zdata)));
zdata(find(zdata~=1))=zeros;
Tm_on_off = zdata;
Z = find(zdata==0);
zdata(Z(diff(Z)==1))=[];
z2data(Z(diff(Z)==1))=[];
J = find(zdata==0);
data2 = data; data2(:,Z(diff(Z)==1))=[];

%% #################### Find Spike and AV sizes ########################
% figure(chunk)
burst = [];
clear burstshape
for i = 1:(length(J)-1)
    fired = sum(z2data(J(i)+1: J(i+1)-1));%-threshold2*(J(i+1)-J(i)-1); % Use this line if AVsize = # of Spikes
    %         fired= sum(sign(nansum(data2(:,J(i):J(i+1)),2)));
    burst = [burst ; fired];
    burstshape{i} = z2data(J(i)+1: J(i+1)-1);
end
size(burst)
L = max(burst);
s = 1:1:round(L);     s = s';
[pdf junk] = histc(full(burst),s);
subplot(1, 3,1)
loglog(s, pdf/sum(pdf),'o','Markersize', 5,'Markerfacecolor',[1 1 1],'MarkerEdgecolor',color); hold on
xlabel('S')
ylabel('PDF(S)')
XLm = [1, 10^ceil(log10(L))];log10(max(pdf))/sum(pdf)
YLm = [10^floor(log10(pdf(end)/sum(pdf))), 1]%10^ceil(log10(max(pdf))/sum(pdf))]
% axis(XLm, YLm)
set(gca,'xtick',XLm, 'ytick', YLm, 'xlim', XLm, 'Ylim', YLm)
box off
if flag == 1
    [burstMax, burstMin] = EXCLUDE(burst,'setmin', 300);
    [alpha ,xmin] = plmle(burst(burst < burstMax & burst > 20),'limit',burstMin);
    Result.burst = burst;
    Result.alpha = alpha;
    Result.xmin = xmin;
    Result.xmax = burstMax;
    loglog(burstMin:burstMax, pdf(burstMin:burstMax)/sum(pdf),'o','Markersize', 5,'Markerfacecolor',color,'MarkerEdgecolor',color); hold on
    x = burstMin:burstMax;
    burstMin
    xmin
    y = length(find(burst == burstMin))/(burstMin^-alpha)*x.^-alpha; y = y/sum(pdf);
    loglog(x,y,'k','linewidth',1)
%     title(['AV size distribution, exp = ', num2str(alpha)])
elseif flag == 5
    [burstMax, burstMin] = EXCLUDE(burst,'setmin', 5);
    burstMin
    [alpha ,xmin] = plmle(burst(burst < burstMax & burst > burstMin),'limit',burstMin);
    Result.burst = burst;
    Result.alpha = alpha;
    Result.xmin = xmin;
    Result.xmax = burstMax;
    %     Result.P.burst  = pvaluenew(burst(burst <burstMax));
    loglog(burstMin:burstMax, pdf(burstMin:burstMax)/sum(pdf),'o','Markersize', 5,'Markerfacecolor',color,'MarkerEdgecolor',color); hold on
    
end

%% ######### Duration Distribution ########################################
T=diff(find(zdata==0))-1;
T = T(T>0);
r = 1:max(T);
[tdf junk] = histc(T, r);
subplot(1, 3,2)
loglog(r/frmrt, tdf/sum(tdf),'o','Markersize', 5,'Markerfacecolor',[1 1 1],'MarkerEdgecolor',color); hold on
xlabel('D (s)')
ylabel('PDF(D)')
box off
XLmt = [10^floor(log10(1/frmrt)), 10^ceil(log10(max(T)/frmrt))];
YLmt = [10^floor(log10(tdf(end)/sum(tdf))), 1];%10^ceil(log10(max(tdf))/sum(tdf))]
% axis(XLm, YLm)
set(gca,'xtick',XLmt, 'ytick', YLmt, 'xlim', XLmt, 'Ylim', YLmt)
if flag == 1
    [tMax, tMin] = EXCLUDE(T,'setmin',80);
    [beta tmin] = plmle(T(T < tMax & T > tMin),'limit',tMin);
    Result.T = T;
    Result.beta = beta;
    Result.tmin = tmin;
    Result.tmax = tMax;
    tMin
    %     Result.P.T = pvaluenew(T(T < tMax));
    loglog([tMin: tMax]/frmrt, tdf(tMin: tMax)/sum(tdf),'o','Markersize', 5,'Markerfacecolor',color,'MarkerEdgecolor',color); hold on
    x = tMin:tMax;
    y = length(find(T == tmin))/(tmin^-beta)*x.^-beta; y = y/sum(tdf);
    loglog(x/frmrt,y,'k','linewidth',1)
%     title(['AV duration distribution, exp = ',num2str(beta)])
elseif flag == 5
    [tMax, tMin] = EXCLUDE(T,'setmin', 10);
    [beta tmin] = plmle(T(T < tMax & T > tMin),'limit',tMin);
    Result.T = T;
    Result.beta = beta;
    Result.tmin = tmin;
    Result.tmax = tMax;
    %     Result.P.T = pvaluenew(T(T < tMax));
    loglog([tMin: tMax]/frmrt, tdf(tMin: tMax)/sum(tdf),'o','Markersize', 5,'Markerfacecolor',color,'MarkerEdgecolor',color); hold on
end

%% scaling relation
if flag == 1 || flag ==5
    TT = 1:max(T); Sm = []; clear Smshape
    for i=1:length(TT)
        Sm=[Sm mean(burst(find(T==TT(i))))];
        Smshape{i} = burstshape(find(T==TT(i)));
    end
    Loc=find(Sm>0);TT=TT(Loc);Sm=Sm(Loc);
    fit_sigma = polyfit(log(TT), log(Sm), 1);
    subplot(1, 3,3)
    sigma = (beta - 1)/(alpha - 1);
        loglog(TT/frmrt, Sm,'o','Markersize', 5,'Markerfacecolor',[1 1 1],'MarkerEdgecolor',color);
    tMin
    tmin
    loglog(TT(TT>tMin & TT<tMax)/frmrt, Sm(TT>tMin & TT<tMax),'o','Markersize', 5,'Markerfacecolor',color,'MarkerEdgecolor',color);

    hold on
    loglog(TT/frmrt, Sm(TT(1))*TT.^sigma,'color',[0 0 0],'linewidth',1);
    box off; hold on;
    xlabel('D (s)')
    ylabel('<S>')
    set(gca,'xtick',XLmt, 'ytick', XLm, 'Xlim', XLmt, 'Ylim', XLm)
    title({'Scaling, Block ',['fit = ',num2str(fit_sigma(1)), 'predict = ', num2str(sigma)]});
end
% figure
% for i = 5:50%max(TT)
%     if length(Smshape{i})>15
%         figure(2)
%         plot(mean(full(cell2mat(Smshape{i}'))),'color', i/30*ones(1,5));
%         mx(i) = max(mean(full(cell2mat(Smshape{i}'))));
%         hold on
%         figure(5)
% %         if mod(i,2) == 0
% %             L = [1:i/2,i/2:-1:1];
% %         else
% %             L = [1:i/2, i/2+0.5, i/2-0.5:-1:1];
% %         end
%         ssp = cell2mat(Smshape{i}');
%         ssp = ssp(max(ssp')>0,:);
%         if size(ssp,1) > 20
%         plot(1/i:1/i:1, (mean(full(ssp)-1).^(1/(sigma-1))/(i+2)),'color', (i-5)/50*[1 0 0]);
%         end
%         hold on
%     end
% end
% box off
% xlabel('Normalized T'); ylabel('Normalized S')
% title('Shape collapse')
%% ######################Calculating the exponents###############

