
% make raster plots from naive half max expt
path = 'naive half max';
rspikes = csvread(strcat(path, '\trial 1\refSpikes.csv'));
ospikes = csvread(strcat(path, '\trial 1\newOptSpikes.csv'));
nspikes = csvread(strcat(path, '\trial 1\naiveSpikes.csv'));
cspikes = csvread(strcat(path, '\trial 1\controlSpikes.csv'));

N = 100;
T = 170;
rspikes = rspikes(1:N, 100:T+100);
ospikes = ospikes(1:N, 100:T+100);
nspikes = nspikes(1:N, 100:T+100);
cspikes = cspikes(1:N, 100:T+100);

fig1 = figure(1);
axis on
set(fig1,'Units','inches','Position',[1,1,8,6])
subplot(2,2,1)
imshow(rspikes)
title('Reference');
xlabel('Time Steps');
axis on
xticks([0:25:T])
yticks([0:25:N])
subplot(2,2,2)
imshow(nspikes)
title('Naive');
xlabel('Time Steps');
axis on
xticks([0:25:T])
yticks([0:25:N])

subplot(2,2,3)
imshow(ospikes)
title('Optimized');
xlabel('Time Steps');
axis on
xticks([0:25:T])
yticks([0:25:N])

subplot(2,2,4)
imshow(cspikes)
xlabel('Time Steps');
title('Control')
axis on
xticks([0:25:T])
yticks([0:25:N])



print(fig1,'-dpdf','spikes.pdf')


% make plots of VR distance

folders = string({'uniform', 'gaussian',  'sparse','naive half max'});

all_local_vals = zeros(4,30,3);
all_total_vals = zeros(4,30,3);
all_l2_ISI = zeros(4,30,3);
all_mean_vals = zeros(4,30,4);
all_std_vals = zeros(4,30,4);

for i=1:1:4
    fname = folders(i);
    for j=1:1:30
        path = fname + '\trial' + ' ' + int2str(j);
        all_local_vals(i,j,:) = csvread(path + '\localVanRossum.csv');
        all_total_vals(i,j,:) = csvread(path + '\totalVanRossum.csv');
        all_l2_ISI(i,j,:) = csvread(path + '\L2_vals.csv');
        all_mean_vals(i,j,:) = csvread(path + '\mean_vals.csv');
        all_std_vals(i,j,:) = csvread(path + '\std_vals.csv');
    end
end

mean_local = mean(all_local_vals, 2)
std_local = std(all_local_vals,0,2)
mean_total = mean(all_total_vals, 2)
std_total = std(all_total_vals, 0,2)
mean_l2 = mean(all_l2_ISI,2)
std_l2 = mean(all_l2_ISI,2)
%mean(all_mean_vals, 2)
%std(all_std_vals,2)




fig2 = figure(2)
for i=1:1:4
    subplot(2,4,2*(i-1) + 1);
    for j = 1:1:30
        all_local_vals(i,j,:);
        plot([1,2,3], reshape(all_local_vals(i,j,:),[1,3]), '-o','Color','black','MarkerFaceColor', 'blue')
        hold on
    end
    ylim([0,200])
    ylabel('$D_P$','Interpreter','latex')
    xticklabels({'N','O','C'})
    xticks([1,2,3])
    % between plots
    hold off
    subplot(2,4,2*i)
    for j = 1:1:30
        plot([1,2,3], reshape(all_total_vals(i,j,:),[1,3]), '-o','Color','black','MarkerFaceColor', 'blue')
        hold on
    end
    ylim([0,140])
    xticklabels({'N','O','C'})
    xticks([1,2,3])
    ylabel('$D_A$','Interpreter','latex')

    hold off
end
set(fig2,'Units','inches','Position',[1,1,8,6])
print(fig2,'-dpdf','van_rossum_dists.pdf')


fig3 = figure(3)

% now you need to do isi, mean spike rate, etc
for i=1:1:4
    subplot(2,6,3*(i-1)+1)
    for j=1:1:30
        plot(reshape(all_l2_ISI(i,j,:),[1,3]), '-o','Color','black','MarkerFaceColor', 'blue');
        hold on
    end
    ylim([0,3000])
    ylabel('ISI $\ell_2$', 'Interpreter', 'latex','fontsize',10)
    xticklabels({'N','O','C'})
    xticks([1,2,3])
    hold off
    
    subplot(2,6,3*(i-1)+2)
    for j=1:1:30
        plot(reshape(all_mean_vals(i,j,:),[1,4]), '-o','Color','black','MarkerFaceColor', 'blue');
        hold on
    end
    ylim([0,75])
    ylabel('$\mu_{sp}$', 'Interpreter', 'latex','fontsize', 12)
    xticklabels({'R','N','O','C'})
    xticks([1,2,3,4])
    hold off
    
    subplot(2,6,3*(i-1)+3)
    for j=1:1:30
        plot(reshape(all_std_vals(i,j,:),[1,4]), '-o','Color','black','MarkerFaceColor', 'blue');
        hold on
    end
    ylim([2,8])
    ylabel('$\sigma_{sp}$', 'Interpreter', 'latex','fontsize',12)
    xticklabels({'R','N','O','C'})
    xticks([1,2,3,4])
    hold off
    
end
set(fig3,'Units','inches','Position',[1,1,8,4])
print(fig3,'-dpdf','isi.pdf')

