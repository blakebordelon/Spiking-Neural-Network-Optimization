function isiDist(path)

ref = csvread(strcat(path, '\refSpikes.csv'));
opt = csvread(path + '\newOptSpikes.csv');
naive = csvread(path + '\naiveSpikes.csv');
control = csvread(path + '\controlSpikes.csv');


isi_ref = [];
isi_opt = [];
isi_naive = [];
isi_control = [];

cum_time_ref = 0;
cum_time_naive = 0;
cum_time_opt = 0;
cum_time_control = 0;

for n=1:1:size(ref,1)
    for t=1:1:size(ref,2)
        if ref(n,t) == 1
            isi_ref = [isi_ref, cum_time_ref];
            cum_time_ref = 0;
        
        else
            cum_time_ref = cum_time_ref + 1;
        end
        if naive(n,t) == 1
            isi_naive = [isi_naive, cum_time_naive];
            cum_time_naive = 0;
        
        else
            cum_time_time_naive = cum_time_naive + 1;
        end
        if opt(n,t) == 1
            isi_opt = [isi_opt, cum_time_opt];
            cum_time_opt = 0;
        
        else
            cum_time_opt = cum_time_opt + 1;
        end
        if control(n,t) == 1
            isi_control = [isi_control, cum_time_ref];
            cum_time_ref = 0;
        
        else
            cum_time_control = cum_time_control + 1;
        end
    end
end


max_isi_count = max( [sum(isi_ref), sum(isi_naive), sum(isi_opt), sum(isi_control)]);
max_isi = max(isi_ref);

edges = linspace(1,max_isi, max_isi_count / 80);
hr = histogram(isi_ref, edges).Values;
hn = histogram(isi_naive, edges).Values;
ho = histogram(isi_opt, edges).Values;
hc = histogram(isi_control, edges).Values;


l2n = norm(hr - hn);
l2o = norm(hr - ho);
l2c = norm(hr - hc);

l2_vals = [l2n, l2o, l2c];
csvwrite(path + '\L2_vals.csv', l2_vals);


ref_count = sum(ref,2);
naive_count = sum(naive,2);
opt_count = sum(opt,2);
control_count = sum(control,2);

mean_vals = [mean(ref_count), mean(naive_count), mean(opt_count), mean(control_count)];
std_vals = [std(ref_count), std(naive_count), std(opt_count), std(control_count)];

csvwrite(path + '\mean_vals.csv', mean_vals);
csvwrite(path + '\std_vals.csv', std_vals);

end