
folders = string({'gaussian', 'naive half max', 'sparse', 'uniform'});

for i=1:1:4
    fname = folders(i);
    for j=1:1:30
        path = fname + '\trial' + ' ' + int2str(j);
        vRossum(path, 0);
        isiDist(path);
    end
end


