
function [] = plotV(v)

nc = size(v,2);
nr = size(v,1);
triSize = 1/1/(nc-1);

% top triangles
for i=1:nc-1
    c = 1-v(1,i);
    patch([(i-1)*triSize triSize/2+(i-1)*triSize triSize+(i-1)*triSize], ... 
        [1 1-triSize/2 1], [c c c]);
end

% bot triangles
for i=1:nc-1
    c = 1-v(nr,i);
    patch([(i-1)*triSize triSize/2+(i-1)*triSize triSize+(i-1)*triSize], ... 
        [0 triSize/2 0], [c c c]);
end

% right triangles
for i=1:nc-1
    c = 1-v(nr-(2+(i-1)*2)+1,nc);
    patch([1 1-triSize/2 1], [(i-1)*triSize triSize/2+(i-1)*triSize ...
        triSize+(i-1)*triSize], [c c c]);
end

% left triangles
for i=1:nc-1
    c = 1-v(nr-(+2+(i-1)*2)+1,1);
    patch([0 triSize/2 0], [(i-1)*triSize triSize/2+(i-1)*triSize ...
        triSize+(i-1)*triSize], [c c c]);
end

% less squared rows
for j=nr-1:-2:2
    for i=1:nc-2
        c = 1-v(j,i+1);
        patch([(i-1)*triSize+triSize/2 (i-1)*triSize+triSize ...
            (i-1)*triSize+3*triSize/2 (i-1)*triSize+triSize], ...
            [1-1/(nr-1)*j+triSize/2 1-1/(nr-1)*j+0 ...
            1-1/(nr-1)*j+triSize/2 1-1/(nr-1)*j+triSize], [c c c]);
    end
end

% more squared rows
for j=nr-2:-2:3
    for i=1:nc-1
        c = 1-v(j,i);
        patch([(i-1)*triSize (i-1)*triSize+triSize/2 ...
            (i-1)*triSize-triSize/2+3*triSize/2 (i-1)*triSize+triSize-triSize/2], ...
            [1-1/(nr-1)*j+triSize/2 1-1/(nr-1)*j+0 ...
            1-1/(nr-1)*j+triSize/2 1-1/(nr-1)*j+triSize], [c c c]);
    end
end

end

