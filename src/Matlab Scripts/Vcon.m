function con = Vcon(v)

mesh_size = size(v, 2);
mesh_size2 = 2*mesh_size-1;
con = 0;

i = 1;
for j = 1:mesh_size-1
    con = con + v(i, j);
end

for i = 2:mesh_size2-1
    if mod(i,2) == 0
        for j = 1:mesh_size-1
            con = con + v(i,j);
        end
        con = con + v(i,mesh_size)/2;
    else
        for j = 1:mesh_size-1
            con = con + v(i,j);
        end
    end
end

i = mesh_size2;
for j = 1:mesh_size-1
    con = con + v(i, j)/2;
end

con = con / ((mesh_size-1)*(mesh_size2-1)+0.5*(mesh_size+mesh_size2));

end

