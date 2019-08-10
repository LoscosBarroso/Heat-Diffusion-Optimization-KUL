
function [val, g] = myconstraint(v,a,b)

mesh_size = 32;
mesh_size2 = 2*mesh_size-1;
v = reshape(v, mesh_size2,mesh_size);
g = v;

val = Vcon(v);
    
i = 1;
for j = 1:mesh_size-1
    g(i, j) = 1;
end

for i = 2:mesh_size2-1
    if mod(i,2) == 0
        for j = 1:mesh_size-1
            g(i, j) = 1;
        end
            g(i, mesh_size) = 0.5;
    else
        for j = 1:mesh_size-1
            g(i, j) = 1;
        end
    end
end

i = mesh_size2;
for j = 1:mesh_size-1
    g(i, j) = 0.5;
end

v = reshape(v, 1,mesh_size2*mesh_size);
g = reshape(g, 1,mesh_size2*mesh_size);

    
end