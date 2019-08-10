function AvT = dervArespvT(v,T,power,k_m,k_p)

mesh_size = sqrt(length(T));
mesh_size2 = 2*mesh_size-1;
p = power;
AvT = sparse(mesh_size*mesh_size,mesh_size2*mesh_size);
k = power * (k_m-k_p);
metal = (2*mesh_size)/5;

vi=0;
while vi < mesh_size-2 % Vs in the top row that dont touch the heat sink
   t1 = vi+1;
   t2 = t1+1;
   [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
   AvT(t1,vi+1) = (T(t2)-T(t1))*k*v(i,j)^(p-1);
   AvT(t2,vi+1) = (T(t1)-T(t2))*k*v(i,j)^(p-1);
   vi = vi + 1;
end
t1 = vi+1;
t2 = t1 + 1;
[i,j] = get2DIndex(mesh_size2, mesh_size, vi);
AvT(t1, vi+1) = (T(t2)-T(t1))*k*v(i,j)^(p-1);
vi = vi + 2; % skip the fictional 

while floor(vi/mesh_size) < 2*metal
    if mod(floor(vi/mesh_size),2) == 0
        while mod(vi,mesh_size) < mesh_size-2 
            t1 = floor(floor(vi/mesh_size)/2)*mesh_size + mod(vi,mesh_size);
            t2 = t1 + 1;
            [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
            AvT(t1+1,vi+1) = (T(t2+1)-T(t1+1))*k*v(i,j)^(p-1);
            AvT(t2+1,vi+1) = (T(t1+1)-T(t2+1))*k*v(i,j)^(p-1);
            vi = vi + 1;
        end
        t1 = floor(floor(vi/mesh_size)/2)*mesh_size + mod(vi,mesh_size);
        t2 = t1 + 1;
        [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
        AvT(t1+1,vi+1) = (T(t2+1)-T(t1+1))*k*v(i,j)^(p-1);  
        vi = vi + 2;
    else
        while mod(vi,mesh_size) < mesh_size-1
            t1 = floor(floor(vi/mesh_size)/2)*mesh_size + mod(vi,mesh_size);
            t2 = t1 + mesh_size;
            [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
            AvT(t1+1,vi+1) = (T(t2+1)-T(t1+1))*k*v(i,j)^(p-1);
            AvT(t2+1,vi+1) = (T(t1+1)-T(t2+1))*k*v(i,j)^(p-1);         
            vi = vi + 1;
        end
        if floor(vi/mesh_size) == 2*metal-1
            t1 = floor(floor(vi/mesh_size)/2)*mesh_size + mod(vi,mesh_size);
            t2 = t1 + mesh_size;
            [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
            AvT(t2+1,vi+1) = ((T(t1+1)-T(t2+1))*k*v(i,j)^(p-1))/2;         
        end
        vi = vi + 1;
    end
end

while floor(vi/mesh_size)<mesh_size2-1
    if mod(floor(vi/mesh_size),2)==0
        while mod(vi,mesh_size) < mesh_size-1
           t1 = floor(floor(vi/mesh_size)/2)*mesh_size+mod(vi,mesh_size);
           t2 = t1 + 1;
           [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
           AvT(t1+1,vi+1) = (T(t2+1)-T(t1+1))*k*v(i,j)^(p-1);
           AvT(t2+1,vi+1) = (T(t1+1)-T(t2+1))*k*v(i,j)^(p-1);
           vi = vi + 1;
        end
        vi = vi + 1;
    else 
       while mod(vi,mesh_size)<mesh_size-1
           t1 = floor(floor(vi/mesh_size)/2)*mesh_size+mod(vi,mesh_size);
           t2 = t1 + mesh_size;
           [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
           AvT(t1+1,vi+1) = (T(t2+1)-T(t1+1))*k*v(i,j)^(p-1);
           AvT(t2+1,vi+1) = (T(t1+1)-T(t2+1))*k*v(i,j)^(p-1);
           vi = vi + 1;
       end
           t1 = floor(floor(vi/mesh_size)/2)*mesh_size+mod(vi,mesh_size);
           t2 = t1 + mesh_size;
           [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
           AvT(t1+1,vi+1) = ((T(t2+1)-T(t1+1))*k*v(i,j)^(p-1))/2;
           AvT(t2+1,vi+1) = ((T(t1+1)-T(t2+1))*k*v(i,j)^(p-1))/2;
           vi = vi + 1;
    end
end

while mod(vi,mesh_size) < mesh_size-1
   t1 = floor(floor(vi/mesh_size)/2)*mesh_size+mod(vi,mesh_size);
   t2 = t1 + 1;
   [i,j] = get2DIndex(mesh_size2, mesh_size, vi);
   AvT(t1+1,vi+1) = ((T(t2+1)-T(t1+1))*k*v(i,j)^(p-1))/2;
   AvT(t2+1,vi+1) = ((T(t1+1)-T(t2+1))*k*v(i,j)^(p-1))/2;
   vi = vi + 1;    
end
vi = vi + 1;
        
end

