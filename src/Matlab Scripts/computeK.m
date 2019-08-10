
function K = computeK(K, v, k_m, k_p, power)

mesh_size = size(v, 2);
mesh_size2 = size(v, 1);

i = 1;
for j = 1:mesh_size-1
    K(i,j) = k_p + v(i,j)^power * (k_m - k_p);
end

for i = 2:mesh_size2-1
    if mod(i,2) == 0
        for j = 1:mesh_size-1
            K(i,j) = k_p + v(i,j)^power * (k_m - k_p);
        end
        j = mesh_size;
        K(i,j) = (k_p + v(i,j)^power * (k_m - k_p))/2;
    else
        for j = 1:mesh_size-1
            K(i,j) = k_p + v(i,j)^power * (k_m - k_p);
        end
    end
end

i = mesh_size2;
for j = 1:mesh_size-1
    K(i,j) = (k_p + v(i,j)^power * (k_m - k_p))/2;
end

% for i = 1:mesh_size
%     for j = 1:mesh_size2
% %         K(i, j) = v(i, j) * k_m + (1 - v(i, j)) * k_p;
%         K(i,j) = k_p + v(i,j)^power * (k_m - k_p);
% %         if i == n_r || j == n_c || (i==1 
% %            K(i,j) = K(i,j)/2; 
% %         end
%     end   
% end

end