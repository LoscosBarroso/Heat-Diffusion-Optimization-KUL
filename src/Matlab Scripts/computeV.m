function v = computeV(v)

n_r = size(v, 1);
n_c = size(v, 2);

for i = 1:n_r
    for j = 1:n_c
        if mod(i, 2) ~= 0 && j == size(v, 2)
            
        else 
            v(i, j) = 1;
        end
%         v(i, j) = ((i-1)*n_c+j-1)/(n_r*n_c); 
    end   
end

end