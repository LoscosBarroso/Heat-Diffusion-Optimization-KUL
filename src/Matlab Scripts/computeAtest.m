
function A = computeAtest(A, K)

n = size(K, 2);

% interior
for k = 2:n-1
    for l = 2:n-1
        ind_row = get1DIndex(n, k, l);
        A(ind_row,get1DIndex(n, k, l-1)) = K(2*k,l-1); % T(i,j-1)*K(i,j-1/2)
        A(ind_row,get1DIndex(n, k, l+1)) = K(2*k,l); % T(i,j+1)*K(i,j+1/2)
        A(ind_row,get1DIndex(n, k-1, l)) = K(2*k-1,l); % T(i-1,j)*K(i-1/2,j)
        A(ind_row,get1DIndex(n, k+1, l)) = K(2*k+1,l); % T(i+1,j)*K(i+1/2,j)
        A(ind_row,get1DIndex(n, k, l)) = -(K(2*k,l-1)+K(2*k,l)+K(2*k-1,l)+...
            K(2*k+1,l)); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i-1/2,j)+K(i+1/2,j))
    end
end

% border left
for k = 1:n
    A(get1DIndex(n, k, 1), get1DIndex(n, k, 1)) = 1;
end

% border right
for k = 1:n
    A(get1DIndex(n, k, n), get1DIndex(n, k, n)) = 1;
end

% border up
for k = 1:n
    A(get1DIndex(n, 1, k), get1DIndex(n, 1, k)) = 1;
end

% border down
for k = 1:n
    A(get1DIndex(n, n, k), get1DIndex(n, n, k)) = 1;
end

end