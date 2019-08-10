function A = computeA(A, K)

n = size(K, 2);
nr = size(K, 1);

% interior
for k = 2:n-1
    for l = 2:n-1
        ind_row = get1DIndex(n, k, l);
        A(ind_row,ind_row-1) = K(2*(k-1)+1,l-1); % T(i,j-1)*K(i,j-1/2)
        A(ind_row,ind_row+1) = K(2*(k-1)+1,l); % T(i,j+1)*K(i,j+1/2)
        A(ind_row,ind_row-n) = K(2*(k-1),l); % T(i-1,j)*K(i-1/2,j)
        A(ind_row,ind_row+n) = K(2*(k-1)+2,l); % T(i+1,j)*K(i+1/2,j)
        A(ind_row,ind_row) = -(K(2*(k-1)+1,l-1)+K(2*(k-1)+1,l)+K(2*(k-1),l)+...
            K(2*(k-1)+2,l)); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i-1/2,j)+K(i+1/2,j))
    end
end


% border left
for k = 2:n-1
    ind_row = get1DIndex(n, k, 1);
    A(ind_row,ind_row-n) = K(2*(k-1),1); % T(i-1,j)*K(i-1/2,j)
    A(ind_row,ind_row+n) = K(2*(k-1)+2,1); % T(i+1,j)*K(i+1/2,j)
    A(ind_row,ind_row+1) = K(2*(k-1)+1,1); % T(i,j+1)*K(i,j+1/2)
    A(ind_row,ind_row) = -(K(2*(k-1),1)+K(2*(k-1)+2,1)+K(2*(k-1)+1,1)); % -T(i,j)*(K(i-1/2,j)+K(i+1/2,j)+K(i,j+1/2))
end

% border right
for k = floor(n*2/5)+1:n-1
    ind_row = get1DIndex(n, k, n);
    A(ind_row,ind_row-n) = K(2*(k-1),n); % T(i-1,j)*K(i-1/2,j)
    A(ind_row,ind_row+n) = K(2*(k-1)+2,n); % T(i+1,j)*K(i+1/2,j)
    A(ind_row,ind_row-1) = K(2*(k-1)+1,n-1); % T(i,j+1)*K(i,j-1/2)
    A(ind_row,get1DIndex(n, k, n)) = -(K(2*(k-1),n)+K(2*(k-1)+2,n)+K(2*(k-1)+1,n-1)); % -T(i,j)*(K(i-1/2,j)+K(i+1/2,j)+K(i,j-1/2))
end

% border up
for k = 2:n-1
    ind_row = get1DIndex(n, 1, k);
    A(ind_row,get1DIndex(n, 1, k-1)) = K(1,k-1); % T(i,j-1)*K(i,j-1/2)
    A(ind_row,get1DIndex(n, 1, k+1)) = K(1,k); % T(i,j+1)*K(i,j+1/2)
    A(ind_row,get1DIndex(n, 2, k)) = K(2,k); % T(i+1,j)*K(i+1/2,j)
    A(ind_row,get1DIndex(n, 1, k)) = -(K(1,k-1)+K(1,k)+...
        K(2,k)); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i+1/2,j))
end

% border down
for k = 2:n-1
    ind_row = get1DIndex(n, n, k);
    A(ind_row,get1DIndex(n, n, k-1)) = K(2*n-1,k-1); % T(i,j-1)*K(i,j-1/2)
    A(ind_row,get1DIndex(n, n, k+1)) = K(2*n-1,k); % T(i,j+1)*K(i,j+1/2)
    A(ind_row,get1DIndex(n, n-1, k)) = K(2*n-2,k); % T(i-1,j)*K(i-1/2,j)
    A(ind_row,get1DIndex(n, n, k)) = -(K(2*n-1,k-1)+K(2*n-1,k)+...
        K(2*n-2,k)); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i-1/2,j))
end

% near heat sink
for k = 1:floor(n*2/5)
    ind = get1DIndex(n, k, n);
    A(ind, ind) = 1;
end

% corners
ind_row = get1DIndex(n, 1, 1);
A(ind_row,get1DIndex(n, 2, 1)) = K(2,1); 
A(ind_row,get1DIndex(n, 1, 2)) = K(1,1); 
A(ind_row,get1DIndex(n, 1, 1)) = -(K(2,1)+K(1,1));

ind_row = get1DIndex(n, n, 1);
A(ind_row,ind_row-n) = K(2*n-2,1); 
A(ind_row,ind_row+1) = K(2*n-1,1); 
A(ind_row,ind_row) = -(K(2*n-2,1)+K(2*n-1,1));

ind_row = get1DIndex(n, n, n);
A(ind_row,ind_row) = -(K(2*n-2,n)+K(2*n-1,n-1));
A(ind_row,ind_row-n) = K(2*n-2,n); 
A(ind_row,ind_row-1) = K(2*n-1,n-1); 


end