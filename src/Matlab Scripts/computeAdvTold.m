function A = computeAdvTold(v, k_p, k_m,ii,jj,power)



nc = size(v, 2);
nr = size(v, 1);
n = nc;
A = sparse(nc ^ 2, nc ^ 2);
m = zeros(nr, nc);
m(ii,jj) = 1; % make sure we do not write uselessly long code
s = (k_m - k_p)/2;
p=power;

% interior
for k = 2:n-1
    for l = 2:n-1
        ind_row = get1DIndex(n, k, l);
        A(ind_row,ind_row-1) = 2*s*p*m(2*(k-1)+1,l-1)*v(2*(k-1)+1,l-1)^(p-1); % T(i,j-1)*K(i,j-1/2)
        A(ind_row,ind_row+1) = 2*s*p*m(2*(k-1)+1,l)*v(2*(k-1)+1,l)^(p-1); % T(i,j+1)*K(i,j+1/2)
        A(ind_row,ind_row-n) = 2*s*p*m(2*(k-1),l)*v(2*(k-1),l)^(p-1); % T(i-1,j)*K(i-1/2,j)
        A(ind_row,ind_row+n) = 2*s*p*m(2*(k-1)+2,l)*v(2*(k-1)+2,l)^(p-1); % T(i+1,j)*K(i+1/2,j)
        A(ind_row,ind_row) = -(2*s*p*m(2*(k-1)+1,l-1)*v(2*(k-1)+1,l-1)^(p-1)+...
            2*s*p*m(2*(k-1)+1,l)*v(2*(k-1)+1,l)^(p-1)+2*s*p*m(2*(k-1),l)*v(2*(k-1),l)^(p-1)+...
            2*s*p*m(2*(k-1)+2,l)*v(2*(k-1)+2,l)^(p-1)); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i-1/2,j)+K(i+1/2,j))
    end
end

% border left
for k = 2:n-1
    ind_row = get1DIndex(n, k, 1);
    A(ind_row,ind_row-n) = 2*s*p*m(2*(k-1),1)*v(2*(k-1),1)^(p-1); % T(i-1,j)*K(i-1/2,j)
    A(ind_row,ind_row+n) = 2*s*p*m(2*(k-1)+2,1)*v(2*(k-1)+2,1)^(p-1); % T(i+1,j)*K(i+1/2,j)
    A(ind_row,ind_row+1) = 2*s*p*m(2*(k-1)+1,1)*v(2*(k-1)+1,1)^(p-1); % T(i,j+1)*K(i,j+1/2)
    A(ind_row,ind_row) = -(2*s*p*m(2*(k-1),1)*v(2*(k-1),1)^(p-1)+...
        2*s*p*m(2*(k-1)+2,1)*v(2*(k-1)+2,1)^(p-1)+2*s*p*m(2*(k-1)+1,1)*v(2*(k-1)+1,1)^(p-1)); % -T(i,j)*(K(i-1/2,j)+K(i+1/2,j)+K(i,j+1/2))
end

% border right
for k = floor(n*2/5)+1:n-1
    if mod(k,2)==2
    ind_row = get1DIndex(n, k, n);
    A(ind_row,ind_row-n) = s*p*m(2*(k-1),n)*v(2*(k-1),n)^(p-1); % T(i-1,j)*K(i-1/2,j)
    A(ind_row,ind_row+n) = s*p*m(2*(k-1)+2,n)*v(2*(k-1)+2,n)^(p-1); % T(i+1,j)*K(i+1/2,j)
    A(ind_row,ind_row-1) = s*p*m(2*(k-1)+1,n-1)*v(2*(k-1)+1,n-1)^(p-1); % T(i,j+1)*K(i,j-1/2)
    A(ind_row,get1DIndex(n, k, n)) = -(s*p*m(2*(k-1),n)*v(2*(k-1),n)^(p-1)+...
        s*p*m(2*(k-1)+2,n)*v(2*(k-1)+2,n)^(p-1)+s*p*m(2*(k-1)+1,n-1)*v(2*(k-1)+1,n-1)^(p-1)); % -T(i,j)*(K(i-1/2,j)+K(i+1/2,j)+K(i,j-1/2))
    end
end

% border up
for k = 2:n-1
    ind_row = get1DIndex(n, 1, k);
    A(ind_row,get1DIndex(n, 1, k-1)) = 2*s*p*m(1,k-1)*v(1,k-1)^(p-1); % T(i,j-1)*K(i,j-1/2)
    A(ind_row,get1DIndex(n, 1, k+1)) = 2*s*p*m(1,k)*v(1,k)^(p-1); % T(i,j+1)*K(i,j+1/2)
    A(ind_row,get1DIndex(n, 2, k)) = 2*s*p*m(2,k)*v(2,k)^(p-1); % T(i+1,j)*K(i+1/2,j)
    A(ind_row,get1DIndex(n, 1, k)) = -(2*s*p*m(1,k-1)*v(1,k-1)^(p-1)+...
        2*s*p*m(1,k)*v(1,k)^(p-1)+2*s*p*m(2,k)*v(2,k)^(p-1)); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i+1/2,j))
end

% border down
for k = 2:n-1
    ind_row = get1DIndex(n, n, k);
    A(ind_row,get1DIndex(n, n, k-1)) = s*p*m(2*n-1,k-1)*v(2*n-1,k-1)^(p-1); % T(i,j-1)*K(i,j-1/2)
    A(ind_row,get1DIndex(n, n, k+1)) = s*p*m(2*n-1,k)*v(2*n-1,k)^(p-1); % T(i,j+1)*K(i,j+1/2)
    A(ind_row,get1DIndex(n, n-1, k)) = s*p*m(2*n-2,k)*v(2*n-2,k)^(p-1); % T(i-1,j)*K(i-1/2,j)
    A(ind_row,get1DIndex(n, n, k)) = -(s*p*m(2*n-1,k-1)*v(2*n-1,k-1)^(p-1)+...
        s*p*m(2*n-1,k)*v(2*n-1,k)^(p-1)+s*p*m(2*n-2,k)*v(2*n-2,k)^(p-1)); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i-1/2,j))
end

% corners
ind_row = get1DIndex(n, 1, 1);
A(ind_row,get1DIndex(n, 2, 1)) = 2*s*p*m(2,1)*v(2,1)^(p-1); 
A(ind_row,get1DIndex(n, 1, 2)) = 2*s*p*m(1,1)*v(1,1)^(p-1); 
A(ind_row,get1DIndex(n, 1, 1)) = -(2*s*p*m(2,1)*v(2,1)^(p-1)+2*s*p*m(1,1)*v(1,1)^(p-1));

ind_row = get1DIndex(n, n, 1);
A(ind_row,ind_row-n) = 2*s*p*m(2*n-2,1)*v(2*n-2,1)^(p-1); 
A(ind_row,ind_row+1) = 2*s*p*m(2*n-1,1)*v(2*n-1,1)^(p-1); 
A(ind_row,ind_row) = -(2*s*p*m(2*n-2,1)*v(2*n-2,1)^(p-1)+2*s*p*m(2*n-1,1)*v(2*n-1,1)^(p-1));

ind_row = get1DIndex(n, n, n);
A(ind_row,ind_row) = -(2*s*p*m(2*n-1,n-1)*v(2*n-1,n-1)^(p-1)+2*s*p*m(2*n-2,n)*v(2*n-2,n)^(p-1));
A(ind_row,ind_row-n) = 2*s*p*m(2*n-2,n)*v(2*n-2,n)^(p-1); 
A(ind_row,ind_row-1) = 2*s*p*m(2*n-1,n-1)*v(2*n-1,n-1)^(p-1); 


% % interior
% for k = 2:n-1
%     for l = 2:n-1
%         ind_row = get1DIndex(n, k, l);
%         A(ind_row,get1DIndex(n, k, l-1)) = m(2*(k-1)+1,l-1)*v(2*(k-1)+1,l-1)^(p-1)*2*s*p; % T(i,j-1)*K(i,j-1/2)
%         A(ind_row,get1DIndex(n, k, l+1)) = m(2*(k-1)+1,l)*v(2*(k-1)+1,l)^(p-1)*2*s*p; % T(i,j+1)*K(i,j+1/2)
%         A(ind_row,get1DIndex(n, k-1, l)) = m(2*(k-1),l)*v(2*(k-1),l)^(p-1)*2*s*p; % T(i-1,j)*K(i-1/2,j)
%         A(ind_row,get1DIndex(n, k+1, l)) = m(2*(k-1)+2,l)*v(2*(k-1)+2,l)^(p-1)*2*s*p; % T(i+1,j)*K(i+1/2,j)
%         A(ind_row,get1DIndex(n, k, l)) = -(m(2*(k-1)+1,l-1)*v(2*(k-1)+1,l-1)^(p-1)*2*s*p+...
%             m(2*(k-1)+1,l)*v(2*(k-1)+1,l)^(p-1)*2*s*p+m(2*(k-1),l)*v(2*(k-1),l)^(p-1)*2*s*p+...
%             m(2*(k-1)+2,l)*v(2*(k-1)+2,l)^(p-1)*2*s*p); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i-1/2,j)+K(i+1/2,j))
%     end
% end
% 
% % border left
% for k = 2:n-1
%     ind_row = get1DIndex(n, k, 1);
%     A(ind_row,ind_row-n) = m(2*(k-1),1)*2*v(2*(k-1),1)^(p-1)*s*p; % T(i-1,j)*K(i-1/2,j)
%     A(ind_row,ind_row+n) = 2*s*m(2*(k-1)+2,1)*v(2*(k-1)+2,1)^(p-1)*p; % T(i+1,j)*K(i+1/2,j)
%     A(ind_row,ind_row+1) = 2*s*m(2*(k-1)+1,1)*v(2*(k-1)+1,1)^(p-1)*p; % T(i,j+1)*K(i,j+1/2)
%     A(ind_row,ind_row) = -(m(2*(k-1),1)*2*v(2*(k-1),1)^(p-1)*s*p+...
%         2*s*m(2*(k-1)+2,1)*v(2*(k-1)+2,1)^(p-1)*p+2*s*m(2*(k-1)+1,1)*v(2*(k-1)+1,1)^(p-1)*p); % -T(i,j)*(K(i-1/2,j)+K(i+1/2,j)+K(i,j+1/2))
% end
% 
% % border right
% for k = floor(n*2/5)+1:n-1
%     ind_row = get1DIndex(n, k, n);
%     A(ind_row,ind_row-n) = 2*s*m(2*(k-1),n)*v(2*(k-1),n)^(p-1)*p; % T(i-1,j)*K(i-1/2,j)
%     A(ind_row,ind_row+n) = 2*s*m(2*(k-1)+2,n)*v(2*(k-1)+2,n)^(p-1)*p; % T(i+1,j)*K(i+1/2,j)
%     A(ind_row,ind_row-1) = 2*s*m(2*(k-1)+1,n-1)*v(2*(k-1)+1,n-1)^(p-1)*p; % T(i,j+1)*K(i,j-1/2)
%     A(ind_row,get1DIndex(n, k, n)) = -(2*s*m(2*(k-1),n)*v(2*(k-1),n)^(p-1)*p+...
%         2*s*m(2*(k-1)+2,n)*v(2*(k-1)+2,n)^(p-1)*p+2*s*m(2*(k-1)+1,n-1)*v(2*(k-1)+1,n-1)^(p-1)*p); % -T(i,j)*(K(i-1/2,j)+K(i+1/2,j)+K(i,j-1/2))
% end
% 
% % border up
% for k = 2:n-1
%     ind_row = get1DIndex(n, 1, k);
%     A(ind_row,get1DIndex(n, 1, k-1)) = 2*s*m(1,k-1)*v(1,k-1)^(p-1)*p; % T(i,j-1)*K(i,j-1/2)
%     A(ind_row,get1DIndex(n, 1, k+1)) = 2*s*m(1,k)*v(1,k)^(p-1)*p; % T(i,j+1)*K(i,j+1/2)
%     A(ind_row,get1DIndex(n, 2, k)) = 2*s*m(2,k)*v(2,k)^(p-1)*p; % T(i+1,j)*K(i+1/2,j)
%     A(ind_row,get1DIndex(n, 1, k)) = -(2*s*m(1,k-1)*v(1,k-1)^(p-1)*p+...
%         2*s*m(1,k)*v(1,k)^(p-1)*p+2*s*m(2,k)*v(2,k)^(p-1)*p); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i+1/2,j))
% end
% 
% % border down
% for k = 2:n-1
%     ind_row = get1DIndex(n, n, k);
%     A(ind_row,get1DIndex(n, n, k-1)) = 2*s*m(2*n-1,k-1)*v(2*n-1,k-1)^(p-1)*p; % T(i,j-1)*K(i,j-1/2)
%     A(ind_row,get1DIndex(n, n, k+1)) = 2*s*m(2*n-1,k)*v(2*n-1,k)^(p-1)*p; % T(i,j+1)*K(i,j+1/2)
%     A(ind_row,get1DIndex(n, n-1, k)) = 2*s*m(2*n-2,k)*v(2*n-2,k)^(p-1)*p; % T(i-1,j)*K(i-1/2,j)
%     A(ind_row,get1DIndex(n, n, k)) = -(2*s*m(2*n-1,k-1)*v(2*n-1,k-1)^(p-1)*p+...
%         2*s*m(2*n-1,k)*v(2*n-1,k)^(p-1)*p+2*s*m(2*n-2,k)*v(2*n-2,k)^(p-1)*p); % -T(i,j)*(K(i,j-1/2)+K(i,j+1/2)+K(i-1/2,j))
% end
% 
% 
% % corners
% ind_row = get1DIndex(n, 1, 1);
% A(ind_row,get1DIndex(n, 2, 1)) = 2*s*m(2,1)*v(2,1)^(p-1)*p; 
% A(ind_row,get1DIndex(n, 1, 2)) = 2*s*m(1,1)*v(1,1)^(p-1)*p; 
% A(ind_row,get1DIndex(n, 1, 1)) = -(2*s*m(2,1)*v(2,1)^(p-1)+2*s*m(1,1)*v(1,1)^(p-1))*p;
% 
% ind_row = get1DIndex(n, n, 1);
% A(ind_row,get1DIndex(n, n-1, 1)) = 2*s*m(2*n-2,1)*v(2*n-2,1)^(p-1)*p; 
% A(ind_row,get1DIndex(n, n, 2)) = 2*s*m(2*n-1,1)*v(2*n-1,1)^(p-1)*p; 
% A(ind_row,get1DIndex(n, n, 1)) = -(2*s*m(2*n-2,1)*v(2*n-2,1)^(p-1)*p+2*s*m(2*n-1,1)*v(2*n-1,1)^(p-1)*p);
% 
% ind_row = get1DIndex(n, n, n);
% A(ind_row,ind_row-n) = 2*s*m(2*n-2,n)*v(2*n-2,n)^(p-1)*p; 
% A(ind_row,ind_row-1) = 2*s*m(2*n-1,n-1)*v(2*n-1,n-1)^(p-1)*p; 
% A(ind_row,ind_row) = -(2*s*m(2*n-2,n)*v(2*n-2,n)^(p-1)*p+2*s*m(2*n-1,n-1)*v(2*n-1,n-1)^(p-1)*p);




end