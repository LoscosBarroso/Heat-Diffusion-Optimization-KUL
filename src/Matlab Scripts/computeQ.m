
function Q = computeQ(n, q, T_s)

Q = q * ones(n * n, 1);

for k = 1:floor(n*2/5)
    Q(get1DIndex(n, k, n)) = T_s;
end


end