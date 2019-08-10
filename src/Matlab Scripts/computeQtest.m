
function Q = computeQtest(n, q, T_s)

Q = q * ones(n * n, 1);

for k = 1:n
    Q(get1DIndex(n, k, n)) = T_s;
end

for k = 1:n
    Q(get1DIndex(n, n, k)) = T_s;
end

for k = 1:n
    Q(get1DIndex(n, 1, k)) = T_s;
end

for k = 1:n
    Q(get1DIndex(n, k, 1)) = T_s;
end

end