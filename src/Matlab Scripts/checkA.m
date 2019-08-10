function checkA(A)

A = full(A);
is_working = true;
for i = 1:size(A,1)
   if sum(A(i, :)) > 1 || sum(A(i, :)) < -0.1 || (sum(A(i, :)) > 0.1 ...
           && sum(A(i, :)) < 0.99)
       is_working = false;
   end
end

if (is_working) 
   fprintf('A matrix is correct!\n')
else
   fprintf('A matrix is wrong!\n')
end

end

