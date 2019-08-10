
function [val, gradient] = myfunc(v)

mesh_size = 32;
mesh_size2 = 2*mesh_size-1;
v = reshape(v, mesh_size2,mesh_size);


plate_width = 5e-3; % name says itself
mesh_width = plate_width / (mesh_size - 1); % name says itself
k_metal = 65; % W*m^-1*K^-1
k_plastic = 0.2; % W*m^-1*K^-1
T_s = 293; % Kelvin
q = - 2e7 * (mesh_width * mesh_width); % heat source as W/m^3
power = 2;


A = sparse(mesh_size ^ 2, mesh_size ^ 2);
K = zeros(2 * mesh_size - 1, mesh_size);
Q = computeQ(mesh_size, q, T_s);
K = computeK(K, v, k_metal, k_plastic, power);
A = computeA(A, K);
T = A\Q;
% % calculating the gradient

AvT = dervArespvT(v,T,power,k_metal,k_plastic);
lambda= transpose(A)\(T-293);
gradient = (-transpose(lambda)*AvT)';
gradient = reshape(gradient, mesh_size, mesh_size2)';

val = 0;
for i=1:length(T)
    val = val + (T(i)-293)^2/2;
end
disp([mean(T) Vcon(v)]);

v = reshape(v, 1,mesh_size2*mesh_size);
gradient = reshape(gradient, 1, mesh_size2*mesh_size)';


end