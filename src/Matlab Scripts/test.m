% clear;
clc;

% CONSTANTS
mesh_size = 32; % sqrt of the number of points in the meshgrid
mesh_size2 = 2*mesh_size-1;
plate_width = 5e-3; % name says itself
mesh_width = plate_width / (mesh_size - 1); % name says itself
k_metal = 65; % W*m^-1*K^-1
k_plastic = 0.2; % W*m^-1*K^-1
T_s = 293; % Kelvin
q = - 2e7 * (mesh_width * mesh_width); % heat source as W/m^3
power=2;
% ---------------------------------------------
A = sparse(mesh_size ^ 2, mesh_size ^ 2);
K = zeros(2 * mesh_size - 1, mesh_size);
v = zeros(2 * mesh_size - 1, mesh_size);
Q = computeQ(mesh_size, q, T_s);
v = computeV(v);
K = computeK(K, v, k_metal, k_plastic,power);
A = computeA(A, K);
T = A\Q;
T_shaped = reshape(T,mesh_size,mesh_size)';
% plotVT(v);
% figure;
% heatmap(T_shaped);

% % calculating the gradient
% AvT = dervArespvT(v,T,power,k_metal,k_plastic);
% lambda= transpose(A)\(T-293);
% gradient = (-transpose(lambda)*AvT)';
% gradient = reshape(gradient, mesh_size, mesh_size2)';


% optimizing
v = zeros(mesh_size2,mesh_size);
[ff,gg] = myconstraint(v,0,0);

% right triangles
v_con = zeros(mesh_size2,mesh_size);
ind = [];
for i=1:mesh_size-1
    if ((i-1)*2)+1 <= 2*mesh_size2/5
        v_con((i-1)*2+1,mesh_size) = 1;
        ind = [ind (mesh_size-1)*mesh_size2+(i-1)*2+2];
    end
end
cons = zeros(length(ind), mesh_size2*mesh_size);
for i = 1:length(ind)
   cons(i, ind(i)) = 1; 
end

tic
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',1000);
fun = @myfunc;
v0 = zeros(mesh_size*mesh_size2, 1);
Ax = gg;
b = 0.4*(mesh_size-1)*(mesh_size2-1);
Aeq = [];
beq = [];
lb = zeros(1, mesh_size*mesh_size2);
lb(ind) = 1;
ub = ones(1, mesh_size*mesh_size2);
nonlcon = [];
v = fmincon(fun,v0,Ax,b,Aeq,beq,lb,ub,nonlcon,options);
v = reshape(v,mesh_size2,mesh_size);
K = computeK(K, v, k_metal, k_plastic,power);
A = computeA(A, K);
T = A\Q;
T_shaped = reshape(T,mesh_size,mesh_size)';
plotV(v)
toc



 



% % calculating the gradient by the finite differences
% delta = 0.001;
% g = zeros(mesh_size2,mesh_size);
% Av = sparse(mesh_size ^ 2, mesh_size2*mesh_size);
% for i = 1:mesh_size2
%     for j = 1:mesh_size
%     Kn = zeros(2 * mesh_size - 1, mesh_size);
%     Kp = zeros(2 * mesh_size - 1, mesh_size);
%     vp = v;
%     vn = v;
%     vp(i,j) = vp(i,j) + delta;
%     vn(i,j) = vn(i,j) - delta;
%     Kp = computeK(Kp, vp, k_metal, k_plastic,power);
%     Kn = computeK(Kn, vn, k_metal, k_plastic,power);
%     Ap = sparse(mesh_size ^ 2, mesh_size ^ 2);
%     An = sparse(mesh_size ^ 2, mesh_size ^ 2);
%     An = computeA(An, Kn);
%     Ap = computeA(Ap, Kp);
%     Tn = An\Q;
%     Tp = Ap\Q;
%     g(i, j) = (computeCost(Tp)-computeCost(Tn))/(2*delta);
%     A = (Ap-An)./(2*delta);
%     end
% end

% (gradient-gc)/sqrt(sum(sum(gradient.^2)))

