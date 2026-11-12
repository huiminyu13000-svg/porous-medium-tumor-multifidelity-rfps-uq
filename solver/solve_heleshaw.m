% 基础函数：保留数值更新公式；由 run_experiment 统一调用。
% Adapted from solver_HeleShaw_vivo_modified.m in the supplied archive.
% Numerical model parameters, including the original PME growth factor 0.5, are preserved.
function [rho_mask, p] = solve_heleshaw(k, T_final, dx, dy, dt, ...
                                                       lambda_k, cB_k, G0_k, rho0)
% Convert integer MAT metadata before dividing time by the step size.
T_final=double(T_final);
a = -2.5; b = 2.5;
L = b;
[X, Y] = meshgrid(a:dx:b, a:dx:b);
Nx=length(a:dx:b);
Ny=length(a:dx:b);
Nt = round(T_final / dt);
assert(Nt>=1 && abs(Nt*dt-T_final)<1e-12, 'Positive T must be a multiple of dt.');
lambda = lambda_k;
cB = cB_k;
G0 = G0_k;  %we need to times 2 to recover the boundary loss
z = 0;  % random/heterogeneity control (kept fixed here)
mask0 = double(rho0>0);
phi = signed_distance_from_mask(mask0, dx);
for n=1:Nt
    rho_mask = double(phi < 0);
    rho_mask(rho_mask>0) = 1;  % ensure interior filled with 1 exactly
    c = solve_c_field_vivo(rho_mask, lambda, z, cB, L, Nx, Ny);
    p = solve_pressure_field(rho_mask, G0, c, dx, dy, Nx, Ny);
    p = gaussian_filter(p,1.0);
    [phi_x, phi_y] = gradient(phi, dx, dy);
    norm_grad = sqrt(phi_x.^2 + phi_y.^2 + 1e-12);
    nx = phi_x ./ norm_grad; ny = phi_y ./ norm_grad;
    [p_x, p_y] = gradient(p, dx, dy);
    Vn = -(p_x .* nx + p_y .* ny);
    V_ext = gradient_extension(phi, Vn, dx, dy, 10);
    [phi_x_p,phi_x_m,phi_y_p,phi_y_m] = compute_phi_derivatives(phi, dx, dy);
    H = H_Godunov(phi_x_m, phi_x_p, phi_y_m, phi_y_p, V_ext);
    phi = phi - dt * H;
    if mod(n,20)==0
        phi = redistance(phi, dx, dy, 10);
    end
end
rho_mask = double(phi<0); rho_mask(rho_mask>0)=1;
end
function p = solve_pressure_field(rho, G0, c, dx, dy, Nx, Ny)
N = Nx*Ny; idx_map = reshape(1:N,Ny,Nx);
A = sparse(N,N); b = zeros(N,1);
for iy=1:Ny
    for ix=1:Nx
        row = idx_map(iy,ix);
        if rho(iy,ix)>0
            A(row,row) = 2/dx^2 + 2/dy^2;
            if ix>1, A(row,idx_map(iy,ix-1))=-1/dx^2; end
            if ix<Nx, A(row,idx_map(iy,ix+1))=-1/dx^2; end
            if iy>1, A(row,idx_map(iy-1,ix))=-1/dy^2; end
            if iy<Ny, A(row,idx_map(iy+1,ix))=-1/dy^2; end
            b(row) = G0*c(iy,ix);
        else
            A(row,row)=1; b(row)=0;
        end
    end
end
p = reshape(A\b,[Ny,Nx]);
end
function phi = signed_distance_from_mask(mask, dx)
Dout = bwdist(mask) * dx; % distance outside region
Din  = bwdist(1-mask) * dx; % distance inside region
phi = Dout - Din; % negative inside
end
function [phi_x_plus,phi_x_minus,phi_y_plus,phi_y_minus] = compute_phi_derivatives(phi, dx, dy)
phi_x_plus  = ([phi(2:end,:); phi(end,:)] - phi)/dx;
phi_x_minus = (phi - [phi(1,:); phi(1:end-1,:)])/dx;
phi_y_plus  = ([phi(:,2:end), phi(:,end)] - phi)/dy;
phi_y_minus = (phi - [phi(:,1), phi(:,1:end-1)])/dy;
end
function H = H_Godunov(phi_x_minus, phi_x_plus, phi_y_minus, phi_y_plus, Vn)
H = zeros(size(phi_x_minus));
pos = Vn>0; neg = ~pos;
H(pos) = Vn(pos).*sqrt(max(phi_x_minus(pos),0).^2 + min(phi_x_plus(pos),0).^2 + ...
                       max(phi_y_minus(pos),0).^2 + min(phi_y_plus(pos),0).^2);
H(neg) = Vn(neg).*sqrt(max(phi_x_plus(neg),0).^2 + min(phi_x_minus(neg),0).^2 + ...
                       max(phi_y_plus(neg),0).^2 + min(phi_y_minus(neg),0).^2);
end
function phi_new = redistance(phi, dx, dy, nSteps)
phi0 = phi; dtau = 0.3*min(dx,dy);
for k=1:nSteps
    [phi_x_plus,phi_x_minus,phi_y_plus,phi_y_minus] = compute_phi_derivatives(phi, dx, dy);
    grad_pos = sqrt(max(phi_x_minus,0).^2 + min(phi_x_plus,0).^2 + ...
                    max(phi_y_minus,0).^2 + min(phi_y_plus,0).^2);
    grad_neg = sqrt(max(phi_x_plus,0).^2 + min(phi_x_minus,0).^2 + ...
                    max(phi_y_plus,0).^2 + min(phi_y_minus,0).^2);
    S = phi0./sqrt(phi0.^2 + dx^2);
    phi = phi - dtau*(S.*((S>=0).*(grad_pos-1) + (S<0).*(grad_neg-1)));
end
phi_new = phi;
end
function v_ext = gradient_extension(phi, v, dx, dy, nSteps)
dtau = 0.3*min(dx,dy);
for k=1:nSteps
    [phi_x_plus,phi_x_minus,phi_y_plus,phi_y_minus] = compute_phi_derivatives(phi, dx, dy);
    phi_x = 0.5*(phi_x_plus+phi_x_minus); phi_y = 0.5*(phi_y_plus+phi_y_minus);
    gradphi_norm = sqrt(phi_x.^2 + phi_y.^2 + 1e-10);
    vx_plus  = ([v(2:end,:); v(end,:)] - v)/dx;
    vx_minus = (v - [v(1,:); v(1:end-1,:)])/dx;
    vy_plus  = ([v(:,2:end), v(:,end)] - v)/dy;
    vy_minus = (v - [v(:,1), v(:,1:end-1)])/dy;
    vx = 0.5*(vx_plus + vx_minus); vy = 0.5*(vy_plus + vy_minus);
    S = phi./sqrt(phi.^2 + dx^2);
    v = v - dtau*(S.*((vx.*phi_x + vy.*phi_y)./gradphi_norm));
end
v_ext = v;
end
function c_sol = solve_c_field_vivo(rho, lambda, z, cB, L, Nx, Ny)
lambda_eff = lambda * (1 + z);
x = linspace(-L, L, Nx);
y = linspace(-L, L, Ny);
dx = x(2) - x(1);
dy = y(2) - y(1);
N = Nx * Ny;
idx_map = reshape(1:N, Ny, Nx); % map (row,column) to 1..N
% Assemble the same five-point operator with sparse triplets.
% Avoid repeated insertion into a sparse matrix inside spatial loops.
interior=idx_map(2:end-1,2:end-1);rows=interior(:);
up=idx_map(1:end-2,2:end-1);down=idx_map(3:end,2:end-1);
left=idx_map(2:end-1,1:end-2);right=idx_map(2:end-1,3:end);
r=rho(rows);inside=r>0;
diagonal=ones(N,1);diagonal(rows)=2/dx^2+2/dy^2+1;
diagonal(rows(inside))=2/dx^2+2/dy^2+lambda_eff*r(inside);
A=sparse([(1:N)';rows;rows;rows;rows], ...
 [(1:N)';up(:);down(:);left(:);right(:)], ...
 [diagonal;repmat(-1/dy^2,2*numel(rows),1);repmat(-1/dx^2,2*numel(rows),1)],N,N);
b=cB*ones(N,1);b(rows(inside))=0;
c_vec = A \ b;
c_sol = reshape(c_vec, [Ny, Nx]);
end

function B=gaussian_filter(A,sigma)
% Match imgaussfilt default support and replicate padding, no toolbox needed.
r=ceil(2*sigma); q=-r:r; h=exp(-q.^2/(2*sigma^2));h=h/sum(h);
rows=[ones(1,r),1:size(A,1),size(A,1)*ones(1,r)];
cols=[ones(1,r),1:size(A,2),size(A,2)*ones(1,r)];
B=conv2(h',h,A(rows,cols),'valid');
end
