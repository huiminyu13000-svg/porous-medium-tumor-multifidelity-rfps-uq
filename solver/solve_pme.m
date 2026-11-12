% 基础函数：保留数值更新公式；由 run_experiment 统一调用。
% Adapted from solver_PME_vivo_modified.m in the supplied archive.
% Numerical model parameters, including the original PME growth factor 0.5, are preserved.
function [rho_final, p_final] = solve_pme(k, T_final, m, dx, dy, dt, ...
                                             lambda_k, cB_k, G0_k, rho0)
% Convert integer MAT metadata before dividing time by the step size.
T_final=double(T_final);m=double(m);
a = -2.5; b = 2.5;
[X, Y] = meshgrid(a:dx:b, a:dx:b);
Nx=length(a:dx:b);
Ny=length(a:dx:b);
Nt = round(T_final / dt);
assert(Nt>=1 && abs(Nt*dt-T_final)<1e-12, 'Positive T must be a multiple of dt.');
lambda = lambda_k;
cB = cB_k;
G0 = G0_k;
z = 0;
if nargin < 10
    error('Please provide initial rho0 as input.');
end
rho = rho0;               % initial density
p   = rho.^(m-1) * m/(m-1);  % initial pressure
Velocity_x = @(X,Y) zeros(size(X));
Velocity_y = @(X,Y) zeros(size(Y));
u_e_full = Velocity_x(X,Y)';  % transpose to match rho orientation
v_e_full = Velocity_y(X,Y)';
nx = Nx-2; ny = Ny-2;
u_e = u_e_full(2:end-1,2:end-1);
v_e = v_e_full(2:end-1,2:end-1);
U_m = zeros(nx,ny); V_m = U_m;
uv_pp = U_m; uv_pm = U_m; uv_mm = U_m; uv_mp = U_m;
vu_pp = U_m; vu_pm = U_m; vu_mm = U_m; vu_mp = U_m;
V_up = U_m; V_low = U_m; U_up = U_m; U_low = U_m;
b_sol = zeros(nx*ny*2,1);
u = -(p(3:end,2:end-1) - p(1:end-2,2:end-1)) / (2*dx);
v = -(p(2:end-1,3:end) - p(2:end-1,1:end-2)) / (2*dy);
for step = 1:Nt
L = (b - a)/2;
[Nx_rho, Ny_rho] = size(rho);
c = solve_c_field_vivo(rho, lambda, z, cB, L, Nx_rho, Ny_rho);
 G_field = 0.5 * G0 * c;       % growth field
    rhoGGG = G_field;
    Build_Matrix();
    [sol,flag] = gmres(@fA, b_sol, min(400,numel(b_sol)), 1e-6, 100);
    if flag~=0, error('tumor:GMRES','PME velocity solve failed: %d',flag); end
    u_star = reshape(sol(1:nx*ny), nx, ny);
    v_star = reshape(sol(nx*ny+1:end), nx, ny);
    assert(all(1-dt*rhoGGG(:)>0),'Reaction denominator is nonpositive. Reduce dt.');
    Convec_Cor();
    assert(all(isfinite(rho(:))), 'Nonfinite PME density.');
end
rho_final = rho;
p_final   = p;
    function Build_Matrix
        r_ip = ((rho(3:end,2:end-1)+rho(2:end-1,2:end-1))/2).^(m-2);
        r_im = ((rho(1:end-2,2:end-1)+rho(2:end-1,2:end-1))/2).^(m-2);
        r_jp = ((rho(2:end-1,3:end)+rho(2:end-1,2:end-1))/2).^(m-2);
        r_jm = ((rho(2:end-1,1:end-2)+rho(2:end-1,2:end-1))/2).^(m-2);
        U_m = 1 + m*(dt/dx^2) * rho(2:end-1,2:end-1) .* (r_ip + r_im);
        U_up = -m*(dt/dx^2) * r_ip .* rho(3:end,2:end-1); U_up = U_up(1:end-1,:);
        U_low = -m*(dt/dx^2) * r_im .* rho(1:end-2,2:end-1); U_low = U_low(2:end,:);
        V_m = 1 + m*(dt/dy^2) * rho(2:end-1,2:end-1) .* (r_jp + r_jm);
        V_up = -m*(dt/dy^2) * r_jp .* rho(2:end-1,3:end); V_up = V_up(:,1:end-1);
        V_low = -m*(dt/dy^2) * r_jm .* rho(2:end-1,1:end-2); V_low = V_low(:,2:end);
        lam = m*dt/dx/dy/4;
        uv_pp = -lam * (rho(3:end,2:end-1).^(m-2)) .* rho(3:end,3:end); uv_pp = uv_pp(1:end-1,1:end-1);
        uv_pm =  lam * (rho(3:end,2:end-1).^(m-2)) .* rho(3:end,1:end-2); uv_pm = uv_pm(1:end-1,2:end);
        uv_mp =  lam * (rho(1:end-2,2:end-1).^(m-2)) .* rho(1:end-2,3:end); uv_mm = -lam*(rho(1:end-2,2:end-1).^(m-2)).*rho(1:end-2,1:end-2);
        uv_mp = uv_mp(2:end,1:end-1); uv_mm = uv_mm(2:end,2:end);
        vu_pp = -lam * (rho(2:end-1,3:end).^(m-2)) .* rho(3:end,3:end);
        vu_pm =  lam * (rho(2:end-1,3:end).^(m-2)) .* rho(3:end,1:end-2);
        vu_pp = vu_pp(1:end-1,1:end-1); vu_pm = vu_pm(1:end-1,2:end);
        vu_mp =  lam * (rho(2:end-1,1:end-2).^(m-2)) .* rho(1:end-2,3:end);
        vu_mm = -lam * (rho(2:end-1,1:end-2).^(m-2)) .* rho(1:end-2,1:end-2);
        vu_mp = vu_mp(2:end,1:end-1); vu_mm = vu_mm(2:end,2:end);
        p = rho.^(m-1) * m/(m-1);
        rhoG = rho.^(m-1) .* rhoGGG;
        bu = u - m*(dt/(2*dx))*(rhoG(3:end,2:end-1) - rhoG(1:end-2,2:end-1));
        bv = v - m*(dt/(2*dy))*(rhoG(2:end-1,3:end) - rhoG(2:end-1,1:end-2));
        bu = bu - (U_m-1).*u_e;
        bu(1:end-1,:) = bu(1:end-1,:) - U_up.*u_e(2:end,:);
        bu(2:end,:)   = bu(2:end,:)   - U_low.*u_e(1:end-1,:);
        bu(1:end-1,1:end-1) = bu(1:end-1,1:end-1) - uv_pp.*v_e(2:end,2:end);
        bu(1:end-1,2:end)   = bu(1:end-1,2:end)   - uv_pm.*v_e(2:end,1:end-1);
        bu(2:end,1:end-1)   = bu(2:end,1:end-1)   - uv_mp.*v_e(1:end-1,2:end);
        bu(2:end,2:end)     = bu(2:end,2:end)     - uv_mm.*v_e(1:end-1,1:end-1);
        bv = bv - (V_m-1).*v_e;
        bv(:,1:end-1) = bv(:,1:end-1) - V_up.*v_e(:,2:end);
        bv(:,2:end)   = bv(:,2:end)   - V_low.*v_e(:,1:end-1);
        bv(1:end-1,1:end-1) = bv(1:end-1,1:end-1) - vu_pp.*u_e(2:end,2:end);
        bv(1:end-1,2:end)   = bv(1:end-1,2:end)   - vu_pm.*u_e(2:end,1:end-1);
        bv(2:end,1:end-1)   = bv(2:end,1:end-1)   - vu_mp.*u_e(1:end-1,2:end);
        bv(2:end,2:end)     = bv(2:end,2:end)     - vu_mm.*u_e(1:end-1,1:end-1);
        b_sol = [reshape(bu,nx*ny,1); reshape(bv,nx*ny,1)];
    end
    function ax = fA(x_tmp)
        u_tmp = reshape(x_tmp(1:nx*ny), nx, ny);
        v_tmp = reshape(x_tmp(nx*ny+1:end), nx, ny);
        ax_u = U_m .* u_tmp;
        ax_u(1:end-1,:) = ax_u(1:end-1,:) + U_up.*u_tmp(2:end,:);
        ax_u(2:end,:)   = ax_u(2:end,:)   + U_low.*u_tmp(1:end-1,:);
        ax_u(1:end-1,1:end-1) = ax_u(1:end-1,1:end-1) + uv_pp.*v_tmp(2:end,2:end);
        ax_u(1:end-1,2:end)   = ax_u(1:end-1,2:end)   + uv_pm.*v_tmp(2:end,1:end-1);
        ax_u(2:end,1:end-1)   = ax_u(2:end,1:end-1)   + uv_mp.*v_tmp(1:end-1,2:end);
        ax_u(2:end,2:end)     = ax_u(2:end,2:end)     + uv_mm.*v_tmp(1:end-1,1:end-1);
        ax_v = V_m .* v_tmp;
        ax_v(:,1:end-1) = ax_v(:,1:end-1) + V_up.*v_tmp(:,2:end);
        ax_v(:,2:end)   = ax_v(:,2:end)   + V_low.*v_tmp(:,1:end-1);
        ax_v(1:end-1,1:end-1) = ax_v(1:end-1,1:end-1) + vu_pp.*u_tmp(2:end,2:end);
        ax_v(1:end-1,2:end)   = ax_v(1:end-1,2:end)   + vu_pm.*u_tmp(2:end,1:end-1);
        ax_v(2:end,1:end-1)   = ax_v(2:end,1:end-1)   + vu_mp.*u_tmp(1:end-1,2:end);
        ax_v(2:end,2:end)     = ax_v(2:end,2:end)     + vu_mm.*u_tmp(1:end-1,1:end-1);
        ax = [reshape(ax_u,nx*ny,1); reshape(ax_v,nx*ny,1)];
    end
function rst=getdFx(rho,u_star) % difference in x, fixed y
rho=rho(:,2:end-1);% 0,1,...Nx
diff_l=(rho(2:end-1,:)-rho(1:end-2,:));%rho_{i}-rho_{i-1} i=1,...Nx-1
diff_r=(rho(3:end,:)-rho(2:end-1,:));%rho_{i+1}-rho_{i} i=1,...Nx-1
   parial_rho=min(diff_l,diff_r).*(diff_l>0).*(diff_r>0)+max(diff_l,diff_r).*(diff_l<0).*(diff_r<0);
   rho_l=rho(2:end-2,:)+1/2*parial_rho(1:end-1,:);
   rho_r=rho(3:end-1,:)-1/2*parial_rho(2:end,:);
   F=1/2*(rho_l.*u_star+rho_r.*u_star-abs(u_star).*(rho_r-rho_l));
   F(end+1,:)=0;F(end+1,:)=0; F(2:end,:)=F(1:end-1,:); F(1,:)=0;
   rst=F(2:end,:)-F(1:end-1,:); 
end
function rst=getdFy(rho,v_star) % difference in y, fixed x
rho=rho(2:end-1,:);% 0,1,...Ny
diff_l=(rho(:,2:end-1)-rho(:,1:end-2));%rho_{i}-rho_{i-1} i=1,...Ny-1
diff_r=(rho(:,3:end)-rho(:,2:end-1));%rho_{i+1}-rho_{i} i=1,...Ny-1
   parial_rho=min(diff_l,diff_r).*(diff_l>0).*(diff_r>0)+max(diff_l,diff_r).*(diff_l<0).*(diff_r<0);
   rho_l=rho(:,2:end-2)+1/2*parial_rho(:,1:end-1);
   rho_r=rho(:,3:end-1)-1/2*parial_rho(:,2:end);
   F=1/2*(rho_l.*v_star+rho_r.*v_star-abs(v_star).*(rho_r-rho_l));
   F(:,end+1)=0;F(:,end+1)=0; F(:,2:end)=F(:,1:end-1); F(:,1)=0;
   rst=F(:,2:end)-F(:,1:end-1); 
end
    function Convec_Cor
        u_sg = (u_star(2:end,:)+u_star(1:end-1,:) + u_e(2:end,:) + u_e(1:end-1,:))/2;
        v_sg = (v_star(:,2:end)+v_star(:,1:end-1) + v_e(:,2:end)+v_e(:,1:end-1))/2;
        diff_Fx = getdFx(rho, u_sg);
        diff_Fy = getdFy(rho, v_sg);
        F_total = dt/dx*diff_Fx + dt/dy*diff_Fy;
        rho(2:end-1,2:end-1) = (rho(2:end-1,2:end-1) - F_total) ./ (1 - dt*rhoGGG(2:end-1,2:end-1));
        p = rho.^(m-1) * m/(m-1);
        u = -(p(3:end,2:end-1)-p(1:end-2,2:end-1))/(2*dx);
        v = -(p(2:end-1,3:end)-p(2:end-1,1:end-2))/(2*dy);
    end
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
