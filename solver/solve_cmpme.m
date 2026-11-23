function [rho,p,info] = solve_cmpme(k,T,m,dx,dt,lambda,cB,G0_k,rho0)
%SOLVE_CMPME The same PME discretization on the coarse grid.
% G0_k is the multiplier; solve_pme applies the factor 0.5 exactly once.
[rho,p,info]=solve_pme(k,T,m,dx,dx,dt,lambda,cB,G0_k,rho0);
end
