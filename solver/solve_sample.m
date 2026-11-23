function [vector,info] = solve_sample(z, fidelity, cfg)
%SOLVE_SAMPLE Use one parameter map and initial condition at every fidelity.
p=model_parameters(z,cfg); grid=cfg.domain(1):cfg.dx:cfg.domain(2);
[X,Y]=meshgrid(grid); rho0=cfg.initial_density*double(sqrt(X.^2+Y.^2)<p.radius);
switch fidelity
    case 'hf'
        [field,~,info]=solve_pme(0,cfg.T,cfg.m,cfg.dx,cfg.dy,cfg.dt_pme, ...
            p.lambda,p.cB,p.growth_multiplier,rho0);
    case 'ls'
        [field,~,info]=solve_heleshaw(0,cfg.T,cfg.dx,cfg.dy,cfg.dt_ls, ...
            p.lambda,p.cB,p.growth_multiplier,rho0);
    case 'coarse'
        gc=cfg.domain(1):cfg.coarse_dx:cfg.domain(2);[Xc,Yc]=meshgrid(gc);
        initial=cfg.initial_density*double(sqrt(Xc.^2+Yc.^2)<p.radius);
        [coarse,~,info]=solve_cmpme(0,cfg.T,cfg.coarse_m,cfg.coarse_dx,cfg.coarse_dt, ...
            p.lambda,p.cB,p.growth_multiplier,initial);
        field=interp2(Xc,Yc,coarse,X,Y,'linear',0);
    otherwise,error('Unknown fidelity.');
end
assert(all(isfinite(field(:))),'Nonfinite PDE solution.');vector=field(:);
end
