function cfg = experiment_config(m, options)
%EXPERIMENT_CONFIG Circular nutrient experiment, PPT slides 14-15.
if nargin < 1, m = 8; end
if nargin < 2, options = struct(); end
assert(ismember(m,[8 64]), 'm must be 8 or 64.');
cfg = struct('version','ppt-20261123','case_name','nutrient_circle', ...
    'm',double(m),'domain',[-2.5 2.5],'dx',0.05,'dy',0.05,'T',1, ...
    'dt_pme',0.0015,'dt_ls',0.03,'coarse_dx',0.1,'coarse_dt',0.006, ...
    'coarse_m',double(m),'parameter_floor',1e-6,'growth_base',0.5, ...
    'initial_density',0.95,'initial_radius',0.45, ...
    'n_train',500,'n_pool',100,'n_eval',500,'n_stats',1200, ...
    'initial_basis',15,'max_add',3,'rfps_weight',0.4,'rfps_alpha',1e-8,'rfps_spatial_sigma',0.15, ...
    'seed',20261123,'make_plots',true,'resume',true);
if m==64, cfg.dt_pme=0.0005;cfg.coarse_dt=0.002;end
allowed = {'T','n_train','n_pool','n_eval','n_stats','seed','make_plots','resume'};
names = fieldnames(options);
for j=1:numel(names)
    assert(ismember(names{j},allowed),'Unknown option: %s',names{j});
    cfg.(names{j})=options.(names{j});
end
assert(isscalar(cfg.T) && isfinite(cfg.T) && cfg.T>0,'T must be positive.');
for name={'n_train','n_pool','n_eval','n_stats'}
    v=cfg.(name{1}); assert(isscalar(v)&&isfinite(v)&&v==fix(v)&&v>=2,'Invalid sample count.');
end
assert(cfg.n_train>=18 && cfg.n_pool>=18,'Use at least 18 training and pool parameters.');
end
