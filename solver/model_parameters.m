function p = model_parameters(z, cfg)
%MODEL_PARAMETERS Physical parameters for five independent U[-1,1] inputs.
z=z(:);assert(numel(z)==5 && all(isfinite(z)) && all(abs(z)<=1),'Expected five inputs in [-1,1].');
s=max(1+sum(z./(2*(1:5)')),cfg.parameter_floor);
a1=max(1+z(1),cfg.parameter_floor);
p=struct('lambda',50*s,'cB',20*s,'growth_multiplier',a1, ...
    'G0',cfg.growth_base*a1,'radius',cfg.initial_radius*a1);
end
