function result = run_reference(m, options)
%RUN_REFERENCE Compute HF and reconstructed statistics on 1200 fresh samples.
% Default N_stats=1200 is set in experiment_config.m. No input MAT is needed.
% A matching basis is reused if present; otherwise it is generated first.
% Matching progress is saved every ten samples and resumed automatically.
if nargin<1,m=8;end
if nargin<2,options=struct();end
result=run_experiment('reference',m,options);
end
