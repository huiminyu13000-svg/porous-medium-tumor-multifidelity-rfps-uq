function results = run_experiment(mode, m, options)
%RUN_EXPERIMENT Compute BF/RFPS/TRF with a common physical growth coefficient.
% options may set T, n_candidate, n_pool, n_eval, seed, or make_plots.
% All errors use independent samples generated before basis construction.
if nargin < 1, mode = 'fresh'; end
if nargin < 2, m = 8; end
if nargin < 3, options = struct(); end
assert(ismember(m, [8 64]), 'm must be 8 or 64.');
assert(any(strcmp(mode, {'fresh', 'fresh_reference', 'cached', 'short', 'pde_check'})), ...
    'Unknown computation mode.');
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'solver'));
out = fullfile(root, 'results');
if ~exist(out, 'dir'), mkdir(out); end
source = load(fullfile(root, 'data', sprintf('bf_compute_m%02d.mat', m)));
settings = struct('schema', 2, 'm', double(m), 'T', double(source.T_final), ...
    'dx', double(source.dx), 'dt_pme', double(source.dt_PME), ...
    'dt_ls', double(source.dt_HS), 'coarse_dx', 0.1, 'coarse_dt', 0.002, ...
    'coarse_m', 8, 'growth_base', 0.5, 'n_candidate', 500, 'n_pool', 100, ...
    'n_eval', 20, 'seed', 20261123, 'initial_basis', 15, 'max_add', 3, ...
    'rfps_alpha', double(source.tri_alpha), 'rfps_weight', double(source.tri_lambda), ...
    'rfps_spatial_weights', source.tri_wvec, 'make_plots', true);
if strcmp(mode, 'short'), settings.T = 0.05; end
allowed = {'T', 'n_candidate', 'n_pool', 'n_eval', 'seed', 'make_plots'};
names = fieldnames(options);
for j = 1:numel(names)
    assert(ismember(names{j}, allowed), 'Unknown option: %s', names{j});
    settings.(names{j}) = options.(names{j});
end
assert(settings.n_candidate >= 18 && settings.n_candidate <= size(source.z, 2), ...
    'Use between 18 and 500 candidates.');
assert(settings.n_pool >= 18 && settings.n_pool <= size(source.z_test, 2), ...
    'Use between 18 and 100 enrichment-pool samples.');
assert(settings.n_eval >= 2 && settings.n_eval == fix(settings.n_eval), 'n_eval must be an integer >= 2.');
assert(all([settings.n_candidate settings.n_pool] == fix([settings.n_candidate settings.n_pool])), ...
    'Sample counts must be integers.');
for dt = [settings.dt_pme settings.dt_ls settings.coarse_dt]
    assert(settings.T > 0 && abs(round(settings.T/dt)*dt-settings.T) < 1e-12, ...
        'T must be a positive multiple of all time steps (use T=0.05 for a short run).');
end
candidate_z = source.z(:, 1:settings.n_candidate);
pool_z = source.z_test(:, 1:settings.n_pool);
previous_rng = rng;
rng(settings.seed, 'twister');
eval_z = 2 * rand(size(candidate_z, 1), settings.n_eval) - 1;
rng(previous_rng);
assert(~any(ismember(eval_z', [candidate_z pool_z]', 'rows')), ...
    'Evaluation parameters overlap with construction parameters.');
[X, Y] = meshgrid(-2.5:settings.dx:2.5);
[Xc, Yc] = meshgrid(-2.5:settings.coarse_dx:2.5);
cache_path = fullfile(out, sprintf('basis_cache_m%02d.mat', m));
timer = tic;
counts = struct('candidate_ls', 0, 'pool_ls', 0, 'coarse_pme', 0, ...
    'candidate_hf', 0, 'pool_hf', 0, 'evaluation_ls', 0, 'evaluation_hf', 0);
if strcmp(mode, 'pde_check')
    z = pool_z(:, 1);
    low = solve_sample(z, 'ls');
    high = solve_sample(z, 'hf');
    coarse = solve_sample(z, 'coarse');
    results = struct('settings', settings, 'z', z, 'L', low, 'H', high, 'C', coarse);
    assert(all(isfinite([low(:); high(:); coarse(:)])), 'Nonfinite PDE result.');
    save(fullfile(out, sprintf('pde_check_m%02d.mat', m)), '-struct', 'results', '-v7');
    return;
end
if strcmp(mode, 'cached')
    assert(exist(cache_path, 'file') == 2, 'Run fresh or short first to create the current-model cache.');
    stored = load(cache_path);
    expected = settings; expected.make_plots = true;
    actual = stored.cache.settings; actual.make_plots = true;
    assert(isequal(expected, actual) && isequal(stored.cache.candidate_z, candidate_z) && ...
        isequal(stored.cache.pool_z, pool_z) && isequal(stored.cache.eval_z, eval_z), ...
        'Cache settings differ. Recompute with fresh or use matching options.');
    cache = stored.cache;
else
    n_space = numel(X);
    cache = struct('settings', settings, 'candidate_z', candidate_z, ...
        'pool_z', pool_z, 'eval_z', eval_z);
    cache.candidate_low = zeros(n_space, settings.n_candidate);
    cache.pool_low = zeros(n_space, settings.n_pool);
    cache.pool_coarse = zeros(n_space, settings.n_pool);
    cache.candidate_high = nan(n_space, settings.n_candidate);
    cache.pool_high = nan(n_space, settings.n_pool);
    for j = 1:settings.n_candidate
        cache.candidate_low(:, j) = solve_sample(candidate_z(:, j), 'ls');
        counts.candidate_ls = counts.candidate_ls + 1;
        if mod(j, 25) == 0, fprintf('Candidate LS %d/%d\n', j, settings.n_candidate); end
    end
    for j = 1:settings.n_pool
        cache.pool_low(:, j) = solve_sample(pool_z(:, j), 'ls');
        cache.pool_coarse(:, j) = solve_sample(pool_z(:, j), 'coarse');
        counts.pool_ls = counts.pool_ls + 1;
        counts.coarse_pme = counts.coarse_pme + 1;
        if strcmp(mode, 'fresh_reference'), get_pool_hf(j); end
        if mod(j, 20) == 0, fprintf('Enrichment pool %d/%d\n', j, settings.n_pool); end
    end
end
% Baseline and RFPS share the same initial Cholesky basis.
[baseline_ids, ~] = select_basis(cache.candidate_low, candidate_z, ...
    settings.initial_basis + settings.max_add);
initial_count = min(settings.initial_basis, numel(baseline_ids));
initial_ids = baseline_ids(1:initial_count);
initial_low = cache.candidate_low(:, initial_ids);
initial_high = zeros(size(initial_low));
for j = 1:initial_count, initial_high(:, j) = get_candidate_hf(initial_ids(j)); end
[rfps_low, rfps_high, rfps_added] = enrich_basis(initial_low, initial_high, ...
    cache.pool_low, @get_pool_hf, settings.max_add, settings.rfps_alpha, true, ...
    settings.rfps_spatial_weights, settings.rfps_weight, []);
[trf_ids, ~] = select_basis(cache.pool_coarse, pool_z, settings.initial_basis);
trf_seed_high = zeros(numel(X), numel(trf_ids));
for j = 1:numel(trf_ids), trf_seed_high(:, j) = get_pool_hf(trf_ids(j)); end
[trf_low, trf_high, trf_added] = enrich_basis(cache.pool_low(:, trf_ids), ...
    trf_seed_high, cache.pool_low, @get_pool_hf, settings.max_add, ...
    settings.rfps_alpha, true, settings.rfps_spatial_weights, settings.rfps_weight, []);
baseline_low = cache.candidate_low(:, baseline_ids);
baseline_high = zeros(size(baseline_low));
for j = 1:numel(baseline_ids), baseline_high(:, j) = get_candidate_hf(baseline_ids(j)); end
construction_seconds = toc(timer);
% The independent evaluation set is first used after all selection is complete.
timer = tic;
if ~strcmp(mode, 'cached')
    cache.eval_low = zeros(numel(X), settings.n_eval);
    cache.eval_high = zeros(numel(X), settings.n_eval);
    for j = 1:settings.n_eval
        cache.eval_low(:, j) = solve_sample(eval_z(:, j), 'ls');
        cache.eval_high(:, j) = solve_sample(eval_z(:, j), 'hf');
        counts.evaluation_ls = counts.evaluation_ls + 1;
        counts.evaluation_hf = counts.evaluation_hf + 1;
        fprintf('Independent evaluation %d/%d\n', j, settings.n_eval);
    end
end
evaluation_seconds = toc(timer);
timer = tic;
results = struct('settings', settings, 'X', X, 'Y', Y, 'evaluation_parameters', eval_z, ...
    'baseline_indices', baseline_ids, 'rfps_initial_indices', initial_ids, ...
    'rfps_enrichment_indices', rfps_added, 'trf_initial_indices', trf_ids, ...
    'trf_enrichment_indices', trf_added, 'reference_mean', mean(cache.eval_high, 2), ...
    'reference_std', std(cache.eval_high, 1, 2));
results.bf = evaluate_basis(baseline_low, baseline_high, cache.eval_low, cache.eval_high);
results.rfps = evaluate_basis(rfps_low, rfps_high, cache.eval_low, cache.eval_high);
results.trf = evaluate_basis(trf_low, trf_high, cache.eval_low, cache.eval_high);
results.common_budget = min([results.bf.n_basis, results.rfps.n_basis, results.trf.n_basis]);
r = results.common_budget;
results.equal_budget_rmse = [results.bf.rmse(r), results.rfps.rmse(r), results.trf.rmse(r)];
results.equal_budget_relative_l2 = [results.bf.relative_l2(r), results.rfps.relative_l2(r), results.trf.relative_l2(r)];
results.basis = struct('bf_low', baseline_low, 'bf_high', baseline_high, ...
    'rfps_low', rfps_low, 'rfps_high', rfps_high, 'trf_low', trf_low, 'trf_high', trf_high);
results.new_solver_calls = counts;
results.construction_seconds = construction_seconds;
results.evaluation_seconds = evaluation_seconds;
results.projection_seconds = toc(timer);
assert(all(isfinite(results.equal_budget_rmse)), 'Nonfinite reconstruction error.');
% Atomic replacement prevents an interrupted write from becoming a valid cache.
temporary = [tempname(out) '.mat'];
save(temporary, 'cache', '-v7');
movefile(temporary, cache_path, 'f');
result_path = fullfile(out, sprintf('computed_%s_m%02d.mat', mode, m));
save(result_path, '-struct', 'results', '-v7');
if settings.make_plots, plot_computed(results, out); end
fprintf('m=%d, T=%g, independent N=%d, common HF basis count=%d\n', m, settings.T, settings.n_eval, r);
fprintf('Equal-budget RMSE (BF/RFPS/TRF): %.7g %.7g %.7g\n', results.equal_budget_rmse);

    function vector = get_candidate_hf(index)
        if any(isnan(cache.candidate_high(:, index)))
            cache.candidate_high(:, index) = solve_sample(candidate_z(:, index), 'hf');
            counts.candidate_hf = counts.candidate_hf + 1;
        end
        vector = cache.candidate_high(:, index);
    end

    function vector = get_pool_hf(index)
        if any(isnan(cache.pool_high(:, index)))
            cache.pool_high(:, index) = solve_sample(pool_z(:, index), 'hf');
            counts.pool_hf = counts.pool_hf + 1;
        end
        vector = cache.pool_high(:, index);
    end

    function vector = solve_sample(z, fidelity)
        factor = 1 + sum(z(:) ./ (2 * (1:numel(z))'));
        lambda = 50 * factor;
        cB = 20 * factor;
        growth_multiplier = 1 + z(1);
        radius = 0.45 * (1 + z(1));
        initial = 0.95 * double(X.^2 + Y.^2 < radius^2);
        switch fidelity
            case 'ls'
                field = solve_heleshaw(0, settings.T, settings.dx, settings.dx, ...
                    settings.dt_ls, lambda, cB, growth_multiplier, initial);
            case 'hf'
                field = solve_pme(0, settings.T, settings.m, settings.dx, settings.dx, ...
                    settings.dt_pme, lambda, cB, growth_multiplier, initial);
            case 'coarse'
                initial_coarse = 0.95 * double(Xc.^2 + Yc.^2 < radius^2);
                coarse = solve_cmpme(0, settings.T, settings.coarse_m, settings.coarse_dx, ...
                    settings.coarse_dt, lambda, cB, growth_multiplier, initial_coarse);
                field = interp2(Xc, Yc, coarse, X, Y, 'linear', 0);
        end
        assert(all(isfinite(field(:))), 'Nonfinite PDE solution.');
        vector = field(:);
    end
end
