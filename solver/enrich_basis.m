function [phiL_new, phiH_new, added_indices] = enrich_basis( ...
    phiL, phiH, rhoL, rhoH, K_add, alpha, iterative, wvec, lambda, radius_list)
%ENRICH_BASIS Enrich a paired basis using LF residuals and LF distances.
% rhoH is either the matching HF matrix or a function handle for a pool index.
% added_indices contains only accepted enrichment indices in the supplied pool.
% radius_list is retained for compatibility with the original call signature.
assert(size(phiL, 2) == size(phiH, 2), 'Paired bases must have equal column counts.');
assert(lambda >= 0 && lambda <= 1 && alpha >= 0, 'Invalid RFPS parameters.');
phiL_new = phiL;
phiH_new = phiH;
n_samples = size(rhoL, 2);
added_indices = [];
unavailable = false(1, n_samples);
distances = inf(1, n_samples);
for j = 1:size(phiL, 2)
    distances = min(distances, vecnorm(rhoL - phiL(:, j), 2, 1));
end
% Exclude pool columns already represented by an identical LF basis column.
scale = max(1, max(vecnorm(rhoL, 2, 1)));
unavailable(distances <= 1e-12 * scale) = true;
if ~isempty(wvec)
    wvec = wvec(:);
    assert(numel(wvec) == size(rhoL, 1) && all(wvec >= 0) && sum(wvec) > 0, ...
        'Invalid spatial RFPS weights.');
end
while numel(added_indices) < K_add && any(~unavailable)
    gram = phiL_new' * phiL_new;
    regularized = gram + alpha * eye(size(gram));
    coefficients = regularized \ (phiL_new' * rhoL);
    residual = rhoL - phiL_new * coefficients;
    if isempty(wvec)
        residual_norm = sqrt(mean(residual.^2, 1));
    else
        residual_norm = sqrt(sum(residual.^2 .* wvec, 1) / sum(wvec));
    end
    residual_score = (residual_norm - min(residual_norm)) / ...
        max(max(residual_norm) - min(residual_norm), 1e-12);
    distance_score = distances / max(max(distances), 1e-12);
    score = lambda * residual_score + (1 - lambda) * distance_score;
    score(unavailable) = -inf;
    [~, index] = max(score);
    unavailable(index) = true;
    new_vector = rhoL(:, index);
    orthogonal = new_vector - phiL_new * (regularized \ (phiL_new' * new_vector));
    % Pearson correlation, implemented without the Statistics Toolbox.
    centered = new_vector - mean(new_vector);
    centered_basis = phiL_new - mean(phiL_new, 1);
    denominator = norm(centered) * vecnorm(centered_basis, 2, 1);
    correlations = zeros(size(denominator));
    valid = denominator > 0;
    correlations(valid) = abs(centered' * centered_basis(:, valid)) ./ denominator(valid);
    if norm(orthogonal) < 1e-8 || max(correlations) > 0.999 || residual_norm(index) < 1e-3
        continue;
    end
    if isa(rhoH, 'function_handle')
        high_vector = rhoH(index);
    else
        high_vector = rhoH(:, index);
    end
    assert(numel(high_vector) == size(phiH, 1) && all(isfinite(high_vector(:))), ...
        'Invalid HF enrichment vector.');
    phiL_new(:, end + 1) = new_vector;
    phiH_new(:, end + 1) = high_vector(:);
    added_indices(end + 1) = index;
    distances = min(distances, vecnorm(rhoL - new_vector, 2, 1));
    fprintf('Accepted RFPS pool index %d (%d/%d)\n', index, numel(added_indices), K_add);
    if ~iterative
        break;
    end
end
end
