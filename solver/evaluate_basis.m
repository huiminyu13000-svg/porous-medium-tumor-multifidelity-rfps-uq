function result = evaluate_basis(low_basis, high_basis, low_eval, high_eval)
%EVALUATE_BASIS Evaluate a fixed paired basis on independent parameter samples.
assert(size(low_basis, 2) == size(high_basis, 2), 'Paired basis mismatch.');
assert(isequal(size(low_eval), size(high_eval)), 'Evaluation matrix mismatch.');
assert(all(isfinite([low_basis(:); high_basis(:); low_eval(:); high_eval(:)])), ...
    'Nonfinite basis or evaluation data.');
count = size(low_basis, 2);
result.rmse = zeros(count, 1);
result.relative_l2 = nan(count, 1);
result.mean_rmse = zeros(count, 1);
result.std_rmse = zeros(count, 1);
reference_norm = norm(high_eval, 'fro');
for r = 1:count
    low = low_basis(:, 1:r);
    gram = low' * low;
    if r > 1 && cond(gram) > 1e12
        gram = gram + 1e-8 * trace(gram) / r * eye(r);
    end
    reconstruction = high_basis(:, 1:r) * (gram \ (low' * low_eval));
    difference = reconstruction - high_eval;
    result.rmse(r) = sqrt(mean(difference(:).^2));
    if reference_norm > 0
        result.relative_l2(r) = norm(difference, 'fro') / reference_norm;
    end
    result.mean_rmse(r) = sqrt(mean(mean(difference, 2).^2));
    result.std_rmse(r) = sqrt(mean((std(reconstruction, 1, 2) - std(high_eval, 1, 2)).^2));
end
result.reconstruction = reconstruction;
result.mean = mean(reconstruction, 2);
result.std = std(reconstruction, 1, 2);
result.n_basis = count;
result.n_eval = size(high_eval, 2);
end
