% 基础函数：保留数值更新公式；由 run_experiment 统一调用。
% Original author solver; executable formulas retained.
% Xu'an Dou attribution is retained where present in the original source.
function [phiL_new, phiH_new, idx_all] = enrich_basis( ...
    phiL, phiH, rhoL, rhoH, ...
    K_add, alpha, iterative, wvec, lambda, radius_list)
phiL_new = phiL;
phiH_new = phiH;
nSamples = size(rhoL,2);
added = false(1,nSamples);
idx_all = [];
selected_idx = [];
min_norm = 1e-8;
max_corr = 0.999;
min_gain = 1e-3;
dist_vec = inf(1,nSamples);
for i = 1:size(phiL,2)
    v = phiL(:,i);
    diff = vecnorm(rhoL - v);
    [~,idx] = min(diff);
    selected_idx = [selected_idx, idx];
end
selected_idx = unique(selected_idx);
idx_all = selected_idx;
for j = 1:nSamples
    if ~isempty(selected_idx)
        dmin = inf;
        for ii = 1:length(selected_idx)
            jj = selected_idx(ii);
            d = norm(rhoL(:,j) - rhoL(:,jj));
            if d < dmin
                dmin = d;
            end
        end
        dist_vec(j) = dmin;
    else
        dist_vec(j) = inf;
    end
end
for k = 1:K_add
    PL = phiL_new;
    % ---------- residual ----------
    A = PL'*PL;
    C = (A + alpha*eye(size(A))) \ (PL' * rhoL);
    rhoL_rec = PL * C;
    D = rhoL - rhoL_rec;
    if ~isempty(wvec)
        res = sqrt(sum((D.^2).*wvec,1) ./ sum(wvec));
    else
        res = sqrt(sum(D.^2,1)/size(D,1));
    end
    res_n = (res - min(res)) / max(max(res)-min(res),1e-12);
    % ---------- FPS ----------
    dist_n = dist_vec / max(max(dist_vec),1e-12);
    % ---------- score ----------
    score = lambda*res_n + (1-lambda)*dist_n;
    score(added) = -inf;
    [~, idx] = max(score);
    if ~isfinite(score(idx))
        break;
    end
    new_vec = rhoL(:,idx);
    proj = PL * ((A + alpha*eye(size(A))) \ (PL'*new_vec));
    ortho = new_vec - proj;
    if norm(ortho) < min_norm || ...
       max(abs(corr(new_vec,PL))) > max_corr || ...
       res(idx) < min_gain
        added(idx) = true;
        continue;
    end
    phiL_new = [phiL_new, new_vec];
    phiH_new = [phiH_new, rhoH(:,idx)];
    added(idx) = true;
    selected_idx = [selected_idx, idx];
    idx_all = [idx_all, idx];
    fprintf("  add %d: res=%.2e  FPS=%.3e\n", idx, res(idx), dist_vec(idx));
    for j = 1:nSamples
        if ~added(j)
            d = norm(rhoL(:,j) - rhoL(:,idx));
            if d < dist_vec(j)
                dist_vec(j) = d;
            end
        end
    end
    if ~iterative
        break;
    end
end
idx_all = unique(idx_all, 'stable');
fprintf("Selected indices: ");
fprintf("%d ", idx_all);
fprintf("\n");
end
