function plot_computed(result, out)
%PLOT_COMPUTED Plot independent errors and computed density statistics.
if nargin < 2
    out = fullfile(fileparts(mfilename('fullpath')), 'results');
end
if ~exist(out, 'dir'), mkdir(out); end
names = {'bf', 'rfps', 'trf'};
labels = {'BF', 'BF + RFPS', 'TRF'};
styles = {'ko-', 'bs--', 'rd-.'};
f = figure('visible', 'off', 'Color', 'w', 'position', [100 100 1100 420]);
for panel = 1:2
    subplot(1, 2, panel); hold on;
    for j = 1:3
        values = result.(names{j});
        if panel == 1, errors = values.rmse; else, errors = values.relative_l2; end
        semilogy(1:numel(errors), errors, styles{j}, 'LineWidth', 1.5);
    end
    set(gca, 'YScale', 'log');
    xlabel('HF basis count');
    if panel == 1, ylabel('Independent spatial RMSE'); else, ylabel('Independent relative L2 error'); end
    legend(labels, 'Location', 'best'); grid on;
end
print(f, fullfile(out, sprintf('computed_errors_m%02d.png', result.settings.m)), '-dpng', '-r150');
close(f);
for statistic = {'mean', 'std'}
    stat = statistic{1};
    f = figure('visible', 'off', 'Color', 'w', 'position', [100 100 1400 650]);
    reference = reshape(result.(['reference_' stat]), size(result.X));
    fields = cell(1, 4);
    for j = 1:3, fields{j} = reshape(result.(names{j}).(stat), size(result.X)); end
    fields{4} = reference;
    color_limits = [min(cellfun(@(v) min(v(:)), fields)), max(cellfun(@(v) max(v(:)), fields))];
    if color_limits(1) == color_limits(2), color_limits(2) = color_limits(1) + 1; end
    for j = 1:4
        subplot(2, 4, j);
        imagesc(result.X(1, :), result.Y(:, 1), fields{j});
        axis image; set(gca, 'YDir', 'normal'); caxis(color_limits); colorbar;
        if j <= 3
            title(sprintf('%s %s (%d HF)', labels{j}, stat, result.(names{j}).n_basis));
        else
            title(['HF reference ' stat]);
        end
        xlabel('x'); ylabel('y');
        subplot(2, 4, j + 4);
        row = round(size(result.X, 1) / 2);
        plot(result.X(1, :), fields{j}(row, :), 'b-', result.X(1, :), reference(row, :), 'k--', 'LineWidth', 1.2);
        xlabel('x'); ylabel(stat); grid on;
    end
    print(f, fullfile(out, sprintf('computed_%s_m%02d.png', stat, result.settings.m)), '-dpng', '-r150');
    close(f);
end
end
