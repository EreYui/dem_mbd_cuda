% Draw residual particle fraction curves and the error between the present
% simulation and the digitized experimental curve.
clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

a = readNumericData(fullfile(scriptDir, "experiment.csv")); % Experiment
b = readNumericData(fullfile(scriptDir, "simulation.csv")); % Cited simulation
c = readNumericData(fullfile(scriptDir, "num_re.csv"));     % Present simulation

a = cleanCurve(a);
b = cleanCurve(b);
c = cleanCurve(c);

a = normalizeResidualFraction(a);
b = normalizeResidualFraction(b);
c = normalizeResidualFraction(c);

% Interpolate the experimental curve to the time points of the present
% simulation. The comparison is restricted to the overlapping time interval
% to avoid extrapolation.
tMin = max(min(c(:,1)), min(a(:,1)));
tMax = min(max(c(:,1)), max(a(:,1)));
if tMax < min(max(a(:,1)), max(b(:,1)))
    warning(['The present simulation data only overlap with the experimental ', ...
        'curve from %.4g s to %.4g s. The error curve is evaluated only ', ...
        'over this interval.'], tMin, tMax);
end
idx = c(:,1) >= tMin & c(:,1) <= tMax;
tErr = c(idx,1);
expInterp = interp1(a(:,1), a(:,2), tErr, 'linear');
yOurs = c(idx,2);
valid = isfinite(tErr) & isfinite(expInterp) & isfinite(yOurs);
tErr = tErr(valid);
expInterp = expInterp(valid);
yOurs = yOurs(valid);
err = yOurs - expInterp;

if isempty(err)
    error('No valid overlapping data points are available for error calculation.');
end

rmse = sqrt(mean(err.^2));
mae = mean(abs(err));
maxAbsErr = max(abs(err));
fprintf('Residual fraction error statistics (ours - experiment):\n');
fprintf('  RMSE = %.6f\n', rmse);
fprintf('  MAE  = %.6f\n', mae);
fprintf('  Max absolute error = %.6f\n', maxAbsErr);

fig = figure('Position', [100, 100, 560, 390], 'Color', 'w');
hold on; box on;

yyaxis left
p1 = plot(c(:,1), c(:,2), 'k-', 'LineWidth', 1.6); hold on;
p2 = plot(a(:,1), a(:,2), 'b-.', 'LineWidth', 1.4);
p3 = plot(b(:,1), b(:,2), 'r--', 'LineWidth', 1.4);
ylabel('Residual particle fraction');
ylim([0, 1.05]);
ax = gca;
ax.YAxis(1).Color = 'k';

yyaxis right
p4 = plot(tErr, err, '--', 'Color', [0.00, 0.45, 0.00], 'LineWidth', 1.5);
plot([0, 8], [0, 0], ':', 'Color', [0.35, 0.35, 0.35], ...
    'LineWidth', 0.8, 'HandleVisibility', 'off');
ylabel('Error (ours - experiment)');

errMax = max(abs(err));
if errMax < 1.0e-6
    ylim([-0.01, 0.01]);
else
    ylim(1.15 * [-errMax, errMax]);
end
ax = gca;
ax.YAxis(2).Color = 'k';

yyaxis left
xlim([0, 8]);
xticks([0, 2, 4, 6, 8]);
set(gca, 'XGrid', 'off', 'YGrid', 'off');
xlabel('Time (s)');

lgd = legend([p1, p2, p3, p4], ...
    {'Sim-ours', 'Experiment', 'Sim-cited', 'Error'}, ...
    'Location', 'northeast');
set(lgd, 'Box', 'off');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'LineWidth', 1.0);

outPath = fullfile(scriptDir, "..", "figs", "hopper_discharge_fraction_error.png");
if exist('exportgraphics', 'file')
    exportgraphics(fig, outPath, 'Resolution', 600);
else
    print(fig, outPath, '-dpng', '-r600');
end
fprintf('Figure saved to: %s\n', outPath);

function data = readNumericData(filePath)
    if isZipBasedSpreadsheet(filePath)
        tmpPath = [tempname, '.xlsx'];
        copyfile(char(filePath), tmpPath);
        cleanupObj = onCleanup(@() deleteIfExists(tmpPath));
        data = readmatrix(tmpPath);
    else
        data = readmatrix(filePath);
        if isempty(data) || all(isnan(data), 'all')
            data = load(filePath);
        end
    end

    data = data(:, 1:2);
    data = data(all(~isnan(data), 2), :);
end

function data = cleanCurve(data)
    data = data(:, 1:2);
    data = data(all(isfinite(data), 2), :);
    data = sortrows(data, 1);
    [~, ia] = unique(data(:,1), 'stable');
    data = data(ia, :);
end

function data = normalizeResidualFraction(data)
    if max(data(:,2)) > 1.5
        data(:,2) = data(:,2) ./ data(1,2);
    end
end

function deleteIfExists(filePath)
    if exist(filePath, 'file')
        delete(filePath);
    end
end

function tf = isZipBasedSpreadsheet(filePath)
    fid = fopen(filePath, 'r');
    if fid < 0
        error('Cannot open data file: %s', filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid));
    bytes = fread(fid, 2, 'uint8');
    tf = numel(bytes) == 2 && bytes(1) == 80 && bytes(2) == 75;
end
