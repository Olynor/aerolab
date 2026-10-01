function analyse_lift_curve()
% ANALYSE_LIFT_CURVE  Project AEROLAB
%   Reads ../08-Data/lift_runs.csv, computes lift-coefficient vs angle of
%   attack for each airfoil, fits the lift-curve slope, finds stall, and
%   (optionally) overlays XFOIL polars (../08-Data/polar_<airfoil>.txt).
%
%   Just run:  >> analyse_lift_curve
%
%   CSV columns expected:
%     timestamp, airfoil, airspeed_ms, alpha_deg, lift_g, lift_N, ...

% ---------------- CONFIG ----------------
CHORD   = 0.060;     % m
SPAN    = 0.078;     % m
RHO_AIR = 1.225;     % kg/m^3
MU      = 1.81e-5;   % Pa.s
S       = CHORD*SPAN;
LIN     = [-2 8];    % deg, linear region for slope fit
IDEAL_SLOPE_PER_DEG = 2*pi*pi/180;  % = 0.1097 /deg

here    = fileparts(mfilename('fullpath'));
csvPath = fullfile(here, '..', '08-Data', 'lift_runs.csv');
dataDir = fullfile(here, '..', '08-Data');

if ~isfile(csvPath)
    error('No data file at %s. Run the logger first (or use the Python --demo).', csvPath);
end

T = readtable(csvPath);

% lift in newtons (derive from grams if needed)
if ~ismember('lift_N', T.Properties.VariableNames)
    T.lift_N = T.lift_g/1000*9.81;
end
T.airspeed_ms = double(T.airspeed_ms);

% lift coefficient per row
q     = 0.5*RHO_AIR*T.airspeed_ms.^2;
T.CL  = T.lift_N ./ (q*S);

airspeed = mean(T.airspeed_ms, 'omitnan');
Re = RHO_AIR*airspeed*CHORD/MU;
fprintf('\n=== Lift-curve analysis ===\n');
fprintf('S = %.5f m^2 (chord %.0f mm x span %.0f mm)\n', S, CHORD*1000, SPAN*1000);
fprintf('Mean airspeed ~ %.2f m/s -> Re ~ %.0f\n\n', airspeed, Re);

% average repeats at each (airfoil, alpha)
Gmean = groupsummary(T, {'airfoil','alpha_deg'}, 'mean', 'CL');
Gstd  = groupsummary(T, {'airfoil','alpha_deg'}, 'std',  'CL');
Gmean.std_CL = Gstd.std_CL;
Gmean.std_CL(isnan(Gmean.std_CL)) = 0;

airfoils = unique(Gmean.airfoil, 'stable');
cmap = lines(numel(airfoils));

figure('Color','w'); hold on; grid on;
for k = 1:numel(airfoils)
    name = airfoils{k};
    rows = strcmp(Gmean.airfoil, name);
    a  = Gmean.alpha_deg(rows);
    cl = Gmean.mean_CL(rows);
    er = Gmean.std_CL(rows);
    [a, idx] = sort(a); cl = cl(idx); er = er(idx);

    errorbar(a, cl, er, 'o-', 'Color', cmap(k,:), 'DisplayName', [name ' (measured)']);

    % slope fit in linear region
    m = a >= LIN(1) & a <= LIN(2);
    if nnz(m) >= 2
        p = polyfit(a(m), cl(m), 1);
        slope = p(1); aL0 = -p(2)/p(1);
        xs = linspace(LIN(1), LIN(2), 50);
        plot(xs, polyval(p, xs), '--', 'Color', cmap(k,:), 'HandleVisibility','off');
        [clmax, im] = max(cl);
        fprintf('%s:\n', name);
        fprintf('   slope = %.4f /deg (%.3f /rad)\n', slope, slope*180/pi);
        fprintf('   ideal = %.4f /deg (2*pi /rad)\n', IDEAL_SLOPE_PER_DEG);
        fprintf('   alpha_L0 ~ %.2f deg\n', aL0);
        fprintf('   max CL ~ %.3f at alpha = %.1f deg (stall)\n\n', clmax, a(im));
    end

    % optional XFOIL overlay
    polar = fullfile(dataDir, ['polar_' name '.txt']);
    if isfile(polar)
        [xa, xcl] = read_xfoil_polar(polar);
        if ~isempty(xa)
            plot(xa, xcl, ':', 'Color', cmap(k,:), 'LineWidth', 2, ...
                'DisplayName', [name ' (XFOIL)']);
        end
    end
end

xs = linspace(LIN(1), LIN(2), 50);
plot(xs, IDEAL_SLOPE_PER_DEG*xs, 'k-', 'Color',[0 0 0 0.3], ...
    'DisplayName', 'thin-airfoil ideal');
yline(0,'Color',[.6 .6 .6]); xline(0,'Color',[.6 .6 .6]);
xlabel('Angle of attack  \alpha  (deg)');
ylabel('Lift coefficient  C_L');
title('AEROLAB: lift curves');
legend('Location','northwest','FontSize',8);

outPng = fullfile(dataDir, 'lift_curve_matlab.png');
saveas(gcf, outPng);
fprintf('Figure saved to %s\n', outPng);
end


function [alpha, cl] = read_xfoil_polar(path)
% Robustly read an XFOIL polar: find the header row containing alpha & CL,
% then read the numeric block beneath it.
alpha = []; cl = [];
lines = readlines(path);
hdr = 0; ia = 1; icl = 2;
for i = 1:numel(lines)
    low = lower(lines(i));
    if contains(low,'alpha') && contains(low,'cl')
        cols = lower(split(strtrim(lines(i))));
        ia  = find(cols=="alpha", 1);
        icl = find(cols=="cl", 1);
        if isempty(ia),  ia = 1;  end
        if isempty(icl), icl = 2; end
        hdr = i; break;
    end
end
if hdr == 0, return; end
for i = hdr+1:numel(lines)
    nums = str2double(split(strtrim(lines(i))));
    nums = nums(~isnan(nums));
    if numel(nums) >= max(ia,icl)
        alpha(end+1,1) = nums(ia); %#ok<AGROW>
        cl(end+1,1)    = nums(icl); %#ok<AGROW>
    end
end
[alpha, idx] = sort(alpha); cl = cl(idx);
end
