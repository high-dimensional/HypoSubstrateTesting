%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      %%%     %%%    %%%%%%%   %%%     %
%                                      %%%     %%%   %%%        %%%     %
%                                      %%%     %%%   %%%        %%%     %
%                                      %%%     %%%   %%%        %%%     %
%                                       %%%%%%%%%     %%%%%%%   %%%%%%%%%
%                                                                       %
%                                                                       %
%                      HypoSubstrateTesting                             %
% Testing hypothesized neural substrates with lesions                   %
%                                                                       %
% This script contains the function for plotting the four-panel figure  %
% of one simulated scenario: example studies at the lowest and highest  %
% noise level with linear fits of the H and ~H compartmental volume;    %
%  signed -ln(p) of the Overlap:Compartment interaction against noise.  %
% -ln(p) of all iterations against noise (blue circles: unadjusted,     %
% green triangles: lesion volume-adjusted); split violins of the signed %
% -ln(p) per noise bin with the percentage of significant iterations    %
% that are H-implicated (top) and ~H-implicated (bottom).               %
% H: hypothesized compartment, ~H: outside it. Positive signed -ln(p)   %
% implicates H. Titles and labels are placed explicitly, with text      %
% widths from a Helvetica glyph table, because MATLAB cannot measure    %
% text without a display.                                               %
%                                                                       %
% input:                                                                %
% -matFile: result file <outDir>/hyp-<hypothesized>_true-<generator>.mat%
% output: <matFile>_figure.pdf                                          %
%                                                                       %
%                                                                       %
% Author: Patrik Bey, patrik.bey@ucl.ac.uk                              %
%                                                                       %
%                                                                       %
% last update: 2026/11/10.                                              %
%                                                                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_rldr_figure(matFile)
    set(groot, 'defaultAxesFontName', 'Helvetica', 'defaultTextFontName', 'Helvetica', 'defaultLegendFontName', 'Helvetica');
    matFile = char(matFile);
    if ~endsWith(matFile, '.mat'), error('Expected a .mat result file, got %s', matFile); end
    R = load(matFile);
    hyp = R.params.hypothesized;
    thresh = -log(0.05);
    blue = [0.35 0.65 0.90];
    green = [0.45 0.70 0.25];
    sig = char(963);

    fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 26 22]);
    theme(fig, 'light');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'loose', 'Padding', 'loose');

    ax = gobjects(1, 4);
    info = cell(1, 4);     % per axes: title, xlabel, ylabel, rotated x tick labels (optional)
    for k = 1:numel(R.examples)
        ax(k) = nexttile(tl);
        exampleScatter(R.examples(k), hyp);
        info{k} = struct('title', sprintf('Noise %s = %.1f', sig, R.examples(k).sigma), ...
            'xlabel', 'Compartment damage volume', 'ylabel', 'Simulated cognitive performance', 'rotTicks', {{}});
    end

    ax(3) = nexttile(tl);
    hold on;
    S = signedNlogp(R);
    scatter(R.sigma, S(1, :), 4, blue, 'filled', 'MarkerFaceAlpha', 0.5);
    scatter(R.sigma, S(2, :), 6, green, '^', 'filled', 'MarkerFaceAlpha', 0.5);
    yline(thresh, 'k:', 'LineWidth', 1.5);
    yline(-thresh, 'k:', 'LineWidth', 1.5);
    yline(0, '-', 'Color', [0.6 0.6 0.6]);
    xlim([R.params.nlow R.params.nhigh]);
    legend({'Slope comparison', 'with lesion volume covariate', 'Significance threshold'}, ...
        'Location', 'northeast', 'Box', 'off');
    info{3} = struct('title', 'Slope comparison vs noise', 'xlabel', ['Noise (' sig ')'], ...
        'ylabel', 'signed -ln(p)  (+ = H-implicated)', 'rotTicks', {{}});

    ax(4) = nexttile(tl);
    [bandTop, bandBottom, rotTicks] = splitViolins(R, thresh, blue, green);
    info{4} = struct('title', 'Effect of lesion volume adjustment', 'xlabel', ['Noise bin (' sig ')'], ...
        'ylabel', 'signed -ln(p)', 'rotTicks', {rotTicks});

    drawnow;
    for k = 1:4
        ax(k).FontSize = 10;
        placeAxesText(ax(k), info{k});
    end
    placeViolinPercentages(ax(4), bandTop, bandBottom);

    outFile = [matFile(1:end - 4) '_figure.pdf'];
    exportgraphics(fig, outFile, 'ContentType', 'vector');
    fprintf('Saved %s\n', outFile);
end


function exampleScatter(ex, hyp)
    tbl = ex.tbl;
    isTarget = tbl.Compartment == '1';
    hold on;
    scatter(tbl.Overlap(~isTarget), tbl.Behaviour(~isTarget), 8, 'k', 'filled');
    scatter(tbl.Overlap(isTarget), tbl.Behaviour(isTarget), 8, 'r', 'filled');
    xr = [0 300];
    cT = polyfit(tbl.Overlap(isTarget), tbl.Behaviour(isTarget), 1);
    cN = polyfit(tbl.Overlap(~isTarget), tbl.Behaviour(~isTarget), 1);
    plot(xr, polyval(cT, xr), 'r-', 'LineWidth', 1.5);
    plot(xr, polyval(cN, xr), 'k-', 'LineWidth', 1.5);
    xlim(xr); ylim([-2.5 2.5]);
    text(0.55, 0.45, sprintf('p=%.3g', ex.p(1)), 'Units', 'normalized');
    legend({['non-' hyp ' compartment'], [hyp ' compartment'], [hyp ' fit'], ['non-' hyp ' fit']}, ...
        'Location', 'northeast', 'Box', 'off');
end


function S = signedNlogp(R)
    % -ln(p), positive when the interaction implicates H
    S = R.nlogp .* (sign(R.beta) * R.params.expectedSign);
end


function [bandTop, bandBottom, rotTicks] = splitViolins(R, thresh, blue, green)
    edges = round(R.params.nlow:R.params.binWidth:R.params.nhigh, 10);
    bin = discretize(R.sigma, edges);
    nBins = numel(edges) - 1;
    S = signedNlogp(R);
    hold on;
    % violin range: at least 1.5 beyond the +-threshold lines, with a band for the percentages outside it
    yMax = max(prctile(S(:), 99.9), thresh + 1.5);
    yMin = min(prctile(S(:), 0.1), -thresh - 1.5);
    band = 0.12 * (yMax - yMin);
    bandTop = struct('x', [], 'y', yMax + band / 2, 'str', {{}}, 'color', zeros(0, 3));
    bandBottom = struct('x', [], 'y', yMin - band / 2, 'str', {{}}, 'color', zeros(0, 3));
    for b = 1:nBins
        v1 = S(1, bin == b)';
        v2 = S(2, bin == b)';
        [f1, y1] = ksdensity(v1);
        [f2, y2] = ksdensity(v2);
        f1 = f1(:); y1 = y1(:); f2 = f2(:); y2 = y2(:);
        s = 0.42 / max([f1; f2]);
        m1 = median(v1);
        m2 = median(v2);
        w1 = s * interp1(y1, f1, m1);
        w2 = s * interp1(y2, f2, m2);
        % bodies are cropped to the data range, keeping the percentage bands clear
        k1 = y1 >= yMin & y1 <= yMax;
        k2 = y2 >= yMin & y2 <= yMax;
        % half-transparent body with opaque outline and median in the series colour;
        % the median spans only the body
        patch([b; b - s * f1(k1); b], [y1(find(k1, 1)); y1(k1); y1(find(k1, 1, 'last'))], blue, ...
            'FaceAlpha', 0.5, 'EdgeColor', blue, 'EdgeAlpha', 1, 'LineWidth', 1);
        patch([b; b + s * f2(k2); b], [y2(find(k2, 1)); y2(k2); y2(find(k2, 1, 'last'))], green, ...
            'FaceAlpha', 0.5, 'EdgeColor', green, 'EdgeAlpha', 1, 'LineWidth', 1);
        plot([b - w1, b], m1 * [1 1], '-', 'Color', blue, 'LineWidth', 2.5);
        plot([b, b + w2], m2 * [1 1], '-', 'Color', green, 'LineWidth', 2.5);
        vs = {v1, v2};
        cs = {blue, green};
        for k = 1:2
            x = b + (k - 1.5) * 0.48;
            bandTop.x(end + 1) = x;
            bandTop.str{end + 1} = sprintf('%.0f%%', 100 * mean(vs{k} > thresh));
            bandTop.color(end + 1, :) = cs{k} * 0.7;
            bandBottom.x(end + 1) = x;
            bandBottom.str{end + 1} = sprintf('%.0f%%', 100 * mean(vs{k} < -thresh));
            bandBottom.color(end + 1, :) = cs{k} * 0.7;
        end
    end
    yline(thresh, 'k:', 'LineWidth', 1.5);
    yline(-thresh, 'k:', 'LineWidth', 1.5);
    yline(0, '-', 'Color', [0.6 0.6 0.6]);
    ylim([yMin - band, yMax + band]); xlim([0.4 nBins + 0.6]);
    yt = get(gca, 'YTick');
    set(gca, 'YTick', yt(yt >= yMin & yt <= yMax));
    xticks(1:nBins);
    xticklabels(repmat({''}, 1, nBins));
    rotTicks = arrayfun(@(i) sprintf('%.2f-%.2f', edges(i), edges(i + 1)), 1:nBins, 'UniformOutput', false);
end


function placeViolinPercentages(ax, bandTop, bandBottom)
    fs = 7;
    [wpt, hpt] = axesSizePt(ax);
    dx = diff(xlim(ax)) / wpt;
    dy = diff(ylim(ax)) / hpt;
    for set_ = {bandTop, bandBottom}
        B = set_{1};
        w = textWidths(B.str, fs * ones(size(B.str)), false(size(B.str)));
        for i = 1:numel(B.str)
            text(ax, B.x(i) - w(i) / 2 * dx, B.y - 0.36 * fs * dy, B.str{i}, 'FontSize', fs, 'Color', B.color(i, :), ...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'Clipping', 'off');
        end
    end
end


function placeAxesText(ax, info)
    % title, axis labels and rotated x tick labels
    tsz = 11;  lsz = 11;  ksz = ax.FontSize;
    YLABEL_GAP = 13;   % pt from the widest y tick label to the label baseline
    [wpt, hpt] = axesSizePt(ax);
    xl = xlim(ax);  yl = ylim(ax);
    dx = diff(xl) / wpt;           % data units per point
    dy = diff(yl) / hpt;
    xmid = mean(xl);  ymid = mean(yl);

    yts = ax.YTickLabel;
    if isempty(yts), yts = {'0'}; end
    wy = max(textWidths(yts(:)', ksz * ones(1, numel(yts)), false(1, numel(yts))));
    if isempty(info.rotTicks)
        xDepth = ksz + 4;       % short numeric tick labels: x label sits right below them
    else
        wr = textWidths(info.rotTicks, ksz * ones(1, numel(info.rotTicks)), false(1, numel(info.rotTicks)));
        xDepth = max(wr) * sind(45) + ksz * cosd(45) + 8;
        xt = ax.XTick;
        for i = 1:numel(info.rotTicks)     % labels end just below their tick (45 degrees)
            ex = xt(i) + 2 * dx;
            ey = yl(1) - 8 * dy;
            text(ax, ex - wr(i) * cosd(45) * dx, ey - wr(i) * sind(45) * dy, info.rotTicks{i}, 'FontSize', ksz, ...
                'Rotation', 45, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'Clipping', 'off');
        end
    end

    w = textWidths({info.title, info.xlabel, info.ylabel}, [tsz lsz lsz], [true false false]);
    text(ax, xmid - w(1) / 2 * dx, yl(2) + 8 * dy, info.title, 'FontSize', tsz, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'Clipping', 'off');
    text(ax, xmid - w(2) / 2 * dx, yl(1) - (xDepth + 6 + 0.75 * lsz) * dy, info.xlabel, 'FontSize', lsz, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'Clipping', 'off');
    text(ax, xl(1) - (wy + YLABEL_GAP) * dx, ymid - w(3) / 2 * dy, info.ylabel, 'FontSize', lsz, 'Rotation', 90, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'Clipping', 'off');
end


function [wpt, hpt] = axesSizePt(ax)
    p = getpixelposition(ax, true);
    f = 72 / get(groot, 'ScreenPixelsPerInch');
    wpt = p(3) * f;
    hpt = p(4) * f;
end


function w = textWidths(strs, sizes, bold)
    % widths in points from the Helvetica glyph advance widths (1/1000 em)
    %   characters 32-126, then sigma (963) and the minus sign (8722)
    REG=[278 278 355 556 556 889 667 191 333 333 389 584 278 333 278 278 556 556 556 556 556 556 556 556 556 556 278 278 584 584 584 556 1015 667 667 722 722 667 611 778 722 278 500 667 556 833 722 778 667 778 722 667 611 722 667 944 667 667 611 278 278 278 469 556 333 556 556 500 556 556 278 556 556 222 222 500 222 833 556 556 556 556 333 500 278 556 500 722 500 500 500 334 260 334 584];
    BOLD=[278 333 474 556 556 889 722 238 333 333 389 584 278 333 278 278 556 556 556 556 556 556 556 556 556 556 333 333 584 584 584 611 975 722 722 722 722 667 611 778 722 278 556 722 611 833 722 778 667 778 722 667 611 722 667 944 667 667 611 333 278 333 584 556 333 556 611 556 611 556 333 611 611 278 278 556 278 889 611 611 611 611 389 556 333 611 556 778 556 556 500 389 280 389 584];
    EXTRA_REG=[617 584];
    EXTRA_BOLD=[684 584];
    w = zeros(1, numel(strs));
    for i = 1:numel(strs)
        c = double(strs{i});
        if bold(i), a = BOLD; e = EXTRA_BOLD; else, a = REG; e = EXTRA_REG; end
        adv = zeros(size(c));
        isAscii = c >= 32 & c <= 126;
        adv(isAscii) = a(c(isAscii) - 31);
        adv(c == 963) = e(1);
        adv(c == 8722) = e(2);
        if any(~(isAscii | c == 963 | c == 8722)), error('Unsupported character in "%s"', strs{i}); end
        w(i) = sum(adv) / 1000 * sizes(i);
    end
end
