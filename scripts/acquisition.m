function [u, y, sel, sig, realizations, power_levels] = acquisition(filename)
% [u, y, sel, sig, realizations, power_levels] = acquisition(fs)
%
% Load the acquired input and output signals from .mat files
%
% Inputs:
%   fs - Sampling frequency (Hz)
% Outputs:
%   u - Input signal tensor (samples x realizations x power levels)
%   y - Output signal tensor (samples x realizations x power levels)
%   sel - excited bins
%   sig - input signal
%   realizations - Number of realizations, determined using available files
%   power_levels - Number of power levels, determined using available files
%
% It assumes the measurements name format is:
% 'out{fs in kHz}k_ACQ_R{realization}_P{power level}_E0_M0_F0.mat'


    resultsPath = '../results/';
    [realizations, power_levels] = detect_counts(resultsPath, filename);

    % load a single file to get N
    file = strcat(resultsPath, filename, "ACQ_R0_P0_E0_M0_F0.mat");
    load(file);
    N = length(YR0);

    % u and y are tensor with following dimensions:
    % (sample, realization, power level)
    u = zeros(N, realizations, power_levels);
    y = zeros(N, realizations, power_levels);

    % load measured signals
    for r = 0:realizations-1
        for p = 0:power_levels-1
            file = strcat(resultsPath, filename, "ACQ_R", num2str(r), "_P", num2str(p), "_E0_M0_F0.mat");
            load(file);
            u(:, r+1, p+1) = YR0(:);
            y(:, r+1, p+1) = YR1(:);
        end
    end

    % load sig
    sigFile = strcat(resultsPath, filename, "AWG_R0_E0.mat");
    load(sigFile);
    sig = SigR0;

    % load sel
    selFile = strcat(resultsPath, filename, "AWGEXC_R0_E0.mat");
    load(selFile);
    sel = SelExc0;

end

function [realizations, power_levels] = detect_counts(resultsPath, filename)
% Detect number of realizations and power levels by probing for existing files.
% Returns counts (non-negative integers). Filenames checked:
% [resultsPath]/[filename]ACQ_R{r}_P{p}_E0_M0_F0.mat
%
% Example: [realizations, power_levels] = detect_counts('../results/', 'out2k');

    % detect realizations by probing P0 files starting at R0
    r = 0;
    while true
        fname = fullfile(resultsPath, strcat(filename, "ACQ_R", num2str(r), "_P0_E0_M0_F0.mat"));
        if isfile(fname)
            r = r + 1;
        else
            break;
        end
    end
    realizations = r; % may be 0 if no files found

    % detect power levels by probing R0 files starting at P0 (only if at least one realization)
    p = 0;
    if realizations > 0
        while true
            fname = fullfile(resultsPath, strcat(filename, "ACQ_R0_P", num2str(p), "_E0_M0_F0.mat"));
            if isfile(fname)
                p = p + 1;
            else
                break;
            end
        end
    end
    power_levels = p; % may be 0 if no files found
end

