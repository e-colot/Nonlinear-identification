% Computes the FRF estimate using the robust method.
% Based on multiple realization, each of which have multiple periods.
% The FRF is computed for each period (G_2D) and is averaged over them to build 1 FRF estimate for each realization (G_1D)
% The variance of G_2D equals the noise variance. The FRF are then averaged over the realizations, giving G_ML.
% The variance of G_1D gives a measure for the nonlinear distortions.
% showPlot = false disables figures
%

showPlot = 1;
dataFile = "robustMethod/full_5k";
fs = 5e3;

%% Load data
[u, y, ~, sig, realizations, ~] = acquisition(dataFile);

N = size(u, 1);
periodN = size(sig, 1); % number of samples of the original period
repNumber = N / periodN; % number of periods in the acquired signal

%% Remove transient

if showPlot
% transient visualization
    figure;
    plot(db(y(1:periodN, 1, 1) - y(periodN + 1:2*periodN, 1, 1)), "LineWidth", 2);
    xlabel('Samples');
    ylabel('Amplitude [dB]');
    grid on;
end

% transient removal
    transientPeriods = 1;
    assert(transientPeriods < repNumber, 'Transient periods to remove exceed total number of periods.');

    u = u(transientPeriods*periodN + 1:end, :, :);
    y = y(transientPeriods*periodN + 1:end, :, :);

% update sizes
    N = size(u, 1);
    repNumber = repNumber - transientPeriods;

%% --------------- FRF estimate ---------------

% preallocation
    G_2D = zeros(realizations, repNumber, periodN);
    G_1D = zeros(realizations, periodN);
    noiseVar_2D = zeros(realizations, periodN);

% for each realization
for r = 1:realizations
    % for each period
    for p = 1:repNumber
        u_per = u((p-1)*periodN + 1:p*periodN, r, 1);
        y_per = y((p-1)*periodN + 1:p*periodN, r, 1);
        
        % compute the FFT of both input and output
        U = fft(u_per);
        Y = fft(y_per);

        % Estimate the FRF
        G_2D(r, p, :) = Y ./ U;
    end

    % average noise out
    G_1D(r, :) = squeeze(mean(G_2D(r, :, :), 2));
    noiseVar_2D(r, :) = squeeze(var(G_2D(r, :, :), 0, 2)) / repNumber;
end

% average nonlinear distortions out
G_ML = squeeze(mean(G_1D, 1)).';
noise_var = mean(noiseVar_2D, 1).' / realizations;
distortion_var = var(G_1D, 0, 1).' / realizations;

total_var = noise_var + distortion_var;
f = (0:periodN-1).'*(fs/periodN);

%% Plots

if showPlot
    figure;
    plot(f, db(G_ML), 'o', 'LineWidth', 3);
    hold on;

    plot(f, db(sqrt(noise_var)), 'r', 'LineWidth', 1.5);
    plot(f, db(sqrt(distortion_var)), 'b', 'LineWidth', 1.5);
    plot(f, db(sqrt(total_var)), 'g', 'LineWidth', 1.5);
    %title('Variance Estimates');
    xlabel('Frequency [Hz]');
    ylabel('FRF [dB]');
    xlim([1/fs fs/4]);
    grid on;
    legend('G_{BLA}', 'Noise stddev', 'Distortion stddev', 'Total stddev');

end
