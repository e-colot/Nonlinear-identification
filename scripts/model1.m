clear; close all; clc;

fs = 5e3;

  %% Data aquisition
    [uFull, yFull, ~, sig, ~, ~] = acquisition('fastMethod1/odd_5k');

    N = size(yFull, 1);
    periodN = size(sig, 1); % number of samples of the original period
    repNumber = N / periodN; % number of periods in the acquired signal

  %% Load a single period

    y = yFull(end - periodN + 1:end, 1);
    u = uFull(end - periodN + 1:end, 1);

    Norder = 30;
    Nlag = 150;

N = length(u);

%% 2. Construct the Regression Matrix (H)
% We need to solve Y = H * A
% H columns correspond to basis functions: u(n-q) * |u(n-q)|^(k-1)

H = memoryPolynomial(u, Norder, Nlag);

%% 3. LLS Estimation (Find Parameters)
% Use MATLAB's backslash operator for Linear Least Squares
coeffs = H \ y;

%% 4. Model Validation
% Calculate the estimated output
y_est = H * coeffs;

% Calculate Normalized Mean Square Error (NMSE) in dB
error_signal = y - y_est;
nmse_db = 10 * log10(sum(abs(error_signal).^2) / sum(abs(y).^2));

fprintf('Model Identification Complete.\n');
fprintf('NMSE: %.2f dB\n', nmse_db);

%% 5. Visualization


NL_model_comp(y_est, y, 5e3);

