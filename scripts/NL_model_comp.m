function NL_model_comp(y_true, y_est, fs)
    % y_true: Measured output signal
    % y_est:  Model predicted output signal
    % fs:     Sampling frequency (Hz)
    
    N = length(y_true);
    t = (0:N-1)/fs;
    
    % --- 1. Calculate PSD and Coherence-like metric ---
    [pxx_true, f] = periodogram(y_true, rectwin(N), N, fs);
    [pxx_est, ~]  = periodogram(y_est, rectwin(N), N, fs);
    
    % Calculate the Error Spectrum (Residuals)
    res = y_true - y_est;
    [pxx_err, ~]  = periodogram(res, rectwin(N), N, fs);
    
    % --- 2. Calculate BLA-style Gain ---
    % In a true BLA we'd use Cross-PSD (Pxy/Pxx). 
    % Here we compare the Estimated Transfer Function magnitude.
    H_true = fft(y_true) ./ fft(y_true); % Reference
    % We estimate the "Model Gain" relative to the truth
    H_model = fft(y_est) ./ fft(y_true); 
    
    nmse_db = 10 * log10(sum(res.^2) / sum(y_true.^2));

    figure('Color', 'w', 'Position', [100, 100, 900, 800]);

    % TOP: Time Domain Zoom
    subplot(211);
    plot(t, y_true, 'k', 'LineWidth', 1.2); hold on;
    plot(t, y_est, 'r--', 'LineWidth', 1);
    title(['Time Domain Comparison (NMSE: ', num2str(nmse_db, '%.2f'), ' dB)']);
    ylabel('Amplitude'); xlabel('Time (s)');
    legend('True', 'Model'); grid on;
    %xlim([t(1), t(min(N, round(fs*0.05)))]); % Zoom into first 50ms by default

    % BOTTOM: PSD Comparison
    subplot(212);
    plot(f, 10*log10(pxx_true), 'k', 'LineWidth', 1.5); hold on;
    plot(f, 10*log10(pxx_est), 'r--', 'LineWidth', 1);
    plot(f, 10*log10(pxx_err), 'g', 'Color', [0.4, 0.7, 0.4]);
    title('Power Spectral Density');
    ylabel('PSD (dB/Hz)'); xlabel('Frequency (Hz)');
    legend('True', 'Model', 'Residual Error');
    grid on; axis tight;
end