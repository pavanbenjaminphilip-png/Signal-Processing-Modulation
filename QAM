%% <<<<<<<<<<<<<<<<<<< 8-QAM Software Simulation >>>>>>>>>>>>>>>>>>>
%The provided MATLAB code is a complete software-based 8-QAM simulation. 
% No comments or implementation blocks reference Raspberry Pi or hardware-based transmission.

% This script implements a complete 8-QAM digital communication system
% including modulation, AWGN channel, demodulation, ML detection,
% BER analysis, and extensive visualization.
clc
clear all;
close all;

%% ******************* Digital/Binary input information ********************
fprintf('========================================\n');
fprintf('   8-QAM Software Simulation System    \n');
fprintf('========================================\n\n');

% Accept binary input sequence from user (must be multiple of 3 bits)
x = input('Enter Digital Input Information in an array form (multiple of 3 bits);default:[ 0 0 0 0 0 1 0 1 0  0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 ] = ');

% Check if first input is empty
if isempty(x)
    fprintf('The input is empty. Using default input.\n');
    x = [ 0 0 0 0 0 1 0 1 0  0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 ];
end

L = length(x);

% Loop until valid input (multiple of 3) is provided
while (mod(L,3))
    disp('The number of bits should be divisible by 3 for 8-QAM.');
    x = input('Enter Digital Input Information in an array form (multiple of 3 bits);default:[ 0 0 0 0 0 1 0 1 0  0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 ] = ');
    
    % Check if re-entered input is empty (user pressed Enter)
    if isempty(x)
        fprintf('The input is empty. Using default input.\n');
        x = [ 0 0 0 0 0 1 0 1 0  0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 ];
    end
    
    L = length(x);
end
% Display transmitter input bits
disp('Binary Input Information at Transmitter: ');
disp(x);

%% Timing Parameters
% Define bit and symbol timing parameters
Tb = 0.000001;          % Bit period (time for one bit)
Rb = 1/Tb;              % Bit rate (bits per second)
Rs = Rb/3;              % Symbol rate (symbols per second)
                        % 8-QAM carries 3 bits per symbol (log2(8) = 3)
Ts = 3*Tb;              % Symbol period [time taken for 1 symbol (3 bits) ]
                 
% Display system timing parameters
fprintf('============================================\n');
fprintf('SYSTEM PARAMETERS\n');
fprintf('============================================\n');
fprintf('Bit Period (Tb):   %.2e seconds\n', Tb);
fprintf('Bit Rate (Rb):        %.2e bps (%.1f Mbps)\n', Rb, Rb/1e6);
fprintf('Symbol Rate (Rs):     %.2e symbols/sec (%.2f Msps)\n', Rs, Rs/1e6);
fprintf('Symbol Period (Ts): %.2e seconds\n', Ts);
fprintf('============================================\n\n');

%% ****************************Vector Encoder*******************************
%******  Defining the 8-point constellation ***************************
X_0 = [-1; -3];  % bits [0 0 0] map to X_0
X_1 = [-1; -1];  % bits [0 0 1] map to X_1
X_2 = [-1; +3];  % bits [0 1 0] map to X_2
X_3 = [-1; +1];  % bits [0 1 1] map to X_3
X_4 = [+1; -3];  % bits [1 0 0] map to X_4
X_5 = [+1; -1];  % bits [1 0 1] map to X_5
X_6 = [+1; +3];  % bits [1 1 0] map to X_6
X_7 = [+1; +1];  % bits [1 1 1] map to X_7

% Defining the 8-point constellation
C = [X_0, X_1, X_2, X_3, X_4, X_5, X_6, X_7];

% Normalize constellation to unit average energy
C_normalized = C / sqrt(6);

%% ******************* Bit to Symbol Mapping(Modulation) *******************
% Convert groups of 3 bits into corresponding 8-QAM symbols
X = [];
for n2 = 1:3:L  % Since we take 3 bits to map for 8-QAM, the step size is 3
    if x(n2:n2+2) == [0 0 0]
        X = [X, C_normalized(:,1)]; % Map to X_0
    elseif x(n2:n2+2) == [0 0 1]
        X = [X, C_normalized(:,2)]; % Map to X_1
    elseif x(n2:n2+2) == [0 1 0]
        X = [X, C_normalized(:,3)]; % Map to X_2
    elseif x(n2:n2+2) == [0 1 1]
        X = [X, C_normalized(:,4)]; % Map to X_3
    elseif x(n2:n2+2) == [1 0 0]
        X = [X, C_normalized(:,5)]; % Map to X_4
    elseif x(n2:n2+2) == [1 0 1]
        X = [X, C_normalized(:,6)]; % Map to X_5
    elseif x(n2:n2+2) == [1 1 0]
        X = [X, C_normalized(:,7)]; % Map to X_6
    elseif x(n2:n2+2) == [1 1 1]
        X = [X, C_normalized(:,8)]; % Map to X_7
    end
end

% Number of transmitted symbols
[r2, c2] = size(X);  % c2 is the number of symbols
N_symbols = c2;

disp('Symbol Vector at Transmitter:');
if N_symbols <= 10
    disp(X);
else
    fprintf('First 10 symbols (I,Q):\n');
    disp(X(:,1:10));
    fprintf('... (%d total symbols)\n', N_symbols);
end
fprintf('✓ Mapped %d bits to %d symbols\n', L, N_symbols);

% Calculate transmission time
T_transmission = N_symbols * Ts;
fprintf('✓ Total transmission time: %.2e seconds (%.2f µs)\n', T_transmission, T_transmission*1e6);

%% ==================== CARRIER WAVE MODULATION ====================
% Perform coherent carrier modulation using orthogonal cosine and sine carriers
fprintf('\n--- Applying Carrier Wave Modulation ---\n');

% Carrier parameters
Ac = 1;           % Carrier amplitude
Fc = 10^7;        % Carrier frequency (10 MHz)
nb = 100;         % Samples per symbol for time resolution

% Time vector for one symbol period
t_symbol = Ts/nb : Ts/nb : Ts;

% Basis functions (I and Q carriers)
phi_1 = cos(2*pi*Fc*t_symbol);      % I-carrier (cosine)
phi_2 = sin(2*pi*Fc*t_symbol);      % Q-carrier (sine)

% Modulate symbols with carrier
mod = [];
for k1 = 1:1:c2
    y = X(1,k1)*Ac*phi_1 + X(2,k1)*Ac*phi_2;   % Modulation: I*cos + Q*sin
    mod = [mod y];
end

t_mod = Ts/nb : Ts/nb : Ts*c2;      % Time vector for modulated signal

fprintf('✓ Carrier frequency (Fc): %.2e Hz (%.1f MHz)\n', Fc, Fc/1e6);
fprintf('✓ Carrier amplitude (Ac): %.2f\n', Ac);
fprintf('✓ Modulated signal generated (%d samples)\n', length(mod));

% Transmitted signal
tx_signal = mod;

%% ==================== CHANNEL: AWGN (Software) + ATTENUATION ====================
fprintf('\n--- Applying AWGN Channel (Software) ---\n');

% Ask for SNR
SNR_dB = input('Enter SNR in dB (e.g., 6): ');
if isempty(SNR_dB)
    SNR_dB = 6;
end

% Add AWGN noise to modulated signal ausing the specified SNR
sigma = sqrt(1.0 / (2.0 * 10^(SNR_dB/10.0)));
noise = sigma * randn(1, length(tx_signal));
rx_signal_noisy = tx_signal + noise;

fprintf('✓ Applied AWGN noise (SNR = %.1f dB) [Software]\n', SNR_dB);

%% ********************* Channel attenuation *****************************
h = input('Enter channel attenuation h (0 < h ≤ 1; default 0.7): ');
if isempty(h)
    h = 0.7;
end

% Apply channel attenuation to modulated signal
rx_signal_attenuated = h .* rx_signal_noisy;

fprintf('✓ Applied channel attenuation (h = %.2f)\n', h);

%% ==================== DEMODULATION ====================
% Correlator-based (refers to integration using Trapezoidal Functions) coherent demodulation
fprintf('\n--- Demodulation ---\n');

% Demodulate the received signal
s = length(t_symbol);
demod = [];
for n1 = s:s:length(rx_signal_attenuated)
    z1 = trapz(t_symbol, phi_1 .* rx_signal_attenuated((n1-(s-1)):n1));  % Correlate with I-carrier
    z2 = trapz(t_symbol, phi_2 .* rx_signal_attenuated((n1-(s-1)):n1));  % Correlate with Q-carrier
    
    % Normalize to remove carrier amplitude and integration scaling
    rz1 = ((2*z1/Ts)) / Ac;  % Normalize by Ts and Ac to recover I component
    rz2 = ((2*z2/Ts)) / Ac;  % Normalize by Ts and Ac to recover Q component
    demod = [demod [rz1; rz2]];
end

disp('Demodulated Symbol Vector at Receiver (after h and AWGN):');
if N_symbols <= 10
    disp(demod);
else
    fprintf('First 10 symbols (I,Q):\n');
    disp(demod(:,1:10));
    fprintf('... (%d total symbols)\n', N_symbols);
end

%% ==================== ML DETECTION ====================
% Maximum Likelihood detection using Euclidean distance
fprintf('--- ML Detection ---\n');

m = [];
[r4, c4] = size(demod);

for l1 = 1:c4
    % Call ML detection function for each received symbol
    m = [m, ML_detection_8QAM(demod(:,l1), C_normalized)];
end

disp('Binary Output Information at Receiver:');
if length(m) <= 30
    disp(m);
else
    fprintf('First 30 bits: ');
    disp(m(1:30));
    fprintf('... (%d total bits)\n', length(m));
end

%% ==================== PERFORMANCE ANALYSIS ====================
% Calculate Bit Error Rate
x = x(:)';  % Ensure row vector
m = m(:)';  % Ensure row vector

num_errors = sum(x ~= m);
ber = num_errors / length(x);

fprintf('\n========================================\n');
fprintf('PERFORMANCE RESULTS\n');
fprintf('========================================\n');
fprintf('Modulation:           8-QAM\n');
fprintf('Symbol Rate (Rs):     %.2f Msps\n', Rs/1e6);
fprintf('Bit Rate (Rb):        %.1f Mbps\n', Rb/1e6);
fprintf('Carrier Frequency:    %.2f MHz\n', Fc/1e6);
fprintf('Channel attenuation:  h = %.2f\n', h);
fprintf('SNR:                  %.1f dB\n', SNR_dB);
fprintf('Symbols transmitted:  %d\n', N_symbols);
fprintf('Bits transmitted:     %d\n', length(x));
fprintf('Transmission time:    %.2f µs\n', T_transmission*1e6);
fprintf('Bit errors:           %d\n', num_errors);
fprintf('BER:                  %.6f (%.4e)\n', ber, ber);
fprintf('System Accuracy:      %.2f%%\n', 100*(1-ber));
fprintf('========================================\n');

%% ==================== TIME DOMAIN SIGNALS WITH DIGITAL I/O ====================
% Here we have created stair plots for displaying the Digital Inputs and
% Outputs Signals
% Analog Waveforms of Received and Transmitted Signals has also been shown as
% well
figure('Position', [50 50 1400 900], 'Color', 'w');

t_display = t_mod * 1e6; % Convert to microseconds

% Create digital input signal (for display purposes - show limited bits)
num_bits_display = min(length(x), 24); % Limit to 24 bits for clarity
t_bit = linspace(0, num_bits_display * Tb * 1e6, num_bits_display * 10);
digital_input = zeros(size(t_bit));
for i = 1:num_bits_display
    idx_start = round((i-1) * 10 + 1);
    idx_end = round(i * 10);
    digital_input(idx_start:idx_end) = x(i);
end

% Create digital output signal (for display purposes)
digital_output = zeros(size(t_bit));
for i = 1:num_bits_display
    idx_start = round((i-1) * 10 + 1);
    idx_end = round(i * 10);
    digital_output(idx_start:idx_end) = m(i);
end

% Subplot 1: Digital Input Signal
subplot(4, 1, 1);
stairs(t_bit, digital_input, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 2);
grid on;
ylim([-0.2 1.2]);
xlabel('Time (µs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Value', 'FontSize', 12, 'FontWeight', 'bold');
title('Digital Input Signal (First 24 bits)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'YTick', [0 1]);
box on;

% Subplot 2: Transmitted Analog Signal
subplot(4, 1, 2);
plot(t_display, tx_signal, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 1.5);
grid on;
xlabel('Time (µs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
title('Transmitted Analog Signal (Modulated Carrier)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max(t_display)]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

% Subplot 3: Received Analog Signal
subplot(4, 1, 3);
plot(t_display, rx_signal_attenuated, 'Color', [0.9, 0.6, 0.2], 'LineWidth', 1.5);
grid on;
xlabel('Time (µs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Received Analog Signal (h=%.2f, SNR=%ddB)', h, SNR_dB), ...
    'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max(t_display)]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

% Subplot 4: Digital Output Signal
subplot(4, 1, 4);
stairs(t_bit, digital_output, 'Color', [0.2, 0.7, 0.3], 'LineWidth', 2);
hold on;
% Overlay input signal in lighter color for comparison
stairs(t_bit, digital_input, 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1, 'LineStyle', '--');
grid on;
ylim([-0.2 1.2]);
xlabel('Time (µs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Value', 'FontSize', 12, 'FontWeight', 'bold');
title('Digital Output Signal (Detected - First 24 bits)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Output', 'Input (ref)', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.5, 'YTick', [0 1]);
box on;

sgtitle('Complete Signal Flow: Digital Input → Analog Transmission → Digital Output', ...
    'FontSize', 16, 'FontWeight', 'bold');

%% ==================== VISUALIZATION ====================
% Color scheme
color_ideal = [0.8, 0.2, 0.2];         % Red for ideal points
color_tx = [0.2, 0.4, 0.8];            % Blue for transmitted
color_rx = [0.9, 0.6, 0.2];            % Orange for received
color_detected = [0.2, 0.7, 0.3];      % Green for detected
color_attenuated = [0.5, 0.5, 0.5];    % Gray for attenuated

figure('Position', [50 700 1400 900], 'Color', 'w');

%% Plot 1: Ideal 8-QAM Constellation
subplot(2,3,1);
hold on;

% Ideal constellation points
scatter(C_normalized(1,:), C_normalized(2,:), 200, 'o', 'filled', ...
    'MarkerFaceColor', color_ideal, 'MarkerEdgeColor', 'k', 'LineWidth', 2);

% Add bit labels
bit_patterns = {'000', '001', '010', '011', '100', '101', '110', '111'};
for k = 1:8
    text(C_normalized(1,k)+0.12, C_normalized(2,k)+0.12, bit_patterns{k}, ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0 0]);
end

% Reference axes
plot([-1.2 0.2], [0 0], 'k--', 'LineWidth', 1);
plot([0 0], [-1.2 1.2], 'k--', 'LineWidth', 1);

grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);
xlabel('In-Phase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quadrature (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ideal 8-QAM Constellation', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

%% Plot 2: Transmitted Symbols
subplot(2,3,2);
hold on;

% Transmitted symbols
scatter(X(1,:), X(2,:), 20, 'o', 'filled', ...
    'MarkerFaceColor', color_tx, 'MarkerFaceAlpha', 0.5);

% Ideal constellation
scatter(C_normalized(1,:), C_normalized(2,:), 150, 'p', 'filled', ...
    'MarkerFaceColor', color_ideal, 'MarkerEdgeColor', 'k', 'LineWidth', 2);

% Reference axes
plot([-1.2 0.2], [0 0], 'k--', 'LineWidth', 1);
plot([0 0], [-1.2 1.2], 'k--', 'LineWidth', 1);

grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);
xlabel('In-Phase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quadrature (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Transmitted Symbols (%d symbols)', N_symbols), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Transmitted', 'Ideal Points', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

%% Plot 3: Received Symbols (After AWGN and Attenuation)
subplot(2,3,3);
hold on;

% Received symbols (with noise and attenuation)
scatter(demod(1,:), demod(2,:), 20, 'o', 'filled', ...
    'MarkerFaceColor', color_rx, 'MarkerFaceAlpha', 0.5);

% Ideal constellation (for reference)
scatter(C_normalized(1,:), C_normalized(2,:), 150, 'p', 'filled', ...
    'MarkerFaceColor', color_ideal, 'MarkerEdgeColor', 'k', 'LineWidth', 2);

% Attenuated constellation (show effect of h)
C_attenuated = h * C_normalized;
scatter(C_attenuated(1,:), C_attenuated(2,:), 120, 'x', ...
    'MarkerEdgeColor', color_attenuated, 'LineWidth', 3);

% Reference axes
plot([-1.5 0.5], [0 0], 'k--', 'LineWidth', 1);
plot([0 0], [-1.5 1.5], 'k--', 'LineWidth', 1);

grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);
xlabel('In-Phase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quadrature (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Received Symbols (h=%.2f, SNR=%ddB)', h, SNR_dB), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Received', 'Ideal', 'Attenuated', 'Location', 'best', 'FontSize', 9);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

%% Plot 4: Detected Constellations
subplot(2,3,4);
hold on;

% Reconstruct detected symbols for visualization
detected_symbols = zeros(2, length(m)/3);
for k = 1:length(m)/3
    bits = m((k-1)*3+1 : k*3);
    % Find matching constellation point based on bit pattern
    if isequal(bits, [0 0 0])
        detected_symbols(:,k) = C_normalized(:,1);
    elseif isequal(bits, [0 0 1])
        detected_symbols(:,k) = C_normalized(:,2);
    elseif isequal(bits, [0 1 0])
        detected_symbols(:,k) = C_normalized(:,3);
    elseif isequal(bits, [0 1 1])
        detected_symbols(:,k) = C_normalized(:,4);
    elseif isequal(bits, [1 0 0])
        detected_symbols(:,k) = C_normalized(:,5);
    elseif isequal(bits, [1 0 1])
        detected_symbols(:,k) = C_normalized(:,6);
    elseif isequal(bits, [1 1 0])
        detected_symbols(:,k) = C_normalized(:,7);
    elseif isequal(bits, [1 1 1])
        detected_symbols(:,k) = C_normalized(:,8);
    end
end

% Received symbols
scatter(demod(1,:), demod(2,:), 25, 'o', 'filled', ...
    'MarkerFaceColor', color_rx, 'MarkerFaceAlpha', 0.4);

% Detected symbols
scatter(detected_symbols(1,:), detected_symbols(2,:), 50, 'd', 'filled', ...
    'MarkerFaceColor', color_detected, 'MarkerEdgeColor', 'k', 'LineWidth', 1);

% Reference axes
plot([-1.5 0.5], [0 0], 'k--', 'LineWidth', 1);
plot([0 0], [-1.5 1.5], 'k--', 'LineWidth', 1);

grid on; axis equal;
xlim([-2 2]); ylim([-2 2]);
xlabel('In-Phase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quadrature (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title('Received 8-QAM Constellation', 'FontSize', 14, 'FontWeight', 'bold');
legend('Received', 'Detected', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

%% Plot 5: Bit Errors
subplot(2,3,5);
correct_bits = length(x) - num_errors;
bar_data = [correct_bits, num_errors];
b = bar(bar_data, 0.6, 'FaceColor', 'flat');
b.CData(1,:) = color_detected;
b.CData(2,:) = color_ideal;
b.EdgeColor = 'k';
b.LineWidth = 1.5;

% Add value labels
if correct_bits > 0
    text(1, correct_bits/2, sprintf('%d\n(%.1f%%)', correct_bits, 100*correct_bits/length(x)), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'w');
end
if num_errors > 0
    text(2, num_errors/2, sprintf('%d\n(%.1f%%)', num_errors, 100*num_errors/length(x)), ...
        'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'w');
end

set(gca, 'XTickLabel', {'Correct', 'Errors'}, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Number of Bits', 'FontSize', 12, 'FontWeight', 'bold');
title('Bit Error Analysis', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim([0 length(x)*1.15]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

%% Plot 6: Performance Summary
subplot(2,3,6);
axis off;

% Background box
rectangle('Position', [0.1 0.1 0.8 0.8], 'FaceColor', [0.95 0.97 1], ...
          'EdgeColor', [0 0 0.6], 'LineWidth', 2);

summary_text = {
    '\fontsize{13}\bfPerformance Summary';
    '';
    sprintf('Modulation: 8-QAM');
    sprintf('Symbol Rate: %.2f Msps', Rs/1e6);
    sprintf('Bit Rate: %.1f Mbps', Rb/1e6);
    sprintf('Carrier Freq: %.1f MHz', Fc/1e6);
    sprintf('Symbols: %d', N_symbols);
    sprintf('Bits: %d', length(x));
    sprintf('Tx Time: %.2f µs', T_transmission*1e6);
    '';
    sprintf('Channel: AWGN (Software)');
    sprintf('Attenuation: h = %.2f', h);
    sprintf('SNR: %d dB', SNR_dB);
    '';
    sprintf('Bit Errors: %d', num_errors);
    sprintf('BER: %.5e', ber);
    sprintf('Accuracy: %.2f%%', 100*(1-ber));
};

% Ensure text is displayed properly
text(0.5, 0.5, summary_text, 'FontSize', 8, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Units', 'normalized', 'Interpreter', 'tex', ...
    'Color', [0 0 0.6], 'FontName', 'Arial');

%% Overall title
sgtitle(sprintf('8-QAM System '), 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\n✓ Visualization complete!\n');
fprintf('========================================\n');

%% ==================== ML DETECTION FUNCTION ====================
% Performs ML detection for 8-QAM using Euclidean distance
function m = ML_detection_8QAM(X, C)
% ML Detection for 8-QAM
% C is the constellation which is a 2x8 matrix
% X is the received symbol as [I; Q]
% Returns 3 bits m

[r, c] = size(C);  % C is 2x8

% Calculate distances between received symbol X and all constellation points
for k1 = 1:c
    d(k1) = sqrt((X(1) - C(1,k1))^2 + (X(2) - C(2,k1))^2);
end

% Find the location of the minimum distance
a = find(d == min(d));
a = a(1);  % Take first if there are ties

% Map the detected constellation point to 3 bits
if a == 1
    m = [0 0 0];  % X_0: bits [0 0 0]
elseif a == 2
    m = [0 0 1];  % X_1: bits [0 0 1]
elseif a == 3
    m = [0 1 0];  % X_2: bits [0 1 0]
elseif a == 4
    m = [0 1 1];  % X_3: bits [0 1 1]
elseif a == 5
    m = [1 0 0];  % X_4: bits [1 0 0]
elseif a == 6
    m = [1 0 1];  % X_5: bits [1 0 1]
elseif a == 7
    m = [1 1 0];  % X_6: bits [1 1 0]
elseif a == 8
    m = [1 1 1];  % X_7: bits [1 1 1]
end

end