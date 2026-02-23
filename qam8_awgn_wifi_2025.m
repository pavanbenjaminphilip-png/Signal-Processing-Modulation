% qam8_awgn_raspberry_PI
% 8-QAM Communication System with Raspberry Pi AWGN Channel
% ======================================================================
% 8-QAM Communication System with Raspberry Pi AWGN Channel
% ======================================================================
% PURPOSE:
% This script implements a complete 8-QAM digital communication system
% with carrier modulation. MATLAB performs modulation, demodulation,
% detection, visualization, and BER analysis, while a Raspberry Pi acts
% as a real-time AWGN channel via TCP/IP.
%
% SYSTEM FLOW:
% Digital bits → 8-QAM mapping → Carrier modulation → TCP transmission
% → AWGN (Raspberry Pi) → Channel attenuation → Demodulation → ML detection
% → BER analysis → Constellation & time-domain visualization
clc;
clear all;
close all;

%% Configuration
% ----------------------------------------------------------------------
% Raspberry Pi network configuration for TCP communication
PI_IP = '192.168.0.178';  % Raspberry Pi IP address
PI_PORT = 5000;           % TCP port for communication

% User inputs
fprintf('============================================\n');
fprintf('8-QAM Communication System (PDE2103 Format)\n');
fprintf('============================================\n\n');

% Input binary sequence from user (default provided)
x = input('Enter bit sequence (Default:[ 0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 ]): ', 's');
if isempty(x)
    x = [ 0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 ];
else
    x = x - '0'; % Convert char to numeric array
    N_BITS = length(x);
    
    % Pad bits if not divisible by 3 (8-QAM: 3 bits per symbol)
    while (mod(N_BITS,3))
        disp('The number of bits should be divisible by 3 for 8-QAM.');
        x = input('Enter Digital Input Information in an array form (multiple of 3 bits), e.g., [ 0 0 0 0 0 1 0 1 0  0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 ], = ');
        N_BITS = length(x);
    end
end

% SNR input
SNR_dB = input('SNR in dB (default 9): ');
if isempty(SNR_dB), SNR_dB = 9; end

% Channel attenuation input
h_attenuation = input('Channel attenuation h (0 < h ≤ 1, default 0.7): ');
if isempty(h_attenuation), h_attenuation = 0.7; end
if h_attenuation <= 0 || h_attenuation > 1
    error('Channel attenuation h must be in range (0, 1]!');
end

fprintf('\n');

%% Timing Parameters
% ----------------------------------------------------------------------
% Defines bit rate, symbol rate, and symbol duration
Tb = 0.000001;          % Bit period (seconds)
Rb = 1/Tb;              % Bit rate (bps)
Rs = Rb/3;              % Symbol rate (symbols/sec), 8-QAM = 3 bits/symbol
Ts = 3*Tb;              % Symbol period (time for 1 symbol)

% Display system parameters
fprintf('============================================\n');
fprintf('SYSTEM PARAMETERS\n');
fprintf('============================================\n');
fprintf('Bit Period (Tb):   %.2e seconds\n', Tb);
fprintf('Bit Rate (Rb):        %.2e bps (%.1f Mbps)\n', Rb, Rb/1e6);
fprintf('Symbol Rate (Rs):     %.2e symbols/sec (%.2f Msps)\n', Rs, Rs/1e6);
fprintf('Symbol Period (Ts): %.2e seconds\n', Ts);
fprintf('============================================\n\n');

%% Digital Input
% Ensure number of bits is multiple of 3
while mod(length(x), 3) ~= 0
    x = x(1:end-1);
end
L = length(x);
% Displays transmitted binary data
disp('Binary Input Information at Transmitter:');
if L <= 30
    disp(x);
else
    fprintf('First 30 bits: ');
    disp(x(1:30));
    fprintf('... (%d total bits)\n', L);
end

%% Vector Encoder (8-QAM constellation)
% Define 8-QAM symbol points (I-Q plane)
X_0 = [-1; -3];  % bits [0 0 0] → X_0
X_1 = [-1; -1];  % bits [0 0 1] → X_1
X_2 = [-1; +3];  % bits [0 1 0] → X_2
X_3 = [-1; +1];  % bits [0 1 1] → X_3
X_4 = [+1; -3];  % bits [1 0 0] → X_4
X_5 = [+1; -1];  % bits [1 0 1] → X_5
X_6 = [+1; +3];  % bits [1 1 0] → X_6
X_7 = [+1; +1];  % bits [1 1 1] → X_7

% Constellation matrix
C = [X_0, X_1, X_2, X_3, X_4, X_5, X_6, X_7];

% Normalize to unit average energy
C_normalized = C / sqrt(6);

% Bit to Symbol Mapping (8-QAM Modulation)
% Groups bits into symbols and maps them to constellation points
X = [];
for n2 = 1:3:L  % 3 bits per symbol
    if x(n2:n2+2) == [0 0 0]
        X = [X, C_normalized(:,1)]; 
    elseif x(n2:n2+2) == [0 0 1]
        X = [X, C_normalized(:,2)]; 
    elseif x(n2:n2+2) == [0 1 0]
        X = [X, C_normalized(:,3)]; 
    elseif x(n2:n2+2) == [0 1 1]
        X = [X, C_normalized(:,4)]; 
    elseif x(n2:n2+2) == [1 0 0]
        X = [X, C_normalized(:,5)]; 
    elseif x(n2:n2+2) == [1 0 1]
        X = [X, C_normalized(:,6)]; 
    elseif x(n2:n2+2) == [1 1 0]
        X = [X, C_normalized(:,7)]; 
    elseif x(n2:n2+2) == [1 1 1]
        X = [X, C_normalized(:,8)]; 
    end
end

% Number of symbols
[r2, c2] = size(X);  
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

% Total transmission time
T_transmission = N_symbols * Ts;
fprintf('✓ Total transmission time: %.2e seconds (%.2f µs)\n', T_transmission, T_transmission*1e6);

%% Carrier Wave Modulation
% Uses orthogonal carriers(2 basis functions) for I/Q modulation
fprintf('\n--- Carrier Wave Modulation ---\n');

Ac = 1;                % Carrier amplitude
Fc = 10^7;             % Carrier frequency
nb = 100;              % Number of samples per symbol
t2 = Ts/nb:Ts/nb:Ts;   % Time vector for 1 symbol

% Orthogonal carriers
phi_1 = cos(2*pi*Fc*t2); % In-phase carrier (1st basis function)
phi_2 = sin(2*pi*Fc*t2); % Quadrature carrier (2nd basis function)

mod = [];
for k1 = 1:1:c2
    % I*phi1 + Q*phi2 (modulation)
    y = X(1,k1)*Ac*phi_1 + X(2,k1)*Ac*phi_2;
    mod = [mod y];
end

t3 = Ts/nb:Ts/nb:Ts*c2; % Total time vector
fprintf(' Modulated %d symbols with carrier (Fc = %.2e Hz)\n', N_symbols, Fc);

x_modulated = mod;  % Store modulated signal

%% Transmit to Raspberry Pi
% Sends modulated samples to Raspberry Pi and receives noisy samples
fprintf('\n--- Transmitting to Raspberry Pi ---\n');
fprintf('Connecting to %s:%d...\n', PI_IP, PI_PORT);

try
    % Create TCP client and connect
    tcpClient = tcpclient(PI_IP, PI_PORT, 'Timeout', 10);
    fprintf('Connected!\n');
    
    n_samples = length(x_modulated);  % Number of samples
    % Header: n_samples + SNR
    header = [typecast(swapbytes(int32(n_samples)), 'uint8'), ...
              typecast(swapbytes(double(SNR_dB)), 'uint8')];
    write(tcpClient, header);
    
    % Convert modulated signal to byte stream
    tx_data = zeros(n_samples * 8, 1, 'uint8');
    for k = 1:n_samples
        idx = (k-1)*8 + 1;
        tx_data(idx:idx+7) = typecast(swapbytes(double(x_modulated(k))), 'uint8');
    end
    write(tcpClient, tx_data);  % Send signal
    fprintf('Sent %d samples (SNR = %.1f dB)\n', n_samples, SNR_dB);
    
    fprintf('  Waiting for AWGN processing...\n');
    rx_data = read(tcpClient, n_samples * 8, 'uint8'); % Receive noisy signal
    
    % Convert received bytes back to double
    y_received = zeros(n_samples, 1);
    for k = 1:n_samples
        idx = (k-1)*8 + 1;
        y_received(k) = swapbytes(typecast(rx_data(idx:idx+7), 'double'));
    end
    
    clear tcpClient;  % Close connection
    fprintf(' Received %d noisy samples\n', n_samples);
    
catch ME
    % Handle TCP errors and provide troubleshooting steps
    fprintf('\n ERROR: %s\n', ME.message);
    fprintf('\nTroubleshooting Steps:\n');
    fprintf('1. Verify Raspberry Pi IP address: %s\n', PI_IP);
    fprintf('   - On Pi, run: hostname -I\n');
    fprintf('2. Ensure Python server is running on Pi:\n');
    fprintf('   - Run: python3 awgn_server.py\n');
    fprintf('3. Check network connectivity:\n');
    fprintf('   - From your computer, run: ping %s\n', PI_IP);
    fprintf('4. Check firewall on Pi:\n');
    fprintf('   - Run: sudo ufw allow 5000/tcp\n');
    if exist('tcpClient', 'var')
        clear tcpClient;
    end
    fprintf('\n Exiting due to connection error.\n');
    return;
end

% Check if data received
if ~exist('y_received', 'var') || isempty(y_received)
    fprintf(' ERROR: No data received from Raspberry Pi!\n');
    return;
end

%% Apply Channel Attenuation
% Simulates signal path loss
fprintf('\n--- Applying Channel Attenuation ---\n');
y = h_attenuation * y_received;  % Scale received signal
fprintf(' Applied channel attenuation: h = %.2f\n', h_attenuation);

%% Carrier Wave Demodulation
% Demodulation using Matched Filter Method (correlation) that's basically integration
fprintf('\n--- Carrier Wave Demodulation ---\n');
s = length(t2);  % Samples per symbol
demod = [];

for n1 = s:s:length(y)
    % Integrate over one symbol period using the Trapezoidal Function
    z1 = trapz(t2, phi_1 .* y((n1-(s-1)):n1)'); % I-component
    z2 = trapz(t2, phi_2 .* y((n1-(s-1)):n1)'); % Q-component
% Normalization is done here as well because the amplitude undergoes various changes as 
% it goes through Carrier Modulation, AWGN (Raspberry Pi), Channel
% Attenuation and at the end integration through the receiver.So inorder to
% stabilise these changes we do Normalization
    rz1 = ((2*z1/Ts))/Ac;  % Normalize I
    rz2 = ((2*z2/Ts))/Ac;  % Normalize Q
    
    demod = [demod [rz1; rz2]]; % Append demodulated symbol
end

fprintf(' Demodulated %d symbols\n', size(demod, 2));
demod = demod / Ac;  % Further normalization

disp('Demodulated Symbol Vector at Receiver:');
if N_symbols <= 10
    disp(demod);
else
    fprintf('First 10 symbols (I,Q):\n');
    disp(demod(:,1:10));
    fprintf('... (%d total symbols)\n', N_symbols);
end

%% ML Detection
% Computes Euclidean distance to all constellation points 
fprintf('--- ML Detection ---\n');
m = [];
[r4, c4] = size(demod);

% Here the function named ML_detection_8QAM is being used here
%This the function was created and has been predifined below

for l1 = 1:c4
    m = [m, ML_detection_8QAM(demod(:,l1), C_normalized)];  % Map symbol → bits
end

disp('Binary Output Information at Receiver:');
if length(m) <= 30
    disp(m);
else
    fprintf('First 30 bits: '); disp(m(1:30));
    fprintf('... (%d total bits)\n', length(m));
end

%% Performance Analysis
% BER is calculated here
% This gives us an understanding of the effectiveness of this system 
x = x(:)'; m = m(:)';  % Ensure row vectors
num_errors = sum(x ~= m);  % Count bit errors
ber = num_errors / length(x);  % Bit Error Rate

fprintf('\n============================================\n');
fprintf('PERFORMANCE RESULTS\n');
fprintf('============================================\n');
fprintf('Modulation:           8-QAM with Carrier\n');
fprintf('Carrier Frequency:    %.2e Hz (%.1f MHz)\n', Fc, Fc/1e6);
fprintf('Symbol Rate (Rs):     %.2f Msps\n', Rs/1e6);
fprintf('Bit Rate (Rb):        %.1f Mbps\n', Rb/1e6);
fprintf('Channel attenuation:  h = %.2f\n', h_attenuation);
fprintf('SNR:                  %.1f dB\n', SNR_dB);
fprintf('Symbols transmitted:  %d\n', N_symbols);
fprintf('Bits transmitted:     %d\n', length(x));
fprintf('Transmission time:    %.2f µs\n', T_transmission*1e6);
fprintf('Bit errors:           %d\n', num_errors);
fprintf('BER:                  %.6f (%.4e)\n', ber, ber);
fprintf('System Accuracy:      %.2f%%\n', 100*(1-ber));
fprintf('============================================\n');

%% Visualization
% The 4 required plots according to the question that are
% Ideal Constellation
% Transmitted symbols
% Received Symbols
% Detected Constellation
color_ideal = [0.8, 0.2, 0.2];
color_tx = [0.2, 0.4, 0.8];
color_rx = [0.9, 0.6, 0.2];
color_detected = [0.2, 0.7, 0.3];

figure('Position', [50 50 1400 900], 'Color', 'w');

% Plot 1: Ideal Constellation
subplot(2, 3, 1);
hold on;

scatter(C_normalized(1,:), C_normalized(2,:), 200, 'o', 'filled', ...
    'MarkerFaceColor', color_ideal, 'MarkerEdgeColor', 'k', 'LineWidth', 2);

bit_patterns = {'000', '001', '010', '011', '100', '101', '110', '111'};
for k = 1:8
    text(C_normalized(1,k)+0.12, C_normalized(2,k)+0.12, bit_patterns{k}, ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0 0]);
end

plot([-1.2 0.2], [0 0], 'k--', 'LineWidth', 1);
plot([0 0], [-1.2 1.2], 'k--', 'LineWidth', 1);

grid on; axis equal;
xlim([-2.5 2.5]); ylim([-2 2]);
xlabel('In-Phase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quadrature (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ideal 8-QAM Constellation', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

% Plot 2: Transmitted Symbols
subplot(2, 3, 2);
hold on;

scatter(X(1,:), X(2,:), 20, 'o', 'filled', ...
    'MarkerFaceColor', color_tx, 'MarkerFaceAlpha', 0.5);

scatter(C_normalized(1,:), C_normalized(2,:), 150, 'p', 'filled', ...
    'MarkerFaceColor', color_ideal, 'MarkerEdgeColor', 'k', 'LineWidth', 2);

plot([-1.2 0.2], [0 0], 'k--', 'LineWidth', 1);
plot([0 0], [-1.2 1.2], 'k--', 'LineWidth', 1);

grid on; axis equal;
xlim([-2.5 2.5]); ylim([-2 2]);
xlabel('In-Phase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quadrature (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Transmitted Symbols (%d symbols)', N_symbols), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Transmitted', 'Ideal Points', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

% Plot 3: Received Symbols
subplot(2, 3, 3);
hold on;

scatter(demod(1,:), demod(2,:), 20, 'o', 'filled', ...
    'MarkerFaceColor', color_rx, 'MarkerFaceAlpha', 0.5);

scatter(C_normalized(1,:), C_normalized(2,:), 150, 'p', 'filled', ...
    'MarkerFaceColor', color_ideal, 'MarkerEdgeColor', 'k', 'LineWidth', 2);

C_attenuated = h_attenuation * C_normalized;
scatter(C_attenuated(1,:), C_attenuated(2,:), 100, 'x', ...
    'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);

plot([-1.5 0.5], [0 0], 'k--', 'LineWidth', 1);
plot([0 0], [-1.5 1.5], 'k--', 'LineWidth', 1);

grid on; axis equal;
xlim([-2.5 2.5]); ylim([-2 2]);
xlabel('In-Phase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quadrature (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Received Symbols(h=%.2f, SNR=%ddB)', h_attenuation, SNR_dB), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Received', 'Ideal', 'Attenuated', 'Location', 'best', 'FontSize', 9);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

% Plot 4: Detection Results
subplot(2,3,4);
hold on;

detected_symbols = zeros(2, length(m)/3);
for k = 1:length(m)/3
    bits = m((k-1)*3+1 : k*3);
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

scatter(demod(1,:), demod(2,:), 25, 'o', 'filled', ...
    'MarkerFaceColor', color_rx, 'MarkerFaceAlpha', 0.4);

scatter(detected_symbols(1,:), detected_symbols(2,:), 50, 'd', 'filled', ...
    'MarkerFaceColor', color_detected, 'MarkerEdgeColor', 'k', 'LineWidth', 1);

plot([-1.5 0.5], [0 0], 'k--', 'LineWidth', 1);
plot([0 0], [-1.5 1.5], 'k--', 'LineWidth', 1);

grid on; axis equal;
xlim([-2.5 2.5]); ylim([-2 2]);
xlabel('In-Phase (I)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quadrature (Q)', 'FontSize', 12, 'FontWeight', 'bold');
title('Received 8-QAM Constellation', 'FontSize', 14, 'FontWeight', 'bold');
legend('Received', 'Detected', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

% Plot 5: Bit Errors
subplot(2, 3, 5);
correct_bits = length(x) - num_errors;
bar_data = [correct_bits, num_errors];
b = bar(bar_data, 'FaceColor', 'flat');
b.CData(1,:) = color_detected;
b.CData(2,:) = color_ideal;

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
ylim([0 length(x)*1.1]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

% Plot 6: Performance Summary
subplot(2, 3, 6);
axis off;

summary_text = {
    '\fontsize{13}\bfPerformance Summary';
    '';
    sprintf('Modulation: 8-QAM');
    sprintf('Carrier: %.1f MHz', Fc/1e6);
    sprintf('Symbol Rate: %.2f Msps', Rs/1e6);
    sprintf('Bit Rate: %.1f Mbps', Rb/1e6);
    sprintf('Symbols: %d', N_symbols);
    sprintf('Bits: %d', length(x));
    sprintf('Tx Time: %.2f µs', T_transmission*1e6);
    '';
    sprintf('Channel: Raspberry Pi AWGN');
    sprintf('Attenuation: h = %.2f', h_attenuation);
    sprintf('SNR: %d dB', SNR_dB);
    '';
    sprintf('Bit Errors: %d', num_errors);
    sprintf('BER: %.5e', ber);
    sprintf('Accuracy: %.2f%%', 100*(1-ber));
};

text(0.5, 0.5, summary_text, 'FontSize', 11, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'Units', 'normalized', 'Interpreter', 'tex');

sgtitle(sprintf('8-QAM with Carrier Wave via Raspberry Pi '), 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\n✓ Visualization complete!\n');
fprintf('============================================\n');

%% Time Domain Signals Visualization with Digital I/O
% Here 4 other Plots are being made 
% Digital Input and Output plots using stair plot in matlab
% Analog waveform of Transmitted and Received Signals
figure('Position', [50 50 1400 900], 'Color', 'w');

t_display = t3 * 1e6; % Convert to microseconds

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
plot(t_display, x_modulated, 'Color', [0.2, 0.4, 0.8], 'LineWidth', 1.5);
grid on;
xlabel('Time (µs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
title('Transmitted Analog Signal (Modulated Carrier)', 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 max(t_display)]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
box on;

% Subplot 3: Received Analog Signal
subplot(4, 1, 3);
plot(t_display, y, 'Color', [0.9, 0.6, 0.2], 'LineWidth', 1.5);
grid on;
xlabel('Time (µs)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Received Analog Signal (h=%.2f, SNR=%ddB)', h_attenuation, SNR_dB), ...
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


%% ML Detection Function
% Here the ML_detection_8QAM function is defined
% This function is used to calculate the euclidean distance and map the
% symbols back to the constellations

function m = ML_detection_8QAM(X, C)
% Finds nearest constellation point to received symbol X
[r, c] = size(C);

for k1 = 1:c
    d(k1) = sqrt((X(1) - C(1,k1))^2 + (X(2) - C(2,k1))^2); % Euclidean distance
end

a = find(d == min(d));  % Index of closest constellation point
a = a(1);               % If multiple minima, take first

% Map index to 3-bit symbol
if a == 1, m = [0 0 0];
elseif a == 2, m = [0 0 1];
elseif a == 3, m = [0 1 0];
elseif a == 4, m = [0 1 1];
elseif a == 5, m = [1 0 0];
elseif a == 6, m = [1 0 1];
elseif a == 7, m = [1 1 0];
elseif a == 8, m = [1 1 1];
end

end
