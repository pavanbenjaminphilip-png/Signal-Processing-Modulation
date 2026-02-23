# Signal Processing Modulation Processes

This repository contains simulations and implementations of digital modulation schemes, specifically focusing on **8-QAM (Quadrature Amplitude Modulation)** and channel modeling using **AWGN (Additive White Gaussian Noise)**.

## Project Overview
This project demonstrates the complete signal flow of a digital communication system:
1.  **Digital Input**: Generation of binary data.
2.  **Vector Encoding**: Mapping bits to constellation points.
3.  **Carrier Modulation**: Implementing coherent modulation (I & Q carriers).
4.  **Channel Modeling**: Simulating signal attenuation and AWGN.
5.  **Demodulation**: Correlator-based coherent demodulation.
6.  **Detection**: Maximum Likelihood (ML) detection using Euclidean distance.
7.  **Performance Analysis**: Bit Error Rate (BER) calculation and visualization.

## Contents
-   `QAM.m`: A comprehensive MATLAB script for 8-QAM software simulation, including performance analysis and detailed visualization of constellation diagrams and time-domain signals.
-   `awgn_server.py`: A Python-based server designed to handle carrier-modulated signals, add controllable AWGN noise, and facilitate communication (useful for hardware-in-the-loop simulations, e.g., with Raspberry Pi).
-   `qam8_awgn_wifi_2025.m`: Refined 8-QAM simulation script.

## Key Features
-   **High Resolution**: Simulates high-frequency carriers (e.g., 10 MHz) with high sampling rates.
-   **Realistic Channel**: Includes both attenuation factor ($h$) and SNR-based noise.
-   **Visual Excellence**: Generates multiple plots including:
    -   Digital bitstreams (Input vs. Output).
    -   Modulated carrier waveforms.
    -   Signal constellations (Ideal, Transmitted, Received, and Detected).
    -   Error distribution analysis.

## Getting Started
### MATLAB Simulation
1.  Open `QAM.m` in MATLAB.
2.  Run the script.
3.  Enter your binary input sequence (multiple of 3 bits) when prompted.
4.  Specify the desired SNR (dB) and channel attenuation ($h$).

### Python AWGN Server
1.  Run `python awgn_server.py`.
2.  The server listens on port 5000 for incoming signal samples to process with noise.

## Author
**Pavan Benjamin Philip**
