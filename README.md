# Signal Processing & Communications: 8-QAM Modulation System

![MATLAB](https://img.shields.io/badge/MATLAB-Simulation-orange.svg) 
![Python](https://img.shields.io/badge/Python-AWGN_Server-blue.svg)
![Status](https://img.shields.io/badge/Status-Project_Submission_Complete-green.svg)

## 📡 Project Overview
This repository contains the complete implementation and simulation of a **Digital Communication System** developed for the Signal Processing and Communications module. The core focus is on **8-QAM (Quadrature Amplitude Modulation)**, utilizing coherent detection and Maximum Likelihood (ML) estimation in a noisy channel environment.

### 🛠 System Architecture
The system follows a classic communication pipeline:
1.  **Source Coding**: Binary data generation (user-defined or default).
2.  **Vector Encoding**: Mapping triplets of bits to a normalized 8-point constellation.
3.  **Modulation**: Coherent carrier modulation using orthogonal basis functions ($\phi_1 = \cos(2\pi f_c t)$ and $\phi_2 = \sin(2\pi f_c t)$).
4.  **Channel Modeling**: Simulation of **AWGN (Additive White Gaussian Noise)** and channel attenuation ($h$).
5.  **Demodulation**: Correlation-based coherent demodulation using trapezoidal integration.
6.  **ML Detection**: Optimal decision making based on minimum Euclidean distance in the signal space.
7.  **Analysis**: Comprehensive BER (Bit Error Rate) performance evaluation and signal visualization.

---

## 📂 Repository Contents
The project is structured into simulation scripts and communication utilities:

| File | Type | Description |
| :--- | :--- | :--- |
| `QAM.m` | MATLAB | **Primary Simulation Engine**. Handles the full 8-QAM lifecycle, from bit mapping to BER visualization. |
| `awgn_server.py` | Python | **Networked Noise Server**. A TCP-based utility that simulates a remote AWGN channel, allowing for hardware-aware testing. |
| `qam8_awgn_wifi_2025.m` | MATLAB | Refined simulation variant optimized for the 2025 technical specifications. |
| `Signal Processing And Communications Final Project.docx` | Documentation | Detailed project report covering theoretical background and results (available in submission zip). |

---

## 🚀 Key Technical Features
### 1. 8-QAM Constellation Mapping
The system maps 3-bit sequences to eight distinct points in the I-Q plane, normalized to ensure unit average energy.
*   **Modulation Order ($M$)**: 8
*   **Symbol Depth**: 3 bits/symbol

### 2. Carrier Modulation
Implements high-frequency carrier wave modulation:
*   **Carrier Frequency ($f_c$)**: 10 MHz
*   **Sampling Density**: 100 samples per symbol for high-resolution time-domain analysis.

### 3. Channel Simulation & ML Detection
The receiver implements an optimal **Maximum Likelihood (ML)** detector. The decision rule minimizes the Euclidean distance between the received vector $\mathbf{r}$ and the constellation points $\mathbf{s}_i$:
$$\hat{\mathbf{s}} = \arg \min_{\mathbf{s}_i \in \mathcal{C}} \|\mathbf{r} - \mathbf{s}_i\|^2$$

---

## 📊 Visualizations
The simulation generates a multi-panel dashboard providing insights into:
*   **Time-Domain Analysis**: Input bitstreams vs. reconstructed output signals.
*   **Constellation Diagrams**: Comparative views of Ideal, Transmitted, and Noisy (Received) constellations.
*   **Performance Metrics**: Automated Bit Error Rate (BER) calculation and Accuracy percentage.

---

## 🛠 Usage Instructions
### Running the MATLAB Simulation
1.  Launch `QAM.m` in MATLAB.
2.  Follow the prompts to enter a binary sequence (multiples of 3).
3.  Adjust parameters for **SNR (dB)** and **Attenuation ($h$)** to observe system robustness.

### Utilizing the AWGN Server
For remote channel simulation:
```bash
python awgn_server.py
```
*The server will listen on port 5000, ready to process incoming signal samples with additive noise.*

---

## 🏆 Author
**Pavan Benjamin Philip**
*Middlesex University - Signal Processing and Communications*
