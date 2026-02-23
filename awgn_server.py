Python 3.14.0 (tags/v3.14.0:ebf955d, Oct  7 2025, 10:15:03) [MSC v.1944 64 bit (AMD64)] on win32
Enter "help" below or click "Help" above for more information.
# awgn_server.py
# Updated for carrier wave modulation - handles time-domain samples
import socket
import numpy as np
import struct
... import sys
... 
... HOST = '0.0.0.0'  # Listen on all network interfaces
... PORT = 5000       # Port for communication
... 
... def add_awgn(signal, snr_db):
...     """Add AWGN to real-valued time-domain signal given SNR in dB"""
...     # Calculate signal power
...     signal_power = np.mean(signal ** 2)
...     
...     # Calculate noise power from SNR
...     snr_linear = 10 ** (snr_db / 10.0)
...     noise_power = signal_power / snr_linear
...     noise_std = np.sqrt(noise_power)
...     
...     # Generate AWGN (real-valued for time-domain signal)
...     noise = noise_std * np.random.randn(len(signal))
...     
...     return signal + noise
... 
... def main():
...     # Create TCP socket
...     with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
...         s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
...         s.bind((HOST, PORT))
...         s.listen(1)
...         print(f"Server listening on {HOST}:{PORT}")
...         print(f"Raspberry Pi IP: {socket.gethostbyname(socket.gethostname())}")
...         print("Ready to receive carrier-modulated signals...")
...         
...         while True:
...             try:
...                 conn, addr = s.accept()
...                 print(f"\nConnected by {addr}")
...                 
...                 with conn:
...                     while True:
...                         # Receive header: N_samples (int32) and SNR_dB (float64)
...                         header = conn.recv(12)
...                         if len(header) < 12:
                            print("Connection closed by client")
                            break
                        
                        N_samples = struct.unpack('!i', header[:4])[0]
                        snr_db = struct.unpack('!d', header[4:12])[0]
                        print(f"Received: N_samples={N_samples}, SNR={snr_db} dB")
                        
                        # Receive N_samples real values (each as float64)
                        data_size = N_samples * 8  # 8 bytes per double
                        data = b''
                        while len(data) < data_size:
                            packet = conn.recv(data_size - len(data))
                            if not packet:
                                break
                            data += packet
                        
                        if len(data) < data_size:
                            print("Incomplete data received")
                            break
                        
                        # Unpack real-valued signal samples
                        signal = np.zeros(N_samples, dtype=np.float64)
                        for i in range(N_samples):
                            signal[i] = struct.unpack('!d', data[i*8:(i+1)*8])[0]
                        
                        print(f"Signal power: {np.mean(signal**2):.6f}")
                        
                        # Add AWGN noise
                        rx_signal = add_awgn(signal, snr_db)
                        
                        print(f"Noisy signal power: {np.mean(rx_signal**2):.6f}")
                        
                        # Send back received signal
                        response = b''
                        for i in range(N_samples):
                            response += struct.pack('!d', rx_signal[i])
                        
                        conn.sendall(response)
                        print(f"Sent {N_samples} noisy samples back")
                        print("-" * 50)
                        
            except KeyboardInterrupt:
                print("\nServer shutting down...")
                break
            except Exception as e:
                print(f"Error: {e}")
                import traceback
                traceback.print_exc()
                continue

if __name__ == "__main__":
