"""
This script performs a complete analysis of vibrating plates:
1. Plotting the vibration spectrum (to determine resonance frequencies, e.g., for model testing) using a Gaussian initial condition.
2. Plotting the center's vibration over time (to visually observe resonance).
3. Plotting a Chladni pattern.
4. Animating the 3D dynamics of the plate.
"""
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.fft import fft, fftfreq

# Global parameters
nx, ny = 101, 101
Lx, Ly = 1.0, 1.0
z_lim = 0.0001

# Load FFT & Time series data
data_fft = np.loadtxt("oscillation_centre.txt")
t = data_fft[:, 0]
u = data_fft[:, 1]

# Load Animation data
raw_anim = np.loadtxt("profil_anim.txt")
n_frames = int(raw_anim.shape[0] / nx)
data_anim = raw_anim.reshape((n_frames, nx, ny))

# Load Chladni data
data_chladni = np.loadtxt('enveloppe_finale.txt')
X_chladni = data_chladni[:, 0].reshape(nx, ny)
Y_chladni = data_chladni[:, 1].reshape(nx, ny)
Z_max = data_chladni[:, 2].reshape(nx, ny)

# Spatial display grids
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)


# FIGURE 1: Time Series and FFT
N = len(t)
dt = t[1] - t[0]

yf = fft(u)
xf = fftfreq(N, dt)[:N//2]
amplitude = 2.0 / N * np.abs(yf[0:N//2])

fig_fft, (ax_time, ax_fft) = plt.subplots(2, 1, figsize=(10, 8))

# 1. Time domain plot
ax_time.plot(t, u, color='r', linewidth=1)
ax_time.set_title("Center Displacement over Time")
ax_time.set_xlabel("Time (s)")
ax_time.set_ylabel("Amplitude (m)")
ax_time.grid(True)

# 2. Frequency domain plot
ax_fft.plot(xf, amplitude, color='b', linewidth=1.5)
ax_fft.set_title("Plate Vibration Spectrum")
ax_fft.set_xlabel("Frequency (Hz)")
ax_fft.set_ylabel("Amplitude")
ax_fft.set_xlim(0, 500)
ax_fft.grid(True)

plt.tight_layout()  # Adjust layout to prevent label overlap

# Print dominant frequency to console (ignoring the DC component at 0 Hz)
idx_max = np.argmax(amplitude[1:]) + 1
print(f"--> Main resonance frequency: {xf[idx_max]:.2f} Hz")


# FIGURE 2: CHLADNI PATTERN (ENVELOPE)
fig_chladni, ax_chladni = plt.subplots(figsize=(8, 7))
contour = ax_chladni.contourf(X_chladni, Y_chladni, Z_max, levels=50, cmap='magma')

# Nodal lines (where sand accumulates)
ax_chladni.contour(X_chladni, Y_chladni, Z_max, levels=[np.max(Z_max)*0.05], colors='cyan', linewidths=1.5)

fig_chladni.colorbar(contour, ax=ax_chladni, label="Oscillation Amplitude")
ax_chladni.set_title("Chladni Pattern")
ax_chladni.set_xlabel("X (m)")
ax_chladni.set_ylabel("Y (m)")
ax_chladni.set_aspect('equal')


# FIGURE 3: 3D ANIMATION
fig_anim = plt.figure(figsize=(9, 7))
ax_anim = fig_anim.add_subplot(111, projection='3d')

def update(frame):
    ax_anim.clear()
    ax_anim.set_zlim(-z_lim, z_lim)
    ax_anim.set_title(f"Animation - Frame {frame}/{n_frames}")
    
    Z_t = data_anim[frame, :, :].T
    surf = ax_anim.plot_surface(X, Y, Z_t, cmap='RdBu_r', antialiased=False, shade=True)
    return surf,

ani = FuncAnimation(fig_anim, update, frames=n_frames, interval=20, blit=False)

plt.show()
