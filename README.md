# Chladni Figures Simulation & Modal Analysis

This project simulates the dynamic behavior of a vibrating brass plate using the Finite Difference Time Domain (FDTD) method. It solves the biharmonic wave equation for different boundary conditions and initial setups to generate and visualize beautiful Chladni patterns.

A Python script is also provided to post-process the generated data, performing a Fast Fourier Transform (FFT) to find resonance frequencies, plotting Chladni figures, and rendering a 3D animation of the plate's dynamics.

## Features

The C++ simulator allows you to choose between 3 distinct setups:
1. **Simply Supported Plate**: Sinusoidal initial deformation with fixed boundaries (u=0).
2. **Free Plate (Gaussian Deformation)**: Free boundaries (Zero Bending Moment & Shear Force) with a Gaussian initial impact.
3. **Free Plate (Forced Regime)**: Free boundaries excited by a central periodic force (vibrator) to generate steady-state Chladni patterns.

## Physics & Parameters

The simulation uses standard mechanical properties for brass:
- Young's Modulus ($E$): $130 \times 10^9 \text{ Pa}$
- Plate Thickness ($h$): $0.002 \text{ m}$
- Density ($\rho$): $8500 \text{ kg/m}^3$
- Poisson's Ratio ($\nu$): $0.37$
- Grid Size: $101 \times 101$
- Time step ($dt$): Calculated dynamically using the Von Neumann stability criterion.

## Prerequisites

To compile and run this project, you will need:
- A C++ compiler
- Python 3.x
- Standard Python libraries: `numpy`, `matplotlib`, `scipy`
- Tkinter (for GUI rendering): `python3-tk` (install via your package manager on Linux, e.g., `sudo apt install python3-tk`)

## Installation & Usage

### 1. Clone the repository
```bash
git clone https://github.com/your-username/Chladni-Figures.git
cd Chladni-Figures
```

### 2. Compile the C++ simulation
Use `g++` with optimizations (`-O3`) for faster computation:
```bash
g++ -O3 simulation.cpp -o simulation
```

### 3. Run the simulation
You can run the simulation either with an interactive menu or by passing the setup number directly as an argument.

**Option A: Interactive menu**
```bash
./simulation
```

**Option B: Direct execution (e.g., for Setup 3)**
```bash
./simulation 3
```
*This will generate three data files: `oscillation_centre.txt`, `profil_anim.txt`, and `enveloppe_finale.txt`.*

### 4. Visualize the results (Python)
It is highly recommended to use a Python virtual environment to avoid package conflicts.

**Setup the virtual environment:**
```bash
python3 -m venv venv
source venv/bin/activate
pip install numpy matplotlib scipy
```

**Run the visualization script:**
```bash
python plot_results.py
```

## Outputs

The Python script will generate the following visualizations:
1. **Time Series & FFT Spectrum**: Shows the displacement of the plate's center over time and its frequency spectrum (identifying the main resonance frequency).
2. **Chladni Pattern**: A 2D contour plot of the maximum amplitude envelope, highlighting the nodal lines (where sand would accumulate on a real plate).
3. **3D Animation**: A dynamic 3D surface plot showing the plate's deflection over time.

## License

This project is open-source and available under the MIT License.
