#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

/* Setups available:
1. Model test (simply supported): Sinusoidal initialization, fixed boundaries (u=0).
2. Free plate (Gaussian deformation): Free boundaries, Gaussian initial deflection.
3. Free plate (Forced regime): Free boundaries, excited by a central vibrator.
*/

/* Global parameters */
double Lx = 1.0, Ly = 1.0;
int nx = 101, ny = 101; 
double dh;
double dh4;
double dt; 
int n_steps = 100000;      
int save_interval = 500;

/* Brass characteristics */
double E = 130e9, h_plate = 0.002, rho = 8500, nu = 0.37, D, alpha;

/* Matrix utility functions */
double** allocate_matrix(int n_x, int n_y) {
    double** mat = (double**)malloc(n_x * sizeof(double*));
    for (int i = 0; i < n_x; i++) {
        mat[i] = (double*)calloc(n_y, sizeof(double)); 
    }
    return mat;
}

void free_matrix(double** mat, int n_x) {
    for (int i = 0; i < n_x; i++) free(mat[i]);
    free(mat);
} 

/* Calculation functions */
void compute_biharmonic(double** u, double** Lu, int nx_phys, int ny_phys, double dh4) {
    for (int i = 2; i <= nx_phys + 1; i++) {
        for (int j = 2; j <= ny_phys + 1; j++) {
            double S1 = u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1];
            double S2 = u[i-1][j-1] + u[i-1][j+1] + u[i+1][j-1] + u[i+1][j+1];
            double S3 = u[i-2][j] + u[i+2][j] + u[i][j-2] + u[i][j+2];
            
            Lu[i][j] = (20.0 * u[i][j] - 8.0 * S1 + 2.0 * S2 + S3) / dh4;
        }
    }
}

void apply_free_boundaries(double** u, int nx_phys, int ny_phys, double nu) {
    int i_start = 2;
    int i_end   = nx_phys + 1;

    // 1. First ghost rows (Zero Bending Moment: Mnn = 0)
    for (int j = i_start; j <= i_end; j++) {
        u[1][j] = 2.0*u[2][j] - u[3][j] - nu*(u[2][j-1] - 2.0*u[2][j] + u[2][j+1]);
        u[i_end + 1][j] = 2.0*u[i_end][j] - u[i_end - 1][j] - nu*(u[i_end][j-1] - 2.0*u[i_end][j] + u[i_end][j+1]);
    }
    for (int i = i_start; i <= i_end; i++) {
        u[i][1] = 2.0*u[i][2] - u[i][3] - nu*(u[i-1][2] - 2.0*u[i][2] + u[i+1][2]);
        u[i][i_end + 1] = 2.0*u[i][i_end] - u[i][i_end - 1] - nu*(u[i-1][i_end] - 2.0*u[i][i_end] + u[i+1][i_end]);
    }

    // 2. Second ghost rows (Zero Shear Force: Vn = 0)
    for (int j = i_start; j <= i_end; j++) {
        double d2w_dy2_i3 = u[3][j-1] - 2.0*u[3][j] + u[3][j+1];
        double d2w_dy2_i1 = u[1][j-1] - 2.0*u[1][j] + u[1][j+1];
        u[0][j] = u[4][j] - 2.0*u[3][j] + 2.0*u[1][j] + (2.0 - nu)*(d2w_dy2_i3 - d2w_dy2_i1);

        double d2w_dy2_i_end_minus_1 = u[i_end - 1][j-1] - 2.0*u[i_end - 1][j] + u[i_end - 1][j+1];
        double d2w_dy2_i_end_plus_1  = u[i_end + 1][j-1] - 2.0*u[i_end + 1][j] + u[i_end + 1][j+1];
        u[i_end + 2][j] = u[i_end - 2][j] - 2.0*u[i_end - 1][j] + 2.0*u[i_end + 1][j] - (2.0 - nu)*(d2w_dy2_i_end_plus_1 - d2w_dy2_i_end_minus_1);
    }
    for (int i = i_start; i <= i_end; i++) {
        double d2w_dx2_j3 = u[i-1][3] - 2.0*u[i][3] + u[i+1][3];
        double d2w_dx2_j1 = u[i-1][1] - 2.0*u[i][1] + u[i+1][1];
        u[i][0] = u[i][4] - 2.0*u[i][3] + 2.0*u[i][1] + (2.0 - nu)*(d2w_dx2_j3 - d2w_dx2_j1);

        double d2w_dx2_j_end_minus_1 = u[i-1][i_end - 1] - 2.0*u[i][i_end - 1] + u[i+1][i_end - 1];
        double d2w_dx2_j_end_plus_1  = u[i-1][i_end + 1] - 2.0*u[i][i_end + 1] + u[i+1][i_end + 1];
        u[i][i_end + 2] = u[i][i_end - 2] - 2.0*u[i][i_end - 1] + 2.0*u[i][i_end + 1] - (2.0 - nu)*(d2w_dx2_j_end_plus_1 - d2w_dx2_j_end_minus_1);
    }

    // 3. Crossed ghost corners (Zero Torsion: Mxy = 0)
    u[1][1] = u[1][3] + u[3][1] - u[3][3];
    u[1][i_end + 1] = u[1][i_end - 1] + u[3][i_end + 1] - u[3][i_end - 1];
    u[i_end + 1][1] = u[i_end - 1][1] + u[i_end + 1][3] - u[i_end - 1][3];
    u[i_end + 1][i_end + 1] = u[i_end - 1][i_end + 1] + u[i_end + 1][i_end - 1] - u[i_end - 1][i_end - 1];
}

void initialize_gaussian(double** u_curr, double** u_prev, int nx, int ny, double dh, double Lx, double Ly) {
    double A = 0.00001;     
    double sigma = 0.05;
    double x0 = 0.5 * Lx;
    double y0 = 0.5 * Ly;

    for (int i = 2; i <= nx + 1; i++) {
        for (int j = 2; j <= ny + 1; j++) {
            double x = (i - 2) * dh;
            double y = (j - 2) * dh;
            double dist = (x - x0)*(x - x0) + (y - y0)*(y - y0);
            u_curr[i][j] = A * exp(-dist / (sigma * sigma));
            u_prev[i][j] = u_curr[i][j];
        }
    }
}

void initialize_sinusoidal(double** u_curr, double** u_prev, int nx, int ny, double dh, double Lx, double Ly) {
    double A = 0.000001; 
    int mode_x = 1;
    int mode_y = 1;

    for (int i = 2; i <= nx + 1; i++) {
        for (int j = 2; j <= ny + 1; j++) {
            double x = (i - 2) * dh;
            double y = (j - 2) * dh;

            u_curr[i][j] = A * sin(mode_x * M_PI * x / Lx) * sin(mode_y * M_PI * y / Ly);
            u_prev[i][j] = u_curr[i][j];
        }
    }
}

double vibrator_force(double t) {
    double amplitude_force = 1.0;
    double frequency = 349.97;  
    return amplitude_force * sin(2.0 * M_PI * frequency * t);
}

int main(int argc, char* argv[]) {
    int setup = 0;

    // Check for command line argument
    if (argc > 1) {
        setup = atoi(argv[1]);
    }

    // Interactive menu if no valid argument is passed
    if (setup < 1 || setup > 3) {
        cout << "========================================" << endl;
        cout << " Vibrating Plate Simulation Setup" << endl;
        cout << "========================================" << endl;
        cout << " 1. Simply Supported (Sinusoidal init)" << endl;
        cout << " 2. Free Plate (Gaussian deformation)" << endl;
        cout << " 3. Free Plate (Forced regime)" << endl;
        cout << "========================================" << endl;
        cout << "Select a setup (1-3): ";
        cin >> setup;
        
        if (setup < 1 || setup > 3) {
            cerr << "Invalid choice. Exiting." << endl;
            return 1;
        }
    }

    cout << "Running setup " << setup << "..." << endl;

    /* Plate dynamics calculation */
    dh = Lx / (nx - 1);
    dh4 = pow(dh, 4);
    D = (E * pow(h_plate, 3)) / (12.0 * (1.0 - nu * nu));
    
    fstream file_anim, file_env, file_time;
    
    /* Von Neumann stability criterion */
    dt = 0.05 * (dh * dh / 4.0) * sqrt((rho * h_plate) / D); 
    alpha = (D * dt * dt) / (rho * h_plate);

    int nx_tot = nx + 4;
    int ny_tot = ny + 4;
    double** u_prev     = allocate_matrix(nx_tot, ny_tot); 
    double** u_curr     = allocate_matrix(nx_tot, ny_tot); 
    double** u_next     = allocate_matrix(nx_tot, ny_tot); 
    double** biharmonic = allocate_matrix(nx_tot, ny_tot); 
    double** u_envelope = allocate_matrix(nx_tot, ny_tot); 
    
    int i_c = 2 + nx / 2;
    int j_c = 2 + ny / 2;

    // Apply initial conditions based on selected setup
    if (setup == 1) {
        initialize_sinusoidal(u_curr, u_prev, nx, ny, dh, Lx, Ly);
    } 
    else if (setup == 2) {
        apply_free_boundaries(u_curr, nx, ny, nu);
        apply_free_boundaries(u_prev, nx, ny, nu);
        initialize_gaussian(u_curr, u_prev, nx, ny, dh, Lx, Ly);
    } 
    else if (setup == 3) {
        apply_free_boundaries(u_curr, nx, ny, nu);
        // u_prev remains at zero, representing the plate at rest before the vibrator turns on
    }

    file_time.open("oscillation_centre.txt", ios::out); 
    file_anim.open("profil_anim.txt", ios::out); 

    /* Time loop */
    for (int n = 0; n < n_steps; n++) { 
        double t = n * dt; 
        
        if (n % save_interval == 0) { 
            for (int i = 2; i <= nx + 1; i++) {
                for (int j = 2; j <= ny + 1; j++) {
                    file_anim << u_curr[i][j] << " "; 
                }
                file_anim << "\n"; 
            }
        }

        if (n % 10 == 0) { 
            file_time << t << " " << u_curr[i_c][j_c] << "\n";
        }

        compute_biharmonic(u_curr, biharmonic, nx, ny, dh4); 

        for (int i = 2; i <= nx + 1; i++) { 
            for (int j = 2; j <= ny + 1; j++) {
                u_next[i][j] = 2.0 * u_curr[i][j] - u_prev[i][j] - alpha * biharmonic[i][j];

                // Apply vibrator force for setup 3
                if (setup == 3 && i == i_c && j == j_c) {
                    double F_t = vibrator_force(t);
                    double force_term = (F_t * dt * dt) / (rho * h_plate * dh * dh);
                    u_next[i][j] += force_term;
                }

                if (n > n_steps / 2) { 
                    if (fabs(u_next[i][j]) > u_envelope[i][j]) {
                        u_envelope[i][j] = fabs(u_next[i][j]);
                    }
                }
            }
        }

        // Apply boundary conditions in the loop
        if (setup == 2 || setup == 3) {
            apply_free_boundaries(u_next, nx, ny, nu);
        }
        // For setup 1 (simply supported), ghost cells remain 0, so we do nothing.

        // Time rotation
        double** temp = u_prev; 
        u_prev = u_curr;
        u_curr = u_next;
        u_next = temp; 

        if (n % (n_steps / 10) == 0) { 
            cout << double(n) / n_steps * 100 << "%" << endl;
        }
    }

    file_env.open("enveloppe_finale.txt", ios::out); 
    for (int i = 2; i <= nx + 1; i++) { 
        double x = (i - 2) * dh; 
        for (int j = 2; j <= ny + 1; j++) {
            double y = (j - 2) * dh; 
            file_env << x << " " << y << " " << u_envelope[i][j] << "\n"; 
        }
        file_env << "\n";
    }

    // Free memory
    free_matrix(u_prev, nx_tot);
    free_matrix(u_curr, nx_tot);
    free_matrix(u_next, nx_tot);
    free_matrix(biharmonic, nx_tot);
    free_matrix(u_envelope, nx_tot); 

    cout << "Simulation finished successfully. Used dt = " << dt << endl; 

    return 0;
}
