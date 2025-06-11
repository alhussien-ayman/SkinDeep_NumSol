import numpy as np
import matplotlib.pyplot as plt

def initial_condition(U):
    U[:, 0] = 0.0

def boundary_condition(U, N):
    for j in range(N + 1):
        U[0, j] = 1.0           
        U[-1, j] = U[-2, j]    

def solve(D=2e-9, r0=0.5, p=0, s_c=8.0e-6, T=30*24*3600, N=1500, M=150): #N=time step, M=space division
    dr = r0 / M
    dt_max = dr**2 / (2 * D)
    dt = min(T / N, dt_max * 0.5)  
    N = int(T / dt)  
    u0 = 1.0
    U = np.zeros((M + 1, N + 1)) 

    initial_condition(U)
    boundary_condition(U, N)

    eps = 1e-8  

    for n in range(N):
        for j in range(1, M):
            u = U[j, n]
           
            ux_fwd = (U[j + 1, n] - U[j, n]) / dr
            ux_bwd = (U[j, n] - U[j - 1, n]) / dr

            
            base_p1 = max(1 - 0.5 * ((U[j + 1, n] + U[j, n])/u0), eps)
            base_p2 = max(1 - 0.5 * ((U[j, n] + U[j - 1, n])/u0), eps)

            phi_p = base_p1**p * ux_fwd
            phi_m = base_p2**p * ux_bwd

            diffusion = (phi_p - phi_m) / dr
            reaction = s_c * u * (1 - u / u0)

            U[j, n + 1] = u + dt * (D * diffusion + reaction)

        U[-1, n + 1] = U[-2, n + 1]

    return U, dr, dt, N


U, dr, dt, N = solve()

r = np.linspace(0, 0.5, U.shape[0])
time_days = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
step_days = [min(int(t * 86400 / dt), N) for t in time_days]

plt.figure(figsize=(9, 5))
for step, day in zip(step_days, time_days):
    plt.plot(r, U[:, step], label=f'{day} days')
plt.xlabel('r (cm)')
plt.ylabel('u(r, t)')
plt.title('epidermal Wound Healing (Explicit finite difference)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
