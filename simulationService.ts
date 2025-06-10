
import { SimulationParams, SimulationResult, U0_NORMAL_CELL_DENSITY } from '../types';
import { N_TERMS_ANALYTICAL } from '../constants';

// Helper to generate time points for output
const getTimeOutput = (totalDays: number): { seconds: number[], labels: string[] } => {
  const t_output_seconds: number[] = [];
  const t_output_labels: string[] = [];
  const intervalDays = Math.max(1, Math.floor(totalDays / 10)); 

  for (let day = intervalDays; day <= totalDays; day += intervalDays) {
    if (t_output_seconds.length < 15) { 
        t_output_seconds.push(day * 24 * 3600);
        t_output_labels.push(`${day}d`);
    } else {
        break;
    }
  }
  if (totalDays > 0 && (totalDays % intervalDays !== 0 || t_output_seconds.length === 0) && !t_output_labels.includes(`${totalDays}d`) && t_output_seconds.length < 15) {
     t_output_seconds.push(totalDays * 24 * 3600);
     t_output_labels.push(`${totalDays}d`);
  }
  if(t_output_seconds.length === 0 && totalDays > 0) {
    t_output_seconds.push(totalDays * 24 * 3600);
    t_output_labels.push(`${totalDays}d`);
  }
  return { seconds: t_output_seconds, labels: t_output_labels };
};

// Tridiagonal Matrix Algorithm (Thomas Algorithm) solver
function solveTridiagonalSystem(a_sub: number[], b_main: number[], c_super: number[], d_rhs: number[]): number[] {
    const N = b_main.length;
    if (N === 0) return [];
    if (N === 1) return [d_rhs[0] / b_main[0]];

    const c_prime = new Array(N).fill(0);
    const d_prime = new Array(N).fill(0);
    const x_sol = new Array(N).fill(0);

    c_prime[0] = c_super[0] / b_main[0];
    d_prime[0] = d_rhs[0] / b_main[0];

    for (let i = 1; i < N; i++) {
        const denominator = b_main[i] - a_sub[i-1] * c_prime[i-1];
        if (Math.abs(denominator) < 1e-12) {
             console.error("TDMA: Division by zero or near zero. Denominator:", denominator, "at index", i);
             return new Array(N).fill(NaN); 
        }
        if (i < N - 1) { 
            c_prime[i] = c_super[i] / denominator;
        }
        d_prime[i] = (d_rhs[i] - a_sub[i-1] * d_prime[i-1]) / denominator;
    }

    x_sol[N-1] = d_prime[N-1];
    for (let i = N - 2; i >= 0; i--) {
        x_sol[i] = d_prime[i] - c_prime[i] * x_sol[i+1];
    }
    return x_sol;
}


export const solveFiniteDifferenceExplicit = (params: SimulationParams): { result: SimulationResult, tOutputLabels: string[] } => {
  const startTime = performance.now();
  const { r0, D, sc, p, days } = params;
  const tf = days * 24 * 3600;
  const u0 = U0_NORMAL_CELL_DENSITY;

  const { seconds: t_output_s, labels: t_output_labels } = getTimeOutput(days);

  const Nr = 101;
  const dr = r0 / (Nr - 1);
  const r_coords = Array.from({ length: Nr }, (_, i) => i * dr);

  // Stability condition for linear diffusion: dt <= dr^2 / (2D_max)
  // For nonlinear D(u) = D_param * (1-u/u0)^p, max(D(u)) = D_param (when u=0, assuming p>=0)
  // Add a safety factor (e.g., 0.1 to 0.5). Here 0.1.
  const dt_cfl_stability_limit = dr * dr / (2 * D); // D is D_param here
  const dt_cfl = 0.1 * dt_cfl_stability_limit; // Safety factor of 0.1

  let Nt = Math.max(1, Math.floor(tf / dt_cfl)); // Ensure Nt is at least 1
  const dt = tf / Nt; // Actual time step used

  let u = new Array(Nr).fill(0); 
  const solutions: number[][] = [];
  let nextOutputIdx = 0;

  for (let n = 0; n < Nt; n++) {
    const u_old = [...u];
    u[0] = u0; // Dirichlet boundary condition at r=0

    // Interior points
    for (let i = 1; i < Nr - 1; i++) {
      const u_rr = (u_old[i + 1] - 2 * u_old[i] + u_old[i - 1]) / (dr * dr);
      const D_coeff = D * Math.pow(Math.max(0, 1 - u_old[i] / u0), p); 
      const diffusion = D_coeff * u_rr;
      const source = sc * u_old[i] * (1 - u_old[i] / u0);
      
      u[i] = u_old[i] + dt * (diffusion + source);
      if (u[i] < 0) u[i] = 0; 
      if (u[i] > u0) u[i] = u0; 
    }

    // Neumann boundary condition at r=r0 (i = Nr-1), du/dr = 0
    // This implies u_Nr = u_Nr-2 for a central difference approximation of du/dr
    // PDE at Nr-1: u_t = D_coeff * (u_Nr - 2*u_Nr-1 + u_Nr-2)/(dr^2) + source
    //           => u_t = D_coeff * 2 * (u_Nr-2 - u_Nr-1)/(dr^2) + source
    const D_coeff_boundary = D * Math.pow(Math.max(0, 1 - u_old[Nr-1] / u0), p);
    const diffusion_at_boundary = D_coeff_boundary * 2 * (u_old[Nr-2] - u_old[Nr-1]) / (dr * dr);
    const source_at_boundary = sc * u_old[Nr-1] * (1 - u_old[Nr-1] / u0);
    u[Nr-1] = u_old[Nr-1] + dt * (diffusion_at_boundary + source_at_boundary);
    
    if (u[Nr-1] < 0) u[Nr-1] = 0;
    if (u[Nr-1] > u0) u[Nr-1] = u0;


    const currentTime = (n + 1) * dt;
    if (nextOutputIdx < t_output_s.length && currentTime >= t_output_s[nextOutputIdx]) {
      solutions.push([...u]);
      nextOutputIdx++;
    }
  }
   if (solutions.length < t_output_s.length && Nt > 0) {
     if (solutions.length === 0 && r_coords.length > 0) solutions.push([...u]); 
     while (solutions.length < t_output_s.length && solutions.length > 0) { 
       solutions.push([...solutions[solutions.length-1]]); 
     }
   }


  const endTime = performance.now();
  return {
    result: {
      r: r_coords,
      solutions: solutions,
      time: (endTime - startTime) / 1000,
      gridPoints: Nr,
      timeSteps: Nt,
      methodName: "Explicit FD (Forward Euler)",
    },
    tOutputLabels: t_output_labels
  };
};

export const solveFiniteDifferenceImplicit = (params: SimulationParams): { result: SimulationResult, tOutputLabels: string[] } => {
    const startTime = performance.now();
    const { r0, D, sc, p, days } = params;
    const tf = days * 24 * 3600;
    const u0 = U0_NORMAL_CELL_DENSITY;

    const { seconds: t_output_s, labels: t_output_labels } = getTimeOutput(days);

    const Nr = 101; 
    const dr = r0 / (Nr - 1);
    const r_coords = Array.from({ length: Nr }, (_, i) => i * dr);

    const Nt = 2000; 
    const dt = tf / Nt;

    let u = new Array(Nr).fill(0); 
    const solutions: number[][] = [];
    let nextOutputIdx = 0;

    const N_internal = Nr - 2; 
    const a_sub = new Array(N_internal - 1);
    const b_main = new Array(N_internal);
    const c_super = new Array(N_internal - 1);
    const d_rhs = new Array(N_internal);

    for (let n = 0; n < Nt; n++) {
        const u_old = [...u];
        
        for (let j = 0; j < N_internal; j++) {
            const i = j + 1; 
            const D_coeff_at_old_u = D * Math.pow(Math.max(0, 1 - u_old[i] / u0), p); // D(u) evaluated at u_old
            const lambda_val = D_coeff_at_old_u * dt / (dr * dr);
            
            b_main[j] = 1 + 2 * lambda_val;
            if (j > 0) a_sub[j-1] = -lambda_val;
            if (j < N_internal - 1) c_super[j] = -lambda_val;
            
            const source_term = sc * u_old[i] * (1 - u_old[i] / u0); // Source term evaluated at u_old
            d_rhs[j] = u_old[i] + dt * source_term;

            // Dirichlet BC u(0,t)=u0 incorporated
            if (j === 0) d_rhs[j] += lambda_val * u0; 
            
            // Neumann BC du/dr(r0,t)=0 incorporated (u_N-1 = u_N-2)
            // For point i = Nr-2 (j = N_internal-1), the equation is -lambda*u_{Nr-3} + (1+lambda)*u_{Nr-2} = rhs
            if (j === N_internal - 1) b_main[j] = 1 + lambda_val; 
        }
        
        const u_internal_new = solveTridiagonalSystem(a_sub, b_main, c_super, d_rhs);
        if (u_internal_new.some(val => isNaN(val))) {
            console.error("Implicit FD: TDMA returned NaN. Marking solutions as problematic.");
            const failed_solution = new Array(Nr).fill(NaN);
            while (solutions.length < t_output_s.length) solutions.push(failed_solution);
            break; 
        }

        u[0] = u0; 
        for (let j = 0; j < N_internal; j++) {
            u[j+1] = u_internal_new[j];
            if (u[j+1] < 0) u[j+1] = 0;
            if (u[j+1] > u0) u[j+1] = u0;
        }
        u[Nr - 1] = u[Nr - 2]; // Sets u_N-1 = u_N-2 based on the newly computed u_N-2

        const currentTime = (n + 1) * dt;
        if (nextOutputIdx < t_output_s.length && currentTime >= t_output_s[nextOutputIdx]) {
            solutions.push([...u]);
            nextOutputIdx++;
        }
    }
    if (solutions.length < t_output_s.length && Nt > 0 && !u.some(isNaN)) {
      if (solutions.length === 0 && r_coords.length > 0) solutions.push([...u]);
      while (solutions.length < t_output_s.length && solutions.length > 0) {
        solutions.push([...solutions[solutions.length-1]]);
      }
    }

    const endTime = performance.now();
    return {
        result: {
            r: r_coords,
            solutions: solutions,
            time: (endTime - startTime) / 1000,
            gridPoints: Nr,
            timeSteps: Nt,
            methodName: "Implicit FD (Backward Euler)",
        },
        tOutputLabels: t_output_labels
    };
};

export const solveFiniteDifferenceCrankNicolson = (params: SimulationParams): { result: SimulationResult, tOutputLabels: string[] } => {
    const startTime = performance.now();
    const { r0, D, sc, p, days } = params;
    const tf = days * 24 * 3600;
    const u0 = U0_NORMAL_CELL_DENSITY;

    const { seconds: t_output_s, labels: t_output_labels } = getTimeOutput(days);

    const Nr = 101;
    const dr = r0 / (Nr - 1);
    const r_coords = Array.from({ length: Nr }, (_, i) => i * dr);

    const Nt = 2000; 
    const dt = tf / Nt;

    let u = new Array(Nr).fill(0);
    const solutions: number[][] = [];
    let nextOutputIdx = 0;

    const N_internal = Nr - 2;
    const a_lhs = new Array(N_internal - 1);
    const b_lhs = new Array(N_internal);
    const c_lhs = new Array(N_internal - 1);
    const d_rhs_cn = new Array(N_internal); 

    for (let n = 0; n < Nt; n++) {
        const u_old = [...u];

        for (let j = 0; j < N_internal; j++) {
            const i = j + 1; 
            const D_coeff_at_old_u = D * Math.pow(Math.max(0, 1 - u_old[i] / u0), p); // D(u) at u_old
            const alpha = D_coeff_at_old_u * dt / (2 * dr * dr); // Factor of 0.5 for Crank-Nicolson

            // LHS matrix coefficients (implicit part)
            b_lhs[j] = 1 + 2 * alpha;
            if (j > 0) a_lhs[j-1] = -alpha;
            if (j < N_internal - 1) c_lhs[j] = -alpha;
            
            // RHS vector (explicit part + source term)
            // Diffusion operator at u_old for point i
            const explicit_diffusion_op_i = D_coeff_at_old_u * (u_old[i+1] - 2*u_old[i] + u_old[i-1]) / (dr*dr);
            const source_term_old_i = sc * u_old[i] * (1 - u_old[i] / u0); // Source at u_old
            
            // RHS: u_old_i + 0.5 * dt * (explicit_diffusion_op_i) + dt * source_term_old_i (source is fully explicit)
            d_rhs_cn[j] = u_old[i] + 0.5 * dt * explicit_diffusion_op_i + dt * source_term_old_i;

            // Incorporate Dirichlet BC u(0,t)=u0 into RHS for the first internal point (j=0, i=1)
            if (j === 0) d_rhs_cn[j] += alpha * u0; // This is for the implicit part u_0^{n+1}
            
            // Neumann BC du/dr(r0,t)=0 (u_N-1 = u_N-2) for the last internal point (j=N_internal-1, i=Nr-2)
            if (j === N_internal - 1) {
              // LHS modification for u_{Nr-1} = u_{Nr-2}
              b_lhs[j] = 1 + alpha; 
              // RHS modification:
              // Explicit diffusion for point Nr-2, considering u_{Nr-1} = u_{Nr-2} in the stencil for u_{Nr-2}
              // (u_{Nr-1} - 2*u_{Nr-2} + u_{Nr-3}) -> (u_{Nr-2} - 2*u_{Nr-2} + u_{Nr-3}) = (u_{Nr-3} - u_{Nr-2})
              const explicit_diffusion_op_boundary = D_coeff_at_old_u * (u_old[Nr-3] - u_old[Nr-2]) / (dr*dr); // D(u_old_Nr-2)
              const source_term_boundary = sc * u_old[Nr-2] * (1 - u_old[Nr-2] / u0); // Source(u_old_Nr-2)
              d_rhs_cn[j] = u_old[Nr-2] + 0.5 * dt * explicit_diffusion_op_boundary + dt * source_term_boundary;
              // The alpha * u_N-1^{n+1} (which is alpha * u_N-2^{n+1}) is handled by b_lhs[j] = 1 + alpha
            }
        }

        const u_internal_new = solveTridiagonalSystem(a_lhs, b_lhs, c_lhs, d_rhs_cn);
        if (u_internal_new.some(val => isNaN(val))) {
            console.error("Crank-Nicolson: TDMA returned NaN. Marking solutions as problematic.");
            const failed_solution = new Array(Nr).fill(NaN);
            while (solutions.length < t_output_s.length) solutions.push(failed_solution);
            break;
        }
        
        u[0] = u0;
        for (let j = 0; j < N_internal; j++) {
            u[j+1] = u_internal_new[j];
            if (u[j+1] < 0) u[j+1] = 0;
            if (u[j+1] > u0) u[j+1] = u0;
        }
        u[Nr - 1] = u[Nr - 2]; // Apply Neumann BC

        const currentTime = (n + 1) * dt;
        if (nextOutputIdx < t_output_s.length && currentTime >= t_output_s[nextOutputIdx]) {
            solutions.push([...u]);
            nextOutputIdx++;
        }
    }
    if (solutions.length < t_output_s.length && Nt > 0 && !u.some(isNaN)) {
      if (solutions.length === 0 && r_coords.length > 0) solutions.push([...u]);
      while (solutions.length < t_output_s.length && solutions.length > 0) {
        solutions.push([...solutions[solutions.length-1]]);
      }
    }
    
    const endTime = performance.now();
    return {
        result: {
            r: r_coords,
            solutions: solutions,
            time: (endTime - startTime) / 1000,
            gridPoints: Nr,
            timeSteps: Nt,
            methodName: "Crank-Nicolson FD",
        },
        tOutputLabels: t_output_labels
    };
};

export const solveMethodOfLines = (params: SimulationParams): { result: SimulationResult, tOutputLabels: string[] } => {
  const startTime = performance.now();
  const { r0, D: D_param, sc, p, days } = params; 
  const tf = days * 24 * 3600;
  const u0 = U0_NORMAL_CELL_DENSITY;

  const { seconds: t_output_s, labels: t_output_labels } = getTimeOutput(days);

  const Nr = 101; 
  const dr = r0 / (Nr - 1);
  const r_coords = Array.from({ length: Nr }, (_, i) => i * dr);

  const Nt_mol = 500000; // Fixed, potentially large number of steps for MOL with explicit Euler in time
  const dt_mol = tf / Nt_mol;

  let u = new Array(Nr).fill(0); 
  u[0] = u0; // Initialize Dirichlet boundary condition u(0,t)=u0 BEFORE the loop

  const solutions: number[][] = [];
  let nextOutputIdx = 0;
  let simulationHasFailed = false;

  for (let n = 0; n < Nt_mol; n++) {
    const u_current_iter = [...u]; 
    const dudt = new Array(Nr).fill(0);
    // dudt[0] remains 0 as u[0] is fixed (Dirichlet condition)

    // Interior points: 1 to Nr-2 (inclusive)
    for (let i = 1; i < Nr - 1; i++) {
      const u_mid_plus = (u_current_iter[i+1] + u_current_iter[i]) / 2;
      const u_mid_minus = (u_current_iter[i] + u_current_iter[i-1]) / 2;

      const D_eff_plus = D_param * Math.pow(Math.max(1e-9, 1 - u_mid_plus / u0), p);
      const D_eff_minus = D_param * Math.pow(Math.max(1e-9, 1 - u_mid_minus / u0), p);

      const flux_plus = D_eff_plus * (u_current_iter[i+1] - u_current_iter[i]) / dr;
      const flux_minus = D_eff_minus * (u_current_iter[i] - u_current_iter[i-1]) / dr;
      
      const diffusion_term = (flux_plus - flux_minus) / dr;
      const source_term = sc * u_current_iter[i] * (1 - u_current_iter[i] / u0);
      
      dudt[i] = diffusion_term + source_term;

      if (!isFinite(dudt[i])) {
        console.error(`MOL: dudt[${i}] became non-finite at t=${(n*dt_mol).toExponential(2)}s, u=${u_current_iter[i].toExponential(2)}. Aborting.`);
        simulationHasFailed = true;
        break;
      }
    }
    if (simulationHasFailed) break;

    // Boundary condition at r=r0 (i = Nr-1): du/dr = 0 implies flux is 0 across this boundary.
    // Flux_plus (outward) is zero. We only consider flux_minus (inward).
    const u_mid_minus_boundary = (u_current_iter[Nr-1] + u_current_iter[Nr-2]) / 2;
    const D_eff_minus_boundary = D_param * Math.pow(Math.max(1e-9, 1 - u_mid_minus_boundary / u0), p);
    const flux_minus_boundary = D_eff_minus_boundary * (u_current_iter[Nr-1] - u_current_iter[Nr-2]) / dr; 
    
    const diffusion_boundary_term = (0 - flux_minus_boundary) / dr; 
    const source_boundary_term = sc * u_current_iter[Nr-1] * (1 - u_current_iter[Nr-1] / u0);
    dudt[Nr-1] = diffusion_boundary_term + source_boundary_term;
    
    if (!isFinite(dudt[Nr-1])) {
        console.error(`MOL: dudt[${Nr-1}] (boundary) became non-finite at t=${(n*dt_mol).toExponential(2)}s. Aborting.`);
        simulationHasFailed = true;
        break;
    }

    // Update u using Forward Euler for points i=1 to Nr-1. u[0] is not changed.
    for (let i = 1; i < Nr; i++) { 
      u[i] = u_current_iter[i] + dt_mol * dudt[i];
      if (u[i] < 0) u[i] = 0;
      if (u[i] > u0) u[i] = u0;
    }

    const currentTime = (n + 1) * dt_mol;
    if (nextOutputIdx < t_output_s.length && currentTime >= t_output_s[nextOutputIdx]) {
      solutions.push([...u]);
      nextOutputIdx++;
    }
  }

  if (simulationHasFailed) {
    const nanSolution = new Array(Nr).fill(NaN);
    while (solutions.length < t_output_s.length) {
        solutions.push(nanSolution);
    }
  } else if (solutions.length < t_output_s.length && Nt_mol > 0) {
    if (solutions.length === 0 && r_coords.length > 0) solutions.push([...u]); 
    while (solutions.length < t_output_s.length && solutions.length > 0) { 
      solutions.push([...solutions[solutions.length-1]]);
    }
  }

  const endTime = performance.now();
  return {
    result: {
      r: r_coords,
      solutions: solutions,
      time: (endTime - startTime) / 1000,
      gridPoints: Nr,
      timeSteps: Nt_mol, 
      methodName: "Method of Lines",
    },
    tOutputLabels: t_output_labels
  };
};


export const solveFiniteElementMethod = (params: SimulationParams): { result: SimulationResult, tOutputLabels: string[] } => {
  const startTime = performance.now();
  const { r0, D: D_param, sc, p, days } = params;
  const tf = days * 24 * 3600;
  const u0 = U0_NORMAL_CELL_DENSITY;
  const eps = 1e-9; // Small epsilon to prevent Math.pow(0, p) issues if p < 0 or u_mid = u0

  const { seconds: t_output_s, labels: t_output_labels } = getTimeOutput(days);

  const Nr = 101; // Number of nodes
  const dr = r0 / (Nr - 1); // Element length
  const r_coords = Array.from({ length: Nr }, (_, i) => i * dr);

  // Stability for explicit Euler in time with FEM can be complex.
  // Using a dt similar to MOL as a starting point.
  const Nt_fem = 500000; 
  const dt_fem = tf / Nt_fem;

  let u = new Array(Nr).fill(0);
  u[0] = u0; // Initialize Dirichlet boundary condition u(0,t)=u0

  const solutions: number[][] = [];
  let nextOutputIdx = 0;
  let simulationHasFailed = false;

  for (let n = 0; n < Nt_fem; n++) {
    const u_current_iter = [...u];
    const dudt = new Array(Nr).fill(0);
    // dudt[0] remains 0 as u[0] is fixed (Dirichlet condition)

    // Interior nodes: 1 to Nr-2 (inclusive)
    // Node i is associated with basis function phi_i.
    // It interacts with element (i-1, i) and element (i, i+1).
    for (let i = 1; i < Nr - 1; i++) {
      // Flux from element (i, i+1) into node i (or out of node i) - "flux_plus" for MOL
      const u_avg_elem_plus = (u_current_iter[i+1] + u_current_iter[i]) / 2;
      const D_eff_elem_plus = D_param * Math.pow(Math.max(eps, 1 - u_avg_elem_plus / u0), p);
      const flux_elem_plus = D_eff_elem_plus * (u_current_iter[i+1] - u_current_iter[i]) / dr; // Flux over element (i,i+1)

      // Flux from element (i-1, i) into node i (or out of node i-1) - "flux_minus" for MOL
      const u_avg_elem_minus = (u_current_iter[i] + u_current_iter[i-1]) / 2;
      const D_eff_elem_minus = D_param * Math.pow(Math.max(eps, 1 - u_avg_elem_minus / u0), p);
      const flux_elem_minus = D_eff_elem_minus * (u_current_iter[i] - u_current_iter[i-1]) / dr; // Flux over element (i-1,i)
      
      // Diffusion term for node i: (flux_elem_plus - flux_elem_minus) / dr (due to lumped mass M_ii = dr)
      const diffusion_term_at_node_i = (flux_elem_plus - flux_elem_minus) / dr;

      // Source term for node i: Integral(phi_i * S(u_avg_elem))dr / M_ii
      // M_ii (lumped) = dr. Integral(phi_i over elem_minus) = dr/2. Integral(phi_i over elem_plus) = dr/2.
      const S_val_elem_plus = sc * u_avg_elem_plus * (1 - u_avg_elem_plus / u0);
      const S_val_elem_minus = sc * u_avg_elem_minus * (1 - u_avg_elem_minus / u0);
      const source_term_at_node_i = (S_val_elem_minus + S_val_elem_plus) / 2;
      
      dudt[i] = diffusion_term_at_node_i + source_term_at_node_i;

      if (!isFinite(dudt[i])) {
        console.error(`FEM: dudt[${i}] became non-finite. Aborting.`);
        simulationHasFailed = true;
        break;
      }
    }
    if (simulationHasFailed) break;

    // Boundary node Nr-1 (Neumann du/dr=0)
    // Basis function phi_{Nr-1} only involves element (Nr-2, Nr-1).
    // Flux from "right" (element Nr-1, Nr) is zero.
    const u_avg_elem_boundary = (u_current_iter[Nr-1] + u_current_iter[Nr-2]) / 2;
    const D_eff_elem_boundary = D_param * Math.pow(Math.max(eps, 1 - u_avg_elem_boundary / u0), p);
    const flux_elem_boundary = D_eff_elem_boundary * (u_current_iter[Nr-1] - u_current_iter[Nr-2]) / dr;
    
    // Diffusion for node Nr-1: (0 - flux_elem_boundary) / dr_boundary_mass_factor
    // If M_{Nr-1,Nr-1} (lumped) is dr/2 for boundary, then factor is dr/2.
    // For simplicity, if we treat node Nr-1 as having lumped mass dr (like interior), 
    // then this approximation means the basis function effectively extends.
    // Let's assume lumped mass M_ii = dr for all i=1..Nr-1 for simplicity,
    // consistent with the (flux_plus - flux_minus)/dr scaling.
    const diffusion_boundary_node = (0 - flux_elem_boundary) / dr;

    // Source for node Nr-1: Integral(phi_{Nr-1} * S(u_avg_elem_boundary)) / M_{Nr-1,Nr-1}
    // Integral(phi_{Nr-1} over elem_boundary) = dr/2.
    // M_{Nr-1,Nr-1} (lumped) = dr implies source_term = S_val_elem_boundary * (dr/2) / dr = S_val_elem_boundary / 2.
    // If M_{Nr-1,Nr-1} (lumped) is dr/2 (more standard for boundary element), then source_term = S_val_elem_boundary.
    // Let's be consistent with interior nodes: (S_left + S_right)/2. Here S_right is effectively zero or not well-defined.
    // So, use S_val_elem_boundary (from the only contributing element).
    // This corresponds to M_lumped_{N-1,N-1} = dr/2.
    // So, dudt[Nr-1] needs to be multiplied by 2 if other terms assumed M_lumped = dr.
    // Alternative: Use the average source logic: (S_left + S_right_conceptual_zero)/2
    const S_val_elem_boundary = sc * u_avg_elem_boundary * (1 - u_avg_elem_boundary / u0);
    // Using same structure as interior for source: (S_left + S_right)/2, where S_right component isn't there.
    // So, source is S_val_elem_boundary / 2 if we are strictly following the interior formula.
    // However, a more direct derivation for node Nr-1 (integral of phi_{Nr-1} S_elem / M_ii) with M_ii = dr/2 gives S_elem.
    // And for diffusion: ( (0 - flux_boundary)/dr_element ) / (dr_mass/dr_element).
    // For now, let's keep the logic parallel to MOL for the boundary:
    const source_boundary_node = S_val_elem_boundary; // Assuming M_ii = dr/2 for this node for source integration
                                                    // and diffusion flux scaling (0 - flux_elem_boundary)/(dr/2)
                                                    // To match (flux_plus-flux_minus)/dr, this needs factor of 2.
                                                    // Let's use simpler parallel from MOL for boundary:
    dudt[Nr-1] = diffusion_boundary_node + source_boundary_node; // Matches MOL structure.

    if (!isFinite(dudt[Nr-1])) {
      console.error(`FEM: dudt[${Nr-1}] (boundary) became non-finite. Aborting.`);
      simulationHasFailed = true;
      break;
    }

    // Update u using Forward Euler for points i=1 to Nr-1. u[0] is not changed.
    for (let i = 1; i < Nr; i++) {
      u[i] = u_current_iter[i] + dt_fem * dudt[i];
      if (u[i] < 0) u[i] = 0;
      if (u[i] > u0) u[i] = u0;
    }

    const currentTime = (n + 1) * dt_fem;
    if (nextOutputIdx < t_output_s.length && currentTime >= t_output_s[nextOutputIdx]) {
      solutions.push([...u]);
      nextOutputIdx++;
    }
  }

  if (simulationHasFailed) {
    const nanSolution = new Array(Nr).fill(NaN);
    while (solutions.length < t_output_s.length) {
      solutions.push(nanSolution);
    }
  } else if (solutions.length < t_output_s.length && Nt_fem > 0) {
    if (solutions.length === 0 && r_coords.length > 0) solutions.push([...u]);
    while (solutions.length < t_output_s.length && solutions.length > 0) {
      solutions.push([...solutions[solutions.length - 1]]);
    }
  }

  const endTime = performance.now();
  return {
    result: {
      r: r_coords,
      solutions: solutions,
      time: (endTime - startTime) / 1000,
      gridPoints: Nr,
      timeSteps: Nt_fem,
      methodName: "Finite Element Method",
    },
    tOutputLabels: t_output_labels
  };
};


export const solveAnalytical = (params: SimulationParams): { result: SimulationResult, tOutputLabels: string[] } => {
  const startTime = performance.now();
  const { r0, D, sc, p, days } = params;
  const u0 = U0_NORMAL_CELL_DENSITY;

  if (p !== 0 || sc !== 0) {
    return {
      result: {
        r: [],
        solutions: [],
        time: (performance.now() - startTime) / 1000,
        gridPoints: 0,
        timeSteps: "N/A",
        methodName: "Analytical (N/A: p≠0 or sc≠0)",
      },
      tOutputLabels: []
    };
  }

  const tf = days * 24 * 3600;
  const { seconds: t_output_s, labels: t_output_labels } = getTimeOutput(days);
  
  const Nr = 101; 
  const dr = r0 / (Nr - 1);
  const r_coords = Array.from({ length: Nr }, (_, i) => i * dr);
  const solutions: number[][] = [];

  for (const t of t_output_s) {
    const u_profile_t = new Array(Nr).fill(0);
    for (let i = 0; i < Nr; i++) {
      const r_val = r_coords[i];
      if (r_val === 0) { // Boundary condition u(0,t) = u0
        u_profile_t[i] = u0;
        continue;
      }
      // Solution for v(r,t) = u(r,t) - u0, with v(0,t)=0, dv/dr(r0,t)=0, v(r,0)=-u0
      let v_rt = 0;
      for (let n_term = 1; n_term <= N_TERMS_ANALYTICAL; n_term++) {
        const lambda_n = (2 * n_term - 1) * Math.PI / (2 * r0);
        const Bn = -4 * u0 / ((2 * n_term - 1) * Math.PI);
        const term = Bn * Math.sin(lambda_n * r_val) * Math.exp(-D * lambda_n * lambda_n * t);
        v_rt += term;
      }
      let u_val = u0 + v_rt; // u(r,t) = v(r,t) + u0
      if (u_val < 0) u_val = 0;
      if (u_val > u0) u_val = u0;
      u_profile_t[i] = u_val;
    }
    solutions.push(u_profile_t);
  }
  
   if (solutions.length < t_output_s.length && t_output_s.length > 0) {
     if (solutions.length === 0 && r_coords.length > 0) { 
        const last_u_profile = new Array(Nr).fill(0); 
        console.warn("Analytical solution padding: May indicate issue if t_output_s was non-empty.")
        solutions.push(last_u_profile); 
     }
     while (solutions.length < t_output_s.length && solutions.length > 0) {
       solutions.push([...solutions[solutions.length-1]]);
     }
   }

  const endTime = performance.now();
  return {
    result: {
      r: r_coords,
      solutions: solutions,
      time: (endTime - startTime) / 1000,
      gridPoints: Nr,
      timeSteps: "N/A", 
      methodName: "Analytical Solution",
    },
    tOutputLabels: t_output_labels
  };
};