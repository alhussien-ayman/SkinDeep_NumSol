
export interface SimulationParams {
  r0: number; // wound half-width (cm)
  D: number; // cell diffusivity (cm²/s)
  sc: number; // logistic source coefficient (1/s)
  p: number; // nonlinear diffusion parameter
  days: number; // simulation duration in days
}

export interface SimulationResult {
  r: number[]; // spatial points
  solutions: number[][]; // array of u-profiles at t_output times
  time: number; // computation time in seconds
  gridPoints: number;
  timeSteps: number | string; // Allow string for 'N/A'
  methodName: string;
  relativeError?: number;
  accuracy?: number;
}

export interface PlotDataItem {
  r: number;
  [key: string]: number | null; // Allow null for missing data points
}

export type DisplayMethod = 'all' | 'fd_explicit' | 'fd_implicit' | 'fd_crank_nicolson' | 'mol' | 'analytical' | 'fem';

export const U0_NORMAL_CELL_DENSITY = 1; // Normal cell density (normalized)

export const LINE_COLORS = ['#e74c3c', '#e67e22', '#f39c12', '#f1c40f', '#2ecc71', '#1abc9c', '#3498db', '#9b59b6', '#8e44ad', '#38598b', '#d35400']; // Added one more color for FEM

export interface ParameterConfig {
  id: keyof SimulationParams;
  label: string;
  defaultValue: number;
  step: number;
  min: number;
  max?: number;
  type?: string; // 'number' or 'text' for scientific notation
}