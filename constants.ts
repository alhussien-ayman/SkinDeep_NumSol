
import { SimulationParams, ParameterConfig } from './types';

export const DEFAULT_SIMULATION_PARAMS: SimulationParams = {
  r0: 0.5,
  D: 2.0e-9,
  sc: 8.0e-6,
  p: 0, // Changed from 5 back to 0
  days: 30,
};

export const PARAMETER_CONFIGS: ParameterConfig[] = [
  { id: 'r0', label: 'Wound Half-Width (cm)', defaultValue: DEFAULT_SIMULATION_PARAMS.r0, step: 0.1, min: 0.1, max: 2.0 },
  { id: 'D', label: 'Diffusivity (cm²/s)', defaultValue: DEFAULT_SIMULATION_PARAMS.D, step: 1e-10, min: 1e-10, type: 'text'}, // Use text for e-notation
  { id: 'sc', label: 'Source Coeff. (1/s)', defaultValue: DEFAULT_SIMULATION_PARAMS.sc, step: 1e-7, min: 0, type: 'text'}, // Use text for e-notation, min changed to 0
  { id: 'p', label: 'Nonlinearity Param.', defaultValue: 0, step: 0.1, min: 0, max: 10 }, // Changed defaultValue from DEFAULT_SIMULATION_PARAMS.p to 0 directly. Increased max p to 10
  { id: 'days', label: 'Simulation Days', defaultValue: DEFAULT_SIMULATION_PARAMS.days, step: 1, min: 1, max: 90 },
];

export const N_TERMS_ANALYTICAL = 100; // Number of terms for the analytical series solution