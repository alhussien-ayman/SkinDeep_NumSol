
import React, { useState, useEffect, useCallback, useMemo } from 'react';
import ParameterInput from './components/ParameterInput';
import SimulationPlot from './components/SimulationPlot.tsx'; // Changed this line
import { SimulationParams, SimulationResult, PlotDataItem, DisplayMethod, LINE_COLORS, ParameterConfig } from './types';
import { DEFAULT_SIMULATION_PARAMS, PARAMETER_CONFIGS } from './constants';
import { 
  solveFiniteDifferenceExplicit, 
  solveFiniteDifferenceImplicit,
  solveFiniteDifferenceCrankNicolson,
  solveMethodOfLines,
  solveAnalytical,
  solveFiniteElementMethod
} from './services/simulationService';

interface TransformedPlotData {
  plotData: PlotDataItem[];
  timeLabels: string[];
}

const transformToRechartsData = (r_coords: number[] | undefined, solutions: number[][] | undefined, timeLabels: string[] | undefined): PlotDataItem[] => {
  if (!r_coords || r_coords.length === 0 || !solutions || solutions.length === 0 || !timeLabels || timeLabels.length === 0) {
     return [];
  }
  // Check if every solution profile is invalid (all NaNs) or if lengths mismatch
  if (solutions.some(s => !s || s.length !== r_coords.length)) {
    console.warn("Mismatch between r_coords length and solutions profiles length, or invalid solution profile found. Returning empty data for plotting.");
    return r_coords.map(r_val => ({ r: parseFloat(r_val.toFixed(3)) })); // Return r-coordinates only
  }
  // If all solutions are entirely NaN, effectively no data to plot time-series for.
  if (solutions.every(s => s.every(val => isNaN(val)))) {
    console.warn("All solution profiles are NaN. Returning r-coordinates only for plotting structure.");
    return r_coords.map(r_val => ({ r: parseFloat(r_val.toFixed(3)) }));
  }

  return r_coords.map((r_val, r_idx) => {
    const point: PlotDataItem = { r: parseFloat(r_val.toFixed(3)) };
    solutions.forEach((u_profile, t_idx) => {
      if (t_idx < timeLabels.length) { // Ensure timeLabel exists
        if (u_profile && u_profile[r_idx] !== undefined && !isNaN(u_profile[r_idx])) {
          point[timeLabels[t_idx]] = parseFloat(u_profile[r_idx].toFixed(4));
        } else {
          point[timeLabels[t_idx]] = null; // Set to null for missing, NaN, or undefined data
        }
      }
    });
    return point;
  });
};


const App: React.FC = () => {
  const [params, setParams] = useState<SimulationParams>(DEFAULT_SIMULATION_PARAMS);
  const [paramInputs, setParamInputs] = useState<{[key: string]: string}>(() => {
    const initialInputs: {[key: string]: string} = {};
    PARAMETER_CONFIGS.forEach(p => {
      initialInputs[p.id] = String(p.defaultValue);
    });
    return initialInputs;
  });
  const [paramErrors, setParamErrors] = useState<{[key: string]: string | null}>({});

  // Store raw results and labels
  const [fdExplicitRawResult, setFdExplicitRawResult] = useState<SimulationResult | null>(null);
  const [fdExplicitTOutputLabels, setFdExplicitTOutputLabels] = useState<string[] | null>(null);
  const [fdImplicitRawResult, setFdImplicitRawResult] = useState<SimulationResult | null>(null);
  const [fdImplicitTOutputLabels, setFdImplicitTOutputLabels] = useState<string[] | null>(null);
  const [fdCNRawResult, setFdCNRawResult] = useState<SimulationResult | null>(null);
  const [fdCNTOutputLabels, setFdCNTOutputLabels] = useState<string[] | null>(null);
  const [molRawResult, setMolRawResult] = useState<SimulationResult | null>(null);
  const [molTOutputLabels, setMolTOutputLabels] = useState<string[] | null>(null);
  const [femRawResult, setFemRawResult] = useState<SimulationResult | null>(null);
  const [femTOutputLabels, setFemTOutputLabels] = useState<string[] | null>(null);
  const [analyticalRawResultState, setAnalyticalRawResultState] = useState<SimulationResult | null>(null);
  const [analyticalTOutputLabels, setAnalyticalTOutputLabels] = useState<string[] | null>(null);
  
  const [comparisonPlotData, setComparisonPlotData] = useState<PlotDataItem[] | null>(null);
  const [isLoading, setIsLoading] = useState<boolean>(true);
  const [displayMethod, setDisplayMethod] = useState<DisplayMethod>('all');

  const handleParamChange = useCallback((id: keyof SimulationParams, value: string) => {
    setParamInputs(prev => ({ ...prev, [id]: value }));
    setParamErrors(prev => ({ ...prev, [id]: null }));
  }, []);

  const runSimulation = useCallback(() => {
    setIsLoading(true);
    setFdExplicitRawResult(null); setFdImplicitRawResult(null); setFdCNRawResult(null);
    setMolRawResult(null); setFemRawResult(null); setAnalyticalRawResultState(null);
    setFdExplicitTOutputLabels(null); setFdImplicitTOutputLabels(null); setFdCNTOutputLabels(null);
    setMolTOutputLabels(null); setFemTOutputLabels(null); setAnalyticalTOutputLabels(null);
    setComparisonPlotData(null);
    
    const currentErrors: {[key: string]: string | null} = {};
    let validParamsOverall = true;
    const currentSimParams: SimulationParams = { ...DEFAULT_SIMULATION_PARAMS }; 

    (Object.keys(paramInputs) as Array<keyof SimulationParams>).forEach(key => {
        const valueStr = paramInputs[key].trim();
        const config = PARAMETER_CONFIGS.find(p => p.id === key);
        let fieldError: string | null = null;
        
        if (valueStr === "") {
            fieldError = "Value cannot be empty.";
        } else {
            const valueNum = parseFloat(valueStr);
            if (isNaN(valueNum)) {
                fieldError = "Not a valid number.";
            } else if (config) {
                if (config.min !== undefined && valueNum < config.min) {
                    fieldError = `Min value is ${config.min}.`;
                } else if (config.max !== undefined && valueNum > config.max) {
                    fieldError = `Max value is ${config.max}.`;
                } else {
                    currentSimParams[key] = valueNum;
                }
            } else {
                 currentSimParams[key] = valueNum; 
            }
        }
        
        if (fieldError) {
            currentErrors[key] = fieldError;
            validParamsOverall = false;
        } else {
            currentErrors[key] = null;
        }
    });

    setParamErrors(currentErrors);

    if (!validParamsOverall) {
        setIsLoading(false);
        return; 
    }
    setParams(currentSimParams);

    setTimeout(() => { 
      try {
        let { result: fdExpData, tOutputLabels: fdExpLabels } = solveFiniteDifferenceExplicit(currentSimParams);
        let { result: fdImpData, tOutputLabels: fdImpLabels } = solveFiniteDifferenceImplicit(currentSimParams);
        let { result: fdCNData, tOutputLabels: fdCNLabels } = solveFiniteDifferenceCrankNicolson(currentSimParams);
        let { result: molData, tOutputLabels: molLabels } = solveMethodOfLines(currentSimParams);
        let { result: femData, tOutputLabels: femLabels } = solveFiniteElementMethod(currentSimParams);
        let { result: anData, tOutputLabels: anLabels } = solveAnalytical(currentSimParams);

        if (anData && !anData.methodName.includes("(N/A") && anData.solutions.length > 0 && anData.r.length > 0) {
          const analyticalFinalProfile = anData.solutions[anData.solutions.length - 1];

          const calculateMetrics = (numericalResult: SimulationResult): { relativeError?: number; accuracy?: number } => {
            if (numericalResult.solutions.length > 0 && numericalResult.r.length === analyticalFinalProfile.length) {
              const numericalFinalProfile = numericalResult.solutions[numericalResult.solutions.length - 1];
              let diffSqSum = 0;
              let anSqSum = 0;
              let validPoints = 0;
              let numericalHasNaNs = false;

              for (let i = 0; i < analyticalFinalProfile.length; i++) {
                if (numericalFinalProfile[i] === undefined || isNaN(numericalFinalProfile[i])) {
                  numericalHasNaNs = true; 
                  continue; 
                }
                if (analyticalFinalProfile[i] === undefined || isNaN(analyticalFinalProfile[i])) {
                  continue; 
                }
                const diff = numericalFinalProfile[i] - analyticalFinalProfile[i];
                diffSqSum += diff * diff;
                anSqSum += analyticalFinalProfile[i] * analyticalFinalProfile[i];
                validPoints++;
              }
              
              if (numericalHasNaNs && validPoints < analyticalFinalProfile.length / 2) { 
                return { relativeError: Infinity, accuracy: 0 };
              }
              if (validPoints === 0) return { relativeError: Infinity, accuracy: 0 };

              const l2NormError = Math.sqrt(diffSqSum);
              const l2NormAnalytical = Math.sqrt(anSqSum);

              const relativeError = l2NormAnalytical > 1e-9 ? l2NormError / l2NormAnalytical : (l2NormError > 1e-9 ? Infinity : 0);
              const accuracy = Math.max(0, 1 - relativeError);
              return { relativeError, accuracy };
            }
            return { relativeError: Infinity, accuracy: 0 }; 
          };
          
          if (fdExpData) fdExpData = { ...fdExpData, ...calculateMetrics(fdExpData) };
          if (fdImpData) fdImpData = { ...fdImpData, ...calculateMetrics(fdImpData) };
          if (fdCNData) fdCNData = { ...fdCNData, ...calculateMetrics(fdCNData) };
          if (molData) molData = { ...molData, ...calculateMetrics(molData) };
          if (femData) femData = { ...femData, ...calculateMetrics(femData) };
        }
        
        setFdExplicitRawResult(fdExpData); setFdExplicitTOutputLabels(fdExpLabels);
        setFdImplicitRawResult(fdImpData); setFdImplicitTOutputLabels(fdImpLabels);
        setFdCNRawResult(fdCNData); setFdCNTOutputLabels(fdCNLabels);
        setMolRawResult(molData); setMolTOutputLabels(molLabels);
        setFemRawResult(femData); setFemTOutputLabels(femLabels);
        setAnalyticalRawResultState(anData); setAnalyticalTOutputLabels(anLabels);

        const comparisonDataItems: PlotDataItem[] = [];
        const isAnalyticalApplicableForComparison = anData && !anData.methodName.includes("(N/A") && anData.r.length > 0;
        
        const allResultsForComparison = [
            fdExpData ? { rawResult: fdExpData, name: fdExpData.methodName } : null,
            fdImpData ? { rawResult: fdImpData, name: fdImpData.methodName } : null,
            fdCNData ? { rawResult: fdCNData, name: 'Crank-Nicolson FD' } : null,
            molData ? { rawResult: molData, name: 'MOL' } : null,
            femData ? { rawResult: femData, name: 'FEM' } : null,
        ].filter(item => item !== null) as { rawResult: SimulationResult; name: string }[];

        if (isAnalyticalApplicableForComparison && anData) {
            allResultsForComparison.push({ rawResult: anData, name: 'Analytical' });
        }
        
        const refRCoords = fdExpData?.r || molData?.r || femData?.r || anData?.r || [];

        if (refRCoords.length > 0) {
            refRCoords.forEach((r_val, idx) => {
                const point: PlotDataItem = { r: parseFloat(r_val.toFixed(3)) };
                allResultsForComparison.forEach(resItem => {
                    const res = resItem.rawResult;
                    if (res && res.solutions && res.solutions.length > 0 && idx < res.r.length) {
                        const lastSolutionProfile = res.solutions[res.solutions.length - 1];
                        if (lastSolutionProfile && idx < lastSolutionProfile.length && !isNaN(lastSolutionProfile[idx])) {
                            const val = lastSolutionProfile[idx];
                            point[resItem.name] = parseFloat(val.toFixed(4));
                        } else {
                             point[resItem.name] = null;
                        }
                    } else {
                        point[resItem.name] = null;
                    }
                });
                comparisonDataItems.push(point);
            });
            setComparisonPlotData(comparisonDataItems);
        }

      } catch (error) {
        console.error("Simulation error:", error);
        alert("An error occurred during simulation. Check console for details.");
      } finally {
        setIsLoading(false);
      }
    }, 50); 
  }, [paramInputs]); 

  useEffect(() => {
    runSimulation(); 
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []); 

  const fdExplicitTransformed = useMemo(() => fdExplicitRawResult && fdExplicitTOutputLabels ? { plotData: transformToRechartsData(fdExplicitRawResult.r, fdExplicitRawResult.solutions, fdExplicitTOutputLabels), timeLabels: fdExplicitTOutputLabels } : null, [fdExplicitRawResult, fdExplicitTOutputLabels]);
  const fdImplicitTransformed = useMemo(() => fdImplicitRawResult && fdImplicitTOutputLabels ? { plotData: transformToRechartsData(fdImplicitRawResult.r, fdImplicitRawResult.solutions, fdImplicitTOutputLabels), timeLabels: fdImplicitTOutputLabels } : null, [fdImplicitRawResult, fdImplicitTOutputLabels]);
  const fdCNTransformed = useMemo(() => fdCNRawResult && fdCNTOutputLabels ? { plotData: transformToRechartsData(fdCNRawResult.r, fdCNRawResult.solutions, fdCNTOutputLabels), timeLabels: fdCNTOutputLabels } : null, [fdCNRawResult, fdCNTOutputLabels]);
  const molTransformed = useMemo(() => molRawResult && molTOutputLabels ? { plotData: transformToRechartsData(molRawResult.r, molRawResult.solutions, molTOutputLabels), timeLabels: molTOutputLabels } : null, [molRawResult, molTOutputLabels]);
  const femTransformed = useMemo(() => femRawResult && femTOutputLabels ? { plotData: transformToRechartsData(femRawResult.r, femRawResult.solutions, femTOutputLabels), timeLabels: femTOutputLabels } : null, [femRawResult, femTOutputLabels]);
  const analyticalTransformed = useMemo(() => analyticalRawResultState && analyticalTOutputLabels && !analyticalRawResultState.methodName.includes("(N/A") ? { plotData: transformToRechartsData(analyticalRawResultState.r, analyticalRawResultState.solutions, analyticalTOutputLabels), timeLabels: analyticalTOutputLabels } : null, [analyticalRawResultState, analyticalTOutputLabels]);

  const isAnalyticalApplicable = analyticalRawResultState && !analyticalRawResultState.methodName.includes("(N/A");

  const performanceResultsText = useMemo(() => {
    return [fdExplicitRawResult, fdImplicitRawResult, fdCNRawResult, molRawResult, femRawResult, analyticalRawResultState]
      .filter(res => res !== null)
      .map(res => {
        const name = res!.methodName.padEnd(28); // Adjusted padding for longer name
        const grid = String(
          res!.methodName.startsWith("Analytical") 
            ? "N/A" 
            : res!.gridPoints
        ).padStart(11);
        const steps = String(res!.timeSteps).padStart(10);
        const cpuTime = res!.time.toFixed(3).padStart(11);
        
        let relErrorStr = "N/A".padStart(12);
        let accStr = "N/A".padStart(12);

        if (res!.methodName !== "Analytical Solution" && !res!.methodName.includes("(N/A")) { 
            if (res!.relativeError !== undefined ) { 
                relErrorStr = (res!.relativeError === Infinity ? "Inf" : (isNaN(res!.relativeError) ? "NaN" : res!.relativeError.toExponential(2))).padStart(12);
            }
            if (res!.accuracy !== undefined ) { 
                accStr = (isNaN(res!.accuracy) ? "NaN" : (res!.accuracy * 100).toFixed(1) + "%").padStart(12);
            }
        }
        return `${name} | ${grid} | ${steps} | ${cpuTime} | ${relErrorStr} | ${accStr}`;
      }).join('\n');
  }, [fdExplicitRawResult, fdImplicitRawResult, fdCNRawResult, molRawResult, femRawResult, analyticalRawResultState]);
  
  const comparisonLegendLabels = useMemo(() => {
    const labels = [];
    if (fdExplicitRawResult) labels.push(fdExplicitRawResult.methodName);
    if (fdImplicitRawResult) labels.push(fdImplicitRawResult.methodName);
    if (fdCNRawResult) labels.push('Crank-Nicolson FD');
    if (molRawResult) labels.push('MOL');
    if (femRawResult) labels.push('FEM');
    if (isAnalyticalApplicable && analyticalRawResultState && analyticalRawResultState.r.length > 0) {
        labels.push('Analytical');
    }
    return labels;
  }, [fdExplicitRawResult, fdImplicitRawResult, fdCNRawResult, molRawResult, femRawResult, analyticalRawResultState, isAnalyticalApplicable]);


  return (
    <div className="min-h-screen p-4 md:p-8 text-gray-800">
      <div className="max-w-7xl mx-auto bg-white/90 backdrop-blur-md shadow-2xl rounded-3xl p-6 md:p-10">
        <header className="text-center mb-8">
          <h1 className="text-4xl md:text-5xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-indigo-600 to-purple-600 pb-2">
            Epidermal Wound Healing PDE Solver
          </h1>
          <p className="text-lg text-gray-600">
            Comparing Numerical Methods for PDE Simulation
          </p>
        </header>

        <section className="controls bg-gray-50 p-6 rounded-xl shadow-md mb-8">
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 items-start"> 
            {PARAMETER_CONFIGS.map(config => (
              <ParameterInput
                key={config.id}
                config={config}
                value={paramInputs[config.id]}
                onChange={handleParamChange}
                error={paramErrors[config.id]}
              />
            ))}
            <div className="flex flex-col">
              <label htmlFor="displayMethod" className="mb-1 text-sm font-medium text-gray-700">Display Plots</label>
              <select
                id="displayMethod"
                value={displayMethod}
                onChange={(e) => setDisplayMethod(e.target.value as DisplayMethod)}
                className="px-3 py-2.5 border border-gray-300 rounded-lg shadow-sm focus:outline-none focus:ring-2 focus:ring-indigo-500 focus:border-indigo-500 transition-colors bg-white text-gray-900"
                aria-label="Select plot display method"
              >
                <option value="all">All Methods</option>
                <option value="fd_explicit">Explicit FD (Forward Euler)</option>
                <option value="fd_implicit">Implicit FD (Backward Euler)</option>
                <option value="fd_crank_nicolson">Crank-Nicolson FD</option>
                <option value="mol">Method of Lines</option>
                <option value="fem">Finite Element Method</option>
                <option value="analytical">Analytical Solution</option>
              </select>
            </div>
            <button
              onClick={runSimulation}
              disabled={isLoading}
              className="col-span-1 md:col-span-2 lg:col-span-3 mt-4 h-12 bg-gradient-to-r from-indigo-600 to-purple-600 hover:from-indigo-700 hover:to-purple-700 text-white font-semibold py-2 px-6 rounded-lg shadow-md hover:shadow-lg transform hover:scale-105 transition-all duration-300 ease-in-out disabled:opacity-50 disabled:cursor-not-allowed flex items-center justify-center"
              aria-label="Run simulation and update plots"
            >
              {isLoading ? (
                <>
                  <svg className="animate-spin -ml-1 mr-3 h-5 w-5 text-white" xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" aria-hidden="true">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"></circle>
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"></path>
                  </svg>
                  Simulating...
                </>
              ) : "Run Simulation & Update Plots"}
            </button>
          </div>
        </section>
        
        {isLoading && (
          <div className="text-center py-10 text-xl font-medium text-indigo-600" role="status" aria-live="polite">
            🧬 Computing wound healing dynamics... Please wait.
          </div>
        )}

        {!isLoading && (
          <>
            <section className="plots-container grid grid-cols-1 md:grid-cols-2 gap-6 mb-8">
              {(displayMethod === 'all' || displayMethod === 'fd_explicit') && fdExplicitTransformed && ( 
                <SimulationPlot data={fdExplicitTransformed.plotData} timeLabels={fdExplicitTransformed.timeLabels} title="Explicit FD (Forward Euler) Method" />
              )}
              {(displayMethod === 'all' || displayMethod === 'fd_implicit') && fdImplicitTransformed && (
                <SimulationPlot data={fdImplicitTransformed.plotData} timeLabels={fdImplicitTransformed.timeLabels} title="Implicit FD (Backward Euler) Method" />
              )}
              {(displayMethod === 'all' || displayMethod === 'fd_crank_nicolson') && fdCNTransformed && (
                <SimulationPlot data={fdCNTransformed.plotData} timeLabels={fdCNTransformed.timeLabels} title="Crank-Nicolson FD Method" />
              )}
              {(displayMethod === 'all' || displayMethod === 'mol') && molTransformed && (
                 <SimulationPlot data={molTransformed.plotData} timeLabels={molTOutputLabels || []} title="Method of Lines" />
              )}
              {(displayMethod === 'all' || displayMethod === 'fem') && femTransformed && (
                 <SimulationPlot data={femTransformed.plotData} timeLabels={femTOutputLabels || []} title="Finite Element Method" />
              )}
              {(displayMethod === 'all' || displayMethod === 'analytical') && (
                isAnalyticalApplicable && analyticalTransformed && analyticalTransformed.plotData.length > 0 ? (
                  <SimulationPlot data={analyticalTransformed.plotData} timeLabels={analyticalTransformed.timeLabels} title="Analytical Solution" />
                ) : (
                  <div className="bg-white p-6 rounded-xl shadow-lg h-full flex flex-col items-center justify-center min-h-[300px]">
                    <h3 className="text-lg font-semibold text-center mb-2 text-gray-700">Analytical Solution</h3>
                    <p className="text-gray-500 text-sm text-center">Only available for p=0 and sc=0.</p>
                    {analyticalRawResultState && analyticalRawResultState.methodName.includes("(N/A") && 
                      <p className="text-xs text-gray-400 mt-2 text-center">Current params: p={params.p.toExponential(1)}, sc={params.sc.toExponential(1)}</p>
                    }
                  </div>
                )
              )}
              {displayMethod === 'all' && comparisonPlotData && comparisonPlotData.length > 0 && Object.keys(comparisonPlotData[0]).filter(k => k !== 'r').length > 0 && (
                <div className="md:col-span-2">
                   <SimulationPlot 
                    data={comparisonPlotData} 
                    timeLabels={comparisonLegendLabels} 
                    title={`Comparison at t = ${params.days} days`} 
                  />
                </div>
              )}
            </section>

            <section className="results bg-gray-50 p-6 rounded-xl shadow-md">
              <h3 className="text-xl font-semibold mb-3 text-gray-700">📊 Performance Analysis</h3>
              <pre className="bg-gray-800 text-white text-sm p-4 rounded-lg overflow-x-auto whitespace-pre" aria-label="Performance summary table">
                {`Method                         | Grid Points | Time Steps | CPU Time (s) | Rel. Error   | Accuracy (%) 
--------------------------------|-------------|------------|--------------|--------------|--------------
${performanceResultsText}`}
              </pre>
            </section>
          </>
        )}
      </div>
    </div>
  );
};

export default App;