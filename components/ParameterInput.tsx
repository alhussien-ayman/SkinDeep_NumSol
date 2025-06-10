
import React from 'react';
import { ParameterConfig, SimulationParams } from '../types';

interface ParameterInputProps {
  config: ParameterConfig;
  value: string | number; 
  onChange: (id: keyof SimulationParams, value: string) => void;
  error?: string | null; // New prop for error message
}

const ParameterInput: React.FC<ParameterInputProps> = ({ config, value, onChange, error }) => {
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    onChange(config.id, e.target.value);
  };

  const inputMode = config.type === 'text' ? 'text' : 'decimal';

  const baseClasses = "px-3 py-2 border rounded-lg shadow-sm focus:outline-none focus:ring-2 transition-colors bg-white text-gray-900";
  const errorClasses = "border-red-500 focus:ring-red-500 focus:border-red-500";
  const normalClasses = "border-gray-300 focus:ring-indigo-500 focus:border-indigo-500";

  return (
    <div className="flex flex-col">
      <label htmlFor={config.id} className="mb-1 text-sm font-medium text-gray-700">
        {config.label}
      </label>
      <input
        type="text" 
        inputMode={inputMode} 
        id={config.id}
        name={config.id}
        value={value} 
        onChange={handleChange}
        step={config.step} 
        min={config.min} 
        max={config.max}
        className={`${baseClasses} ${error ? errorClasses : normalClasses}`}
        aria-invalid={!!error}
        aria-describedby={error ? `${config.id}-error` : undefined}
      />
      {error && (
        <p id={`${config.id}-error`} className="mt-1 text-xs text-red-600" role="alert">
          {error}
        </p>
      )}
    </div>
  );
};

export default ParameterInput;