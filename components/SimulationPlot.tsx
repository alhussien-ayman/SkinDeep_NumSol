
import React, { useMemo } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { PlotDataItem, LINE_COLORS } from '../types';

interface SimulationPlotProps {
  data: PlotDataItem[];
  timeLabels: string[];
  title: string;
}

const SimulationPlot: React.FC<SimulationPlotProps> = ({ data, timeLabels, title }) => {
  const isDataEffectivelyEmpty = !data || data.length === 0 || 
    (timeLabels.length > 0 && data.every(item => item[timeLabels[0]] === null || item[timeLabels[0]] === undefined));

  if (isDataEffectivelyEmpty) {
    if (!title.toLowerCase().includes("analytical solution")) {
         return (
            <div className="bg-white p-4 rounded-xl shadow-lg h-full flex flex-col items-center justify-center min-h-[300px] text-center">
                <h3 className="text-lg font-semibold mb-2 text-gray-700">{title}</h3>
                <p className="text-orange-600">
                    Simulation failed or produced unstable results.
                </p>
                <p className="text-xs text-gray-500 mt-1">
                    Parameters might be too extreme or time steps insufficient for this method.
                </p>
            </div>
        );
    }
    return <div className="text-center p-4 text-gray-500">No plottable data available for {title}.</div>;
  }

  const xTicks = useMemo(() => {
    if (!data || data.length === 0) return undefined;

    let maxR = 0;
    data.forEach(item => {
      if (typeof item.r === 'number' && item.r > maxR) {
        maxR = item.r;
      }
    });
    
    if (maxR === 0 && data.some(d => typeof d.r === 'number' && d.r === 0)) {
        return [0]; // Only a zero point data
    }
    if (maxR === 0) return undefined; // No valid r values or maxR couldn't be determined

    const tickInterval = 0.05;
    const generatedTicks = [];
    // Add a small epsilon (tickInterval * 0.1) to include maxR if it's a multiple of tickInterval.
    for (let tickValue = 0; tickValue <= maxR + tickInterval * 0.1; tickValue += tickInterval) {
      generatedTicks.push(parseFloat(tickValue.toFixed(2)));
    }
    
    return generatedTicks.length > 0 ? generatedTicks : undefined;
  }, [data]);

  const yTicks = useMemo(() => [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], []);


  return (
    <div className="bg-white p-4 rounded-xl shadow-lg h-full flex flex-col">
      <h3 className="text-lg font-semibold text-center mb-4 text-gray-700">{title}</h3>
      <div className="flex-grow">
        <ResponsiveContainer width="100%" height="100%">
          <LineChart data={data} margin={{ top: 5, right: 20, left: 15, bottom: 30 }}> {/* Adjusted left margin slightly for longer label */}
            <CartesianGrid strokeDasharray="3 3" stroke="#e0e0e0"/>
            <XAxis 
              dataKey="r" 
              type="number" 
              label={{ value: 'Position r (cm)', position: 'insideBottom', dy:15, fill:"#555" }} 
              stroke="#555"
              ticks={xTicks}
              domain={[0, 'auto']}
              tickFormatter={(tick) => typeof tick === 'number' ? tick.toFixed(2) : tick}
            />
            <YAxis 
              type="number" 
              domain={[0, 1.0]} 
              label={{ value: 'u(r,t) density (cell/cm³)', angle: -90, position: 'insideLeft', dx:-10, dy:40, fill:"#555" }}
              stroke="#555"
              ticks={yTicks}
              tickFormatter={(tick) => typeof tick === 'number' ? tick.toFixed(1) : tick}
            />
            <Tooltip formatter={(value: number | null) => value !== null ? value.toFixed(4) : "N/A"} />
            <Legend verticalAlign="top" wrapperStyle={{paddingBottom: "20px"}}/>
            {timeLabels.map((label, index) => (
              <Line
                key={label}
                type="monotone"
                dataKey={label}
                stroke={LINE_COLORS[index % LINE_COLORS.length]}
                strokeWidth={2}
                dot={false}
                isAnimationActive={false} 
                connectNulls={false} 
              />
            ))}
          </LineChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
};

export default SimulationPlot;
