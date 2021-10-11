close all;

% Data columns:
%   - Test ID
%   - Number of Orbits
%   - Number of Satellites
%   - Sim Processing Time [s]
%   - Vis Processing Time [s]
%   - Vis Runtime per Orbit [s]
%   - Sim Processing Time per Satellite [s]
%   - Vis Processing Time per Satellite [s]
data = [
  01	5	1	10.281	09.65	8.6	10.281	9.650
  02	5	1	10.313	09.98	8.6	10.313	9.980
  03	5	1	10.437	09.53	8.6	10.437	9.530
  04	5	2	12.017	10.73	8.8	06.009	5.365
  05	5	2	11.937	10.28	8.8	05.969	5.140
  06	5	2	12.078	10.56	8.8	06.039	5.280
  07	5	3	17.562	12.73	8.9	05.854	4.243
  08	5	3	17.626	11.74	8.9	05.875	3.913
  09	5	3	16.625	12.23	8.9	05.542	4.077
  10	5	4	23.687	13.49	9.2	05.922	3.373
  11	5	4	24.203	13.29	9.2	06.051	3.323
  12	5	4	23.891	13.59	9.2	05.973	3.398
  13	5	5	29.515	15.17	9.4	05.903	3.034
  14	5	5	29.468	14.66	9.4	05.894	2.932
  15	5	5	29.579	15.91	9.4	05.916	3.182
  16	5	6	35.563	17.35	9.6	05.927	2.892
  17	5	6	35.047	18.14	9.6	05.841	3.023
  18	5	6	36.016	18.38	9.6	06.003	3.063
  19	5	7	79.329	21.34	9.9	11.333	3.049
  20	5	7	75.750  24.21	9.9	10.821	3.459
  21	5	7	70.094	23.96	9.9	10.013	3.423
  22	5	8	185.516	28.35	9.9	23.190	3.544
  23	5	8	178.813	25.96	9.9	22.352	3.245
  24	5	8	181.407	27.64	9.9	22.676	3.455
  25	5	9	232.954	33.23	9.9	25.884	3.692
  26	5	9	234.262	33.08	9.9	26.029	3.676
  26	5	9	234.262	33.08	9.9	26.029	3.676
  27	5	10  281.672	37.65	9.9	28.167	3.765
  27	5	10  281.672	37.65	9.9	28.167	3.765
  27	5	10  281.672	37.65	9.9	28.167	3.765
];

maxNumSats = data(end,3); % Maximum Number of Satellites.

%% Plots
markerSize = 10;
colorsLight = [
  [0.9290, 0.6940, 0.1250]
  [0.4660, 0.6740, 0.1880]
  [0, 0.4470, 0.7410]
  [0.75, 0, 0.75]
  [1, 0, 0]
];
colorsDark = [
  [0.8500, 0.3250, 0.0980]
  [0, 0.5, 0]
  [0, 0, 1]
  [0.4940, 0.1840, 0.5560]
  [0.6350, 0.0780, 0.1840]
];

%% Total Processing Time [s].
figure('Visible','on');
figTitle = 'Total Processing Time';
xAxisLabel = 'Number of Satellites';
yAxisLabel = 'Processing Time [s]';

% Set array for plots and legends.
plots = [];
plegs = {};

% Simulation Processing Time.
xColumn = 3; % xLabel = 'Number of Satellites';
yColumn = 4; % yLabel = 'Simulation Processing Time [s]';
x = data(:,xColumn);
y = data(:,yColumn);
plot(x,y,'.','MarkerSize',markerSize,'Color',colorsLight(1,:));
% Average.
hold on;
x = 1:1:maxNumSats; % Number of Satellites.
y = zeros(size(x));
for n = 1:maxNumSats
  firstDataPoint = 3*n - 2;
  y(n) = ( ...
    data(firstDataPoint  ,yColumn) + ...
    data(firstDataPoint+1,yColumn) + ...
    data(firstDataPoint+2,yColumn) ) / 3;
end
xlim([0 x(end)+1]);
plots(end+1) = plot(x,y,'Color',colorsDark(1,:));
plegs{end+1} = 'Simulation';

% Visualization Processing Time.
hold on;
xColumn = 3; % xLabel = 'Number of Satellites';
yColumn = 5; % yLabel = 'Visualization Processing Time [s]';
x = data(:,xColumn);
y = data(:,yColumn);
plot(x,y,'.','MarkerSize',markerSize,'Color',colorsLight(2,:));
% Average.
hold on;
x = 1:1:maxNumSats; % Number of Satellites.
y = zeros(size(x));
for n = 1:maxNumSats
  firstDataPoint = 3*n - 2;
  y(n) = ( ...
    data(firstDataPoint  ,yColumn) + ...
    data(firstDataPoint+1,yColumn) + ...
    data(firstDataPoint+2,yColumn) ) / 3;
end
xlim([0 x(end)+1]);
plots(end+1) = plot(x,y,'Color',colorsDark(2,:));
plegs{end+1} = 'Visualization';

% Figure Setup.
title(figTitle);
xlabel(xAxisLabel);
ylabel(yAxisLabel);

% Plots & Legends.
legend(plots,plegs,'Location','north');


%% Processing Time per Satellite [s].
figure('Visible','on');
figTitle = 'Processing Time / Number of Satellites';
xAxisLabel = 'Number of Satellites';
yAxisLabel = 'Processing Time / Num of Satellites [s]';

% Set array for plots and legends.
plots = [];
plegs = {};

% Simulation Processing Time per Satellite.
xColumn = 3; % xLabel = 'Number of Satellites';
yColumn = 7; % yLabel = 'Sim Processing Time / Num of Satellites [s]';
x = data(:,xColumn);
y = data(:,yColumn);
plot(x,y,'.','MarkerSize',markerSize,'Color',colorsLight(1,:));
hold on;
% Compute Average.
x = 1:1:maxNumSats; % Number of Satellites.
y = zeros(size(x));
for n = 1:maxNumSats
  firstDataPoint = 3*n - 2;
  y(n) = ( ...
    data(firstDataPoint  ,yColumn) + ...
    data(firstDataPoint+1,yColumn) + ...
    data(firstDataPoint+2,yColumn) ) / 3;
end
xlim([0 x(end)+1]);
% Plot Average.
plots(end+1) = plot(x,y,'Color',colorsDark(1,:));
plegs{end+1} = 'Simulation';

% Visualization Processing Time per Satellite.
hold on;
xColumn = 3; % xLabel = 'Number of Satellites';
yColumn = 8; % yLabel = 'Vis Processing Time / Num of Satellites [s]';
x = data(:,xColumn);
y = data(:,yColumn);
plot(x,y,'.','MarkerSize',markerSize,'Color',colorsLight(2,:));
% Compute Average.
hold on;
x = 1:1:maxNumSats; % Number of Satellites.
y = zeros(size(x));
for n = 1:maxNumSats
  firstDataPoint = 3*n - 2;
  y(n) = ( ...
    data(firstDataPoint  ,yColumn) + ...
    data(firstDataPoint+1,yColumn) + ...
    data(firstDataPoint+2,yColumn) ) / 3;
end
xlim([0 x(end)+1]);
% Plot Average.
plots(end+1) = plot(x,y,'Color',colorsDark(2,:));
plegs{end+1} = 'Visualization';

% Figure Setup.
title(figTitle);
xlabel(xAxisLabel);
ylabel(yAxisLabel);

% Plots & Legends.
legend(plots,plegs,'Location','north');


%% Sim Processing Time per Satellite [s].
figure('Visible','off');
figTitle = 'Simulation Processing Time per Satellite';
xColumn = 3; xLabel = 'Number of Satellites';
yColumn = 7; yLabel = 'Sim Processing Time / Num of Satellites [s]';
x = data(:,xColumn);
y = data(:,yColumn);
plot(x,y,'.','MarkerSize',markerSize);
hold on;
% Average.
x = 1:1:maxNumSats; % Number of Satellites.
y = zeros(size(x));
for n = 1:maxNumSats
  firstDataPoint = 3*n - 2;
  y(n) = ( ...
    data(firstDataPoint  ,yColumn) + ...
    data(firstDataPoint+1,yColumn) + ...
    data(firstDataPoint+2,yColumn) ) / 3;
end
xlim([0 x(end)+1]);
plot(x,y);
title(figTitle);
xlabel(xLabel);
ylabel(yLabel);

%% Vis Processing Time per Satellite [s].
figure('Visible','off');
figTitle = 'Visualization Processing Time per Satellite';
xColumn = 3; xLabel = 'Number of Satellites';
yColumn = 8; yLabel = 'Vis Processing Time / Num of Satellites [s]';
x = data(:,xColumn);
y = data(:,yColumn);
plot(x,y,'.','MarkerSize',markerSize);
hold on;
% Average.
x = 1:1:maxNumSats; % Number of Satellites.
y = zeros(size(x));
for n = 1:maxNumSats
  firstDataPoint = 3*n - 2;
  y(n) = ( ...
    data(firstDataPoint  ,yColumn) + ...
    data(firstDataPoint+1,yColumn) + ...
    data(firstDataPoint+2,yColumn) ) / 3;
end
xlim([0 x(end)+1]);
plot(x,y);
title(figTitle);
xlabel(xLabel);
ylabel(yLabel);















