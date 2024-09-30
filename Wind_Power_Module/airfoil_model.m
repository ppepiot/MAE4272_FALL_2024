% Airfoil lift and drag coefficient model

clc 
clear all
close all

%% Excel file for airfoil specific lift and drag coefficient
excelfilename = 'Excel_file_name.xlsx';

%% Input airfoil data all the way to stall
% Data can be found at airfoiltools

filename = 'AirfoilTools_generated_csv_file.csv';
startRow = 12;
formatSpec = '%25f%f%f%f%f%f%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

fclose(fileID);

output = [dataArray{1:end-1}];
clearvars filename startRow formatSpec fileID dataArray ans;

%% Lift Coefficient 

% The region between stall and flat plate theory requires some sort of 
% smooth curve fit. To do this, use the following:
% data_endpoint - last point from stall region
% data_endslope - slope at last point (using Euler method)
% flat_plate_start - first point from flate plate region (set to 30)
% flate_plate_slope - slope at first point from flat plate region
% Using these four constraints, a cubic function (y = ax^3+bx^2+cx+d that 
% that creates a smooth curve fit between the two regions (data up to stall 
% and flat plate) is generated. The coefficients (a,b,c,d) are found using
% fsolve.

% Note:
% While the region where flat plate theory starts has been hardcoded to 30
% degrees, this can be changed if you feel that the flate plate region 
% should start before/after 30 degrees.

% Get first two constraints from the airfoil data
data_endcl = output(end,2);
data_endclslope = (output(end,2) - output(end-1,2))/(output(end,1) - output(end-1,1));
data_endclangle = output(end,1);

% Get last two constraints from flat plate theory
alphacl = linspace(30,90,60*4);
flat_plate_clapprox = 2.*sind(alphacl).*cosd(alphacl);
flat_plate_clstart = 2.*sind(30)*cosd(30);
flat_plate_clslope = (flat_plate_clapprox(2) - flat_plate_clapprox(1))/0.25;
flat_plate_clangle = 30;

% Calculate the cubic equation coefficients using fsolve
fun = @(x)cubicfit(x,data_endclangle,data_endcl,data_endclslope,flat_plate_clangle,flat_plate_clstart,flat_plate_clslope);
x0 = [0,0,0,0];
x = fsolve(fun,x0);
a = x(1); b = x(2); c = x(3); d = x(4);

% Generate the curve fit
fit_rangecl = output(end,1):0.25:30;
fit_curvecl = a.*fit_rangecl.^3+b.*fit_rangecl.^2+c.*fit_rangecl+d;

% Plot the lift coefficient vs angle of attack plot 
figure(1);
plot(output(:,1),output(:,2),'LineWidth',2.0);
hold on;
plot(alphacl,flat_plate_clapprox,'LineWidth',2.0);
hold on;
plot(fit_rangecl,fit_curvecl,'LineWidth',2.0);
ylabel('$c_l$','interpreter','latex');
xlabel('$\alpha$ (degrees)','interpreter','latex');
legend('Data','Flat Plate','Fit');
axis([0 90 0 1.5]);
hold on;

%% Drag Coefficient

% The drag coefficient will work nearly identically to the lift coefficient

% Get first two constraints from the airfoil data
data_endcd = output(end,3);
data_endcdslope = (output(end,3) - output(end-1,3))/(output(end,1) - output(end-1,1));
data_endcdangle = output(end,1);

% Get last two constraints from flat plate theory
alphacd = linspace(30,90,60*4);
flat_plate_cdapprox = 2.*sind(alphacd).^2;
flat_plate_cdstart = 2.*sind(30).^2;
flat_plate_cdslope = (flat_plate_cdapprox(2) - flat_plate_cdapprox(1))/0.25;
flat_plate_cdangle = 30;

% Calculate the cubic equation coefficients using fsolve
fun = @(x)cubicfit(x,data_endcdangle,data_endcd,data_endcdslope,flat_plate_cdangle,flat_plate_cdstart,flat_plate_cdslope);
x0 = [0,0,0,0];
x = fsolve(fun,x0);
a = x(1); b = x(2); c = x(3); d = x(4);

% Generate the curve fit
fit_rangecd = output(end,1):0.25:30;
fit_curvecd = a.*fit_rangecd.^3+b.*fit_rangecd.^2+c.*fit_rangecd+d;

% Quick transpose of data
fit_rangecd = fit_rangecd';
fit_curvecd = fit_curvecd';
fit_rangecl = fit_rangecl';
fit_curvecl = fit_curvecl';
alphacd = alphacd';
alphacl = alphacl';
flat_plate_cdapprox = flat_plate_cdapprox';
flat_plate_clapprox = flat_plate_clapprox';

% Plot the lift coefficient vs angle of attack plot 
figure(2);
plot(output(:,1),output(:,3),'LineWidth',2.0);
hold on;
plot(alphacd,flat_plate_cdapprox,'LineWidth',2.0);
hold on;
plot(fit_rangecd,fit_curvecd,'LineWidth',2.0);
ylabel('$c_d$','interpreter','latex');
xlabel('$\alpha$ (degrees)','interpreter','latex');
legend('Data','Flat Plate','Fit');
axis([0 90 0 3]);
hold on;

%% Output the entire C_l,C_d vs alpha data into xlsx for Power Calc model

% Group up the angle of attack, C_l, and C_d data from airfoiltools, 
% the curve fit, and the flat plate theory together. This involves creating
% a zeropadded matrix and then concatanating the data. Not the most elegant
% way, but it works

% Obtain the lengths of the vectors
i = length(output(:,1))-1;
j  = length(fit_curvecl);
k = length(alphacd);

% Create the initial zeropadded matrix
M = zeros(i+j+k-1,3);

% Concatonate the vectors into the matrix
Headers = {'alpha' 'cl' 'cd'};
% Start with the angles of attack range
M(1:i,1) = output(1:i,1);
M(i+1:i+j,1) = fit_rangecl;
M(i+j+1:end,1) = alphacl(2:end);
% Move onto the lift coefficient vectors
M(1:i,2) = output(1:i,2);
M(i+1:i+j,2) = fit_curvecl;
M(i+j+1:end,2) = flat_plate_clapprox(2:end);
% Finish with the drag coefficient vectors
M(1:i,3) = output(1:i,3);
M(i+1:i+j,3) = fit_curvecd;
M(i+j+1:end,3) = flat_plate_cdapprox(2:end);

% Create the output with the headers
Z = [Headers; num2cell(M)];

% Write the matrix data into an xlsx file to put in the Power Calc model
writecell(Z,excelfilename);

%% Fsolve function

function F = cubicfit(x,dangle,d_endcl,ds,fpangle,fp_startcl,fp_slope)
F(1) = x(1)*dangle^3 + x(2)*dangle^2 + x(3)*dangle + x(4) - d_endcl;
F(2) = 3*x(1)*dangle^2 + 2*x(2)*dangle + x(3) - ds;
F(3) = x(1)*fpangle^3 + x(2)*fpangle^2 + x(3)*fpangle + x(4) - fp_startcl;
F(4) = 3*x(1)*fpangle^2 + 2*x(2)*fpangle + x(3) - fp_slope;
end
