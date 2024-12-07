%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 
%% init
clear all; clc;close all;
addpath( "Functions","Nasa");
%% Units
mm      = 1e-3;dm=0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;
%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv Pref Tref
Runiv = 8.314472;                 % Universal gas constant
Pref = 1.01235e5;                 % Reference pressure, 1 atm!
Tref = 298.15;                    % Reference Temperature
[SpS,El]        = myload('Nasa\NasaThermalDatabase.mat',{'Diesel','O2','N2','CO2','H2O','NO','NO2','CO'});
TdataBase = fullfile('Nasa/NasaThermalDatabase.mat');
load(TdataBase);
whos ;
%% Engine geom data (check if these are correct)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 180;

%% Engine geom data for ideal
bore = 0.104;                     % Diameter of piston [m]
stroke = 0.085;                   % Full travel of piston [m]
conrodLength = 0.1365;            % Length of rod which connects the piston head and crankshaft [m]
cRatio = 21.5;                    % Ratio between min and max pressure[-]
displacement = 722e-6;            % The volume of the cylinder [m^3]
P0 = 101325;                      % Ambient Pressure [Pa]
r = stroke / 2;

% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355;
CaIVC = -135;
CaEVO = 149;
CaEVC = -344;
CaSOI = -3.2;
% Write a function [V] = CylinderVolume(Ca,Cyl) that will give you Volume
% for the given Cyl geometry. If you can do that you can create pV-diagrams
rho_HVO = 807;          % Density [kg/m^3]
rho_diesel = 836.1;
cn = 63;                        % Cetane ratio, how ignitable it is [-]
eta = 3.338e-6;                  % Viscosity [m2/kg]
%% 

%% Assumed data Ingore these for stoichiometric analyis except lhv
lhv = [43e6;43.292e6];                % Lower heating value, energy content per kg of fuel [J/kg] left pure diesel, right FAME
RPMn = 1000;                      % Unloaded RPM of a typical diesel engine
rotTime = 60 / RPMn;              % Time of a full rotation of 360 degrees [s]
R_specific = 287;                 % Specific gas constant for air [J/(kg·K)] at higher temperatures
md = 0.00005;                     % Fuel (diesel) mass per cycle [kg]
%Cp = 2.75e3;                      % Specific heat capacity at high temperature [J/kg·K]
gamma = 1.4;                      % Gamma for air-fuel mixture
ma = 14.5 * md;                   % Air-fuel mass ratio (14.5:1) 
mtot = ma + md;                   % Total mass (air + fuel)
%% Timing settings (degrees)
intakeLength = 180;
intakeRange = intakeLength;

comprLength = 180;
comprRange = intakeRange + comprLength;

combLength  = 30;
combRange = comprRange + combLength;

expLength = 150;
expRange = combRange + expLength;

exhaustLength = 180;
exhaustRange = expRange + exhaustLength;

if exhaustRange <720
    error("Timing does not sum to a full rotation.")
end

%% Initializing
% Initialize pressure and temperature
pi(1) = Pref;
Ti(1) = Tref;
for i = 1:1:720
    t(i) = i * (rotTime) / 360;
    ca(i) = i;  % Crank angle in degrees, from 1 to 720
end

W_inst = zeros(1, 720);                                         % Instantaneous work done at each crank angle

%% Modelling
for n = 2:1:720
     Vi(n) = Volume(ca(n), cRatio, conrodLength, bore, displacement, r);

    switch true
        case (ca(n) == 1)
            pi(n) = pi(1);
            Ti(n) = Ti(1);
            Vi(n) = Vi(1);

        %% Intake Stroke
        case (ca(n) <= intakeRange)
            pi(n) = pi(1);                                                                                 
            Ti(n) = Ti(1);    
            Vi(n) = Vi(n);

        %% Compression Stroke
        case (ca(n) > intakeRange && ca(n) <= comprRange)
              %Elements 
            iSp = myfind({Sp.Name},{ 'Diesel','O2','N2','CO2','H2O','NO','NO2','CO'});
            
            SpS = Sp(iSp);
            NSp = length(SpS);
            Mi = [SpS.Mass];

            %Universal gas constant int
            R_O2 = Runiv/(Mi(2));
            R_N2 = Runiv/(Mi(3));
            R_Air = (0.79*R_O2)+(0.21*R_N2);
            
            %Cp int
            CpO2 = CpNasa(Ti(n-1),SpS(1));
            CpN2 = CpNasa(Ti(n-1),SpS(2));
            CpAir = (0.79*CpO2)+(0.21*CpN2);
            
            %Gamma int
            gamma = CpAir/(CpAir-R_Air);
            
            %Adiabatic relations int
            Ti(n) = Ti(n-1) * (Vi(n-1) / Vi(n))^(gamma - 1);                                                 
            pi(n) = pi(n-1) * (Vi(n-1) / Vi(n))^gamma;

          

        %% Combustion Phase
        case (ca(n) > comprRange && ca(n) <= combRange)
            
            %Cp int air
            CpO2 = CpNasa(Ti(n-1),SpS(2));
            CpN2 = CpNasa(Ti(n-1),SpS(3));
            CpAir = (0.79*CpO2)+(0.21*CpN2);
            
            %Cp int fuel, stays the same for HVO50
            CpDiesel = CpNasa(Ti(n-1),SpS(1));
            
            %M int
            M_air = (0.79*Mi(2))+(0.21*Mi(3));
            M_Diesel = 167;  % Fossil diesel molar mass
            
            
            %Mass fractions int
            Mf_Air = M_air/(M_Diesel+M_air);
            Mf_Diesel = M_Diesel/(M_Diesel+M_air);

            %dQ int, Assume adiabatic.
            dQ_comb = sum(lhv .* [md * 0.5; md * 0.5]); 
                % Heat released from 50% pure diesel and 50% HVO
                % Note, this is the potential heat release over a complete
                % combustion!
            dQ(n) = 0;                                              % For easier calculations, assume that there is no heat exchange in this small time instance
           
            %Mixture Cp int
            CpComb = (Mf_Air*CpAir) + (Mf_Diesel*CpDiesel);
            
            %Pressure int, assume isobaric.
            pi(n) = pi(n-1);
            Ti(n) = Ti(n-1) + (Vi(n)*(pi(n-1)-pi(n)) + dQ_comb./combLength - dQ(n)) / (mtot * CpComb);        % Based on the first law of thermodynamics
                %Note, dp = 0, dq = 0, so ONLY CpComb, dQ_comb as variables
                %that change.
            
        %% Expansion Stroke
        case (ca(n) > combRange && ca(n) <= expRange)
            
            % Defining moles of each fraction
            w_HVO = displacement * rho_HVO;   
            w_Diesel = displacement * rho_diesel; 
        
            
            M_HVO = 226;      
            M_Diesel = 167;   
            
            
            m_HVO = w_HVO * md;        
            m_Diesel = w_Diesel * md;  
            
            n_HVO = m_HVO / M_HVO;          
            n_Diesel = m_Diesel / M_Diesel;  


            R_CO2 = Runiv/(Mi(4));
            R_H2O= Runiv/(Mi(5));
             
            
            
            O2_to_HVO = 24.5;         
            O2_to_Diesel = 17.75;     
            
            % Total oxygen required for complete combustion
            molO2 = n_HVO * O2_to_HVO + n_Diesel * O2_to_Diesel;
            
            
            % Products
            CO2_to_HVO = 16;        
            CO2_to_Diesel = 12;     
            H2O_to_HVO = 17;        
            H2O_to_Diesel = 11.5;   
            
            molCO2 = n_HVO * CO2_to_HVO + n_Diesel * CO2_to_Diesel;
            molH2O = n_HVO * H2O_to_HVO + n_Diesel * H2O_to_Diesel;

           
            molAir = ma / M_air;  
            molO2Air = 0.79 * molAir;
            molN2 = 0.21 * molAir;
            
            % Unreacted O2
            molO2Final = molO2Air - molO2;

            % Total moles after combustion
            molTot = molCO2 + molH2O + molO2Final + molN2;
            
            % Mole fractions
            mfO2 = molO2Final / molTot;
            mfCO2 = molCO2 / molTot;
            mfH2O = molH2O / molTot;
            mfN2 = molN2 / molTot;
           
            CpO2 = CpNasa(Ti(n-1),SpS(2));
            CpN2 = CpNasa(Ti(n-1),SpS(3));
            CpCO2 = CpNasa(Ti(n-1),SpS(4));
            CpH2O = CpNasa(Ti(n-1),SpS(5));

            RExp = mfO2*R_O2 + mfH2O*R_H2O + mfCO2*R_CO2 + mfN2*R_N2;
            CpExp = mfO2*CpO2 + mfH2O*CpH2O + mfCO2*CpCO2 + mfN2*CpN2;
            gamma = CpExp/(CpExp-RExp);
            
            Ti(n) = Ti(n-1) * (Vi(n-1) / Vi(n))^(gamma - 1);                                                 
            pi(n) = pi(n-1) * (Vi(n-1) / Vi(n))^gamma;

        %% Exhaust Phase
        case (ca(n) > expRange && ca(n) <= exhaustRange)
            pi(n) = pi(1);
            Ti(n) = Ti(1);

        otherwise
            disp('Crank angle out of range (0 to 720 degrees).');
    end

    %% Instantaneous Work Done
    if n > 2
        W_inst(n) = trapz(Vi(n-1:n), pi(n-1:n));                      % Work done between crank angles n-1 and n
    end
end
%% Variables for KPI's
    
    
Mi = [SpS.Mass];  % Get molar masses from nasa database
MiNox = (Mi(7)+9*Mi(6))/10; % calculates the molar mass of NOx with ratio 1:9 of NO2:NO
    
NOXmatrix = [10,12,22,27,30,34,18,30,35,6,20,47,15,13,0,0,0,0,0,0,0,0,0,0,0,0]; %NOX emissions in ppm per test
    
NoxMole = NOXmatrix/1000000; % Calculates the moles of Nox using the matrix
CoMole = 0; % The amount of moles of CO in the emission
    
emissions = [NoxMole,CoMole]; % Input for the KPI function
MiEmission = [Mi(3),Mi(2),Mi(4),Mi(8),MiNox,Mi(5),Mi(1)]; % Molar masses of compounds [N2,O2,CO2,CO,NOx,H2O,Diesel]
    
LHV = 42.6e6; % Lower heating value of the fuel in joule per kilogram
AFRsto = 14.5; % Stoichiometric Air Fuel ratio of the fuel
    

%% Loading all Data

folderPath = 'Data';

% Get a list of all the files in the folder
fileList = dir(fullfile(folderPath, '*.txt'));

% Make structure for data
dataStruct = struct();

% Loop through each file
    % NOTE - One file represents the following:
    % Nr before Test = # test
    % nr before Load= Procent of maximum load
    % nr before Ca = degrees before TDC of injection
        
for k = 1:length(fileList)
    % Get the full file path
    filePath = fullfile(folderPath, fileList(k).name);
    
    % Put file in an array/matrix
    dataMatrix = readmatrix(filePath); % For .txt files with numeric data
    
    fileName = fileList(k).name;
    if ~strcmp(fileName, 'ExampleDataSet.txt')
        CA_value = regexp(fileName, '_(\d+)Ca','tokens'); % Read file name for crankshaft angle of injection
        CA_value = str2double(CA_value{1}{1});
        Load_value = regexp(fileName, '_(\d+)Load','tokens'); % Read file name for Load
        Load_value = str2double(Load_value{1}{1});
    else
        CA_value = 'idk';
        Load_value = 'idk';
    end

    
    % Putting the name, data matrix, crankshaft angle of injection, Load in the data structure
    dataStruct(k).fileName = fileList(k).name;
    dataStruct(k).data = dataMatrix;
    dataStruct(k).CA = CA_value;  % This is actually the negative crankshaft angle of injection, but value is positive here
    dataStruct(k).Load = Load_value;

end

%% Loop for calculating

ProcessedStruct = struct();

for k = 1:length(fileList)


    % Load data (if txt file)
    
    [Nrows,Ncols]   = size(dataStruct(k).data);                    % Determine size of array
    NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
    Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
    Ca              = reshape(dataStruct(k).data(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
    p               = reshape(dataStruct(k).data(:,2),[],Ncycles)*bara; % type 'help reshape' in the command window if you want to know what it does (reshape is a Matlab buit-in command
    FuelFlow        = reshape(dataStruct(k).data(:,4),[],Ncycles);
    
    % if k ==1
    %     plot(Ca,FuelFlow)
    % end
    % Calculate total work and fuel flow
    V = CylinderVolume(Ca,Cyl);

    for i = 1:Ncycles
        WorkpCycle(i,:) = trapz(V(:,i),p(:,i));  % Calculates the work per cycle
    end
    
    Worktot = sum(WorkpCycle);  % Calculates the total work
    % Fueltot= sum(FuelFlow,[1 2]); %Calculates the total work
    
    workavg = Worktot/Ncycles; % Average work per cycle
    fuelavg = 0.17; %Average grams per second fuel flow rate
    

    % Processing the Data
    
    [efficiency,bsfc,bsCo2,bsNox,equivalence] = KPIs(emissions,fuelavg/1000,workavg,Cyl,AFRsto,LHV,MiEmission,[10.8,18.7]);

    % Putting the Data in a structure
    ProcessedStruct(k).WorkperCycle = workavg;
    ProcessedStruct(k).Efficiency = efficiency;
    ProcessedStruct(k).bsfc = bsfc;
    ProcessedStruct(k).bsCO2 = bsCo2;
    ProcessedStruct(k).bsNOx = bsNox;
    ProcessedStruct(k).Equivalence = equivalence;

end














%% Plotting 
    %% Average Work, Power, and BSFC over the cycle
    m_fuel_dot = md  /RPMn *60*2; % Fuel mass flow rate
    W_total = sum(W_inst); % Total work done over the cycle
    eta = W_total/(sum(lhv.*[md*0.5;md*0.5]));
    

%% Display results
    fprintf('The mass fuel flow rate %.10f \n', m_fuel_dot);
    fprintf('The total work done during the cycle is %.2f J\n', W_total);
    fprintf('The efficiency of the cycle is %.2f \n', eta);

    % Plotting Instantaneous Work vs Crank Angle
    figure(3);
    subplot(1, 2, 1);
    plot(ca, W_inst, 'LineWidth', 2);
    xlabel('Crank Angle (degrees)');
    ylabel('Instantaneous Work (J)');
    title('Instantaneous Work vs. Crank Angle');
    grid on;

    %% Plotting P-V Diagram



    %% Plotting Temperature vs Crank Angle
    subplot(1, 2, 2);
    plot(ca, Ti, 'LineWidth', 2);
    hold on;

    xlabel('Crank Angle (degrees)');
    ylabel('Temperature (K)');
    title('Temperature vs Crank Angle');
    grid on;


%f1=figure(1);
pp = plot(Ca,p/bara,'LineWidth',1);                 % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                    % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('All cycles in one plot.')
%% pV-diagram

f1 = figure(1); 

% Plot the P-V Diagram (first plot)
grid on;
plot(Vi(2:end)/dm^3, pi(2:end)/bara, 'LineWidth', 2);  % P-V diagram

% Hold the plot to overlay the next plot (pV-diagram)
hold on;

% Plot the pV-diagram (second plot)
plot(V/dm^3, p(:,iselect)/bara, 'LineWidth', 2);  % pV-diagram

% Add labels and title
xlabel('Volume (dm^3)');
ylabel('Pressure (bar)');
title('P-V Diagram and pV-Diagram');
legend('P-V curve', 'pV-Diagram', 'Location', 'best');

% Set grid and limits for both plots
grid on;
%xlim([0 0.8]);  % Adjust x-axis range for both plots
%ylim([0.5 50]); % Adjust y-axis range for both plots
set(gca, 'XTick', [0:0.1:0.8], 'YGrid', 'on', 'XGrid', 'on');

% Remove subplot and figure(2) since both plots are now combined in figure(1)
%% Log-Log Plot of Volume vs Pressure
% Create a figure for the Log-Log P-V curve and pV-diagram
figure(2);

% Plot the Log-Log P-V Curve in the first subplot
 
loglog(Vi(2:end)/dm^3, pi(2:end)/bara, 'LineWidth', 2);  % Log-Log P-V curve
hold on
grid on;

% Plot the pV-diagram in the second subplot
loglog(V/dm^3, p(:,iselect)/bara, 'LineWidth', 2);  % pV-diagram
xlabel('V [dm^3]');
ylabel('p [bar]');  % Always add axis labels
%xlim([0.02 0.8]);
%ylim([0 50]);  % Adjust as needed
set(gca, 'XTick', [0.02 0.05 0.1 0.2 0.5 0.8], 'YTick', [0.5 1 2 5 10 20 50], 'XGrid', 'on', 'YGrid', 'on');
title('loglog pV-diagram');


%% Plotting, Interesting results

%We have the following data:
%%Fuel type (test nr.), Crank angle, Load case.
%%-> The following derived information: Work per cycle, Efficiency,
    % BrakeSpecificFuelConsumption, Brake Specific CO2, Brake Specific NOx,
    % Equivalence (Air-Fuel)


%% All the below code has a significant drawback. I need to manually enter the indices that I want to consult.
% This will make data analysis really slow, and I am not quite sure how to
% implement this in an improved manner.

fECU = figure("Name", "ECU observations");
%%Observing the ECU doing its thing.
%Equivalence vs. Crank angle for various loads is interesting
    % Does the engine run leaner (i.e. underfueled) under high load?

    
subplot(2, 2, 1);
hold on
scatter(dataStruct(1).Load(1),ProcessedStruct(1).Equivalence(1), 'r')
scatter(dataStruct(2).Load(1),ProcessedStruct(2).Equivalence(1), 'g')
scatter(dataStruct(3).Load(1),ProcessedStruct(3).Equivalence(1), 'b')
hold off; 
xlabel("Load, percentage of maximum")
ylabel("Air/Fuel ratio?")
title('Effect of load on air/fuel ratio');

%Equivalence vs. Crank angle for various ignition angles is interesting
    % Does it get less stochiometric due to ignition timing?

subplot(2, 2, 2);
hold on
scatter(dataStruct(2).CA(1),ProcessedStruct(2).Equivalence(1), 'r')
scatter(dataStruct(12).CA(1),ProcessedStruct(12).Equivalence(1), 'g')
scatter(dataStruct(14).CA(1),ProcessedStruct(14).Equivalence(1), 'b')
hold off; 
xlabel("Injection at n degrees of crank angle before TDC")
ylabel("Air/Fuel ratio?")
title('Effect of Crank Angle on air/fuel ratio.');

%%Observing the influence of parameters on efficiency -> Set engine parameters.

%Efficiency vs. Crank angle for various different injection timings.
    % Does retarding timing actually increase efficiency?

subplot(2, 2, 3);
hold on
scatter(dataStruct(2).CA(1),ProcessedStruct(2).Efficiency(1), 'r')
scatter(dataStruct(12).CA(1),ProcessedStruct(12).Efficiency(1), 'g')
scatter(dataStruct(14).CA(1),ProcessedStruct(14).Efficiency(1), 'b')
hold off; 
ylabel("Efficiency")
xlabel("Injection at n degrees of crank angle before TDC")
title('Effect of Crank Angle on efficiency.');

%Efficiency vs. Crank angle for various different loads.
    % Does the engine run more efficiently at higher power?

subplot(2, 2, 4);
hold on
scatter(dataStruct(1).Load(1),ProcessedStruct(1).Efficiency(1), 'r')
scatter(dataStruct(2).Load(1),ProcessedStruct(2).Efficiency(1), 'g')
scatter(dataStruct(3).Load(1),ProcessedStruct(3).Efficiency(1), 'b')
ylabel("Efficiency")
xlabel("Load, percentage of maximum")
hold off; 
title('Effect of load on efficiency.');

%%Observing the influence of parameters on emissions
%NOx, CO2 vs. Load
    % Expected behavior is that NOx should increase, CO2 constant for a
    % good ECU controller

fEmissions = figure("Name", "Emissions vs. engine parameters");

subplot(2, 2, 1);
hold on

scatter(dataStruct(1).Load(1),ProcessedStruct(1).bsCO2(1), 'g', "DisplayName", "CO2")
scatter(dataStruct(2).Load(1),ProcessedStruct(2).bsCO2(1), 'g', "DisplayName", "CO2")
scatter(dataStruct(3).Load(1),ProcessedStruct(3).bsCO2(1), 'g', "DisplayName", "CO2")

legend
ylabel("Emissions [unit].")
xlabel("Load, percentage of maximum")
hold off; 
title('Effect of load on emissions.');

subplot(2, 2, 2);
hold on


scatter(dataStruct(1).Load(1),ProcessedStruct(1).bsNOx(1), 'r', "DisplayName", "NOx")
scatter(dataStruct(2).Load(1),ProcessedStruct(2).bsNOx(1), 'r', "DisplayName", "NOx")
scatter(dataStruct(3).Load(1),ProcessedStruct(3).bsNOx(1), 'r', "DisplayName", "NOx")

legend
ylabel("Emissions [unit].")
xlabel("Load, percentage of maximum")
hold off; 
title('Effect of load on emissions.');


%NOx, CO2 vs. Crank angle is interesting
    % Can we see the predicted behavior of a lower CA (before TDC) causing
    % a drop in NOx?
% subplot(2, 2, 2);
% hold on
% 
% scatter(dataStruct(2).CA(1),ProcessedStruct(2).bsCO2(1), 'g', "DisplayName", "CO2")
% scatter(dataStruct(12).CA(1),ProcessedStruct(12).bsCO2(1), 'g', "DisplayName", "CO2")
% scatter(dataStruct(14).CA(1),ProcessedStruct(14).bsCO2(1), 'g', "DisplayName", "CO2")
% 
% scatter(dataStruct(2).CA(1),ProcessedStruct(2).bsNOx(1), 'r', "DisplayName", "NOx")
% scatter(dataStruct(12).CA(1),ProcessedStruct(12).bsNOx(1), 'r', "DisplayName", "NOx")
% scatter(dataStruct(14).CA(1),ProcessedStruct(14).bsNOx(1), 'r', "DisplayName", "NOx")
% 
% legend
% ylabel("Emissions [unit].")
% xlabel("Injection at n degrees of crank angle before TDC")
% hold off; 
% title('Pre-injection effect on emissions.');


%Experiment using Copilot to create a polynomial fit.

subplot(2, 2, 3);
hold on

% CO2 scatter points
scatter(dataStruct(2).CA(1), ProcessedStruct(2).bsCO2(1), 'g', "DisplayName", "CO2");
scatter(dataStruct(12).CA(1), ProcessedStruct(12).bsCO2(1), 'g', "DisplayName", "CO2");
scatter(dataStruct(14).CA(1), ProcessedStruct(14).bsCO2(1), 'g', "DisplayName", "CO2");

% CO2 polynomial fit (degree 2)
x_CO2 = [dataStruct(2).CA(1), dataStruct(12).CA(1), dataStruct(14).CA(1)];
y_CO2 = [ProcessedStruct(2).bsCO2(1), ProcessedStruct(12).bsCO2(1), ProcessedStruct(14).bsCO2(1)];
p_CO2 = polyfit(x_CO2, y_CO2, 2); % Quadratic fit
xfit_CO2 = linspace(min(x_CO2), max(x_CO2), 100); % Generate x values for smooth curve
yfit_CO2 = polyval(p_CO2, xfit_CO2);
plot(xfit_CO2, yfit_CO2, 'g--', "DisplayName", "CO2 Fit");

legend;
ylabel("Emissions [unit]");
xlabel("Injection at n degrees of crank angle before TDC");
title('Pre-injection effect on emissions');
hold off;

% NOx scatter points
subplot(2, 2, 4);
hold on
scatter(dataStruct(2).CA(1), ProcessedStruct(2).bsNOx(1), 'r', "DisplayName", "NOx");
scatter(dataStruct(12).CA(1), ProcessedStruct(12).bsNOx(1), 'r', "DisplayName", "NOx");
scatter(dataStruct(14).CA(1), ProcessedStruct(14).bsNOx(1), 'r', "DisplayName", "NOx");

% NOx polynomial fit (degree 2)
x_NOx = [dataStruct(2).CA(1), dataStruct(12).CA(1), dataStruct(14).CA(1)];
y_NOx = [ProcessedStruct(2).bsNOx(1), ProcessedStruct(12).bsNOx(1), ProcessedStruct(14).bsNOx(1)];
p_NOx = polyfit(x_NOx, y_NOx, 2); % Quadratic fit
xfit_NOx = linspace(min(x_NOx), max(x_NOx), 100); % Generate x values for smooth curve
yfit_NOx = polyval(p_NOx, xfit_NOx);
plot(xfit_NOx, yfit_NOx, 'r--', "DisplayName", "NOx Fit");

legend;
ylabel("Emissions [unit]");
xlabel("Injection at n degrees of crank angle before TDC");
title('Pre-injection effect on emissions');
hold off;

function Vi = Volume(ca, cRatio, conrodLength, bore, displacement, r)
    % Input parameters:
    % ca: Crank angle (in degrees)
    % cRatio: Compression ratio
    % conrodLength: Length of the connecting rod (m)
    % bore: Bore diameter of the engine (m)
    % displacement: Total displacement volume (m^3)

    % Convert crank angle to radians
    ca = deg2rad(ca);

    % Clearance volume based on compression ratio
    Vc = displacement / (cRatio - 1);

    % Volume calculation at a given crank angle
    Vi = Vc + (pi / 4) * bore.^2 * (r + conrodLength - (r .* cos(ca) + sqrt(conrodLength .^2 - r^2 * sin(ca) .^2)));
end