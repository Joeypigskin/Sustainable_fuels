clc, clear, close all
%% Initial Parameters
bore = 0.104;                     % Diameter of piston [m]
stroke = 0.085;                   % Full travel of piston [m]
conrodLength = 0.1365;            % Length of rod which connects the piston head and crankshaft [m]
cRatio = 21.5;                    % Ratio between min and max pressure[-]
displacement = 722e-6;            % The volume of the cylinder [m^3]
P0 = 101325;                      % Ambient Pressure [Pa]
r = stroke / 2;

TdataBase = fullfile('Nasa/NasaThermalDatabase.mat');
load(TdataBase);
whos ;

global Runiv Pref Tref
Runiv = 8.314472;                 % Universal gas constant
Pref = 1.01235e5;                 % Reference pressure, 1 atm!
Tref = 298.15;                    % Reference Temperature

% Diesel properties : HVO50 will be used

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

% Crank angle and time setup
for i = 1:1:720
    t(i) = i * (rotTime) / 360;
    ca(i) = i;  % Crank angle in degrees, from 1 to 720
end

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
p(1) = P0;
T(1) = Tref;


W_inst = zeros(1, 720);                                         % Instantaneous work done at each crank angle

%% Modelling
for n = 2:1:720
    V(n) = Volume(ca(n), cRatio, conrodLength, bore, displacement, r);      % Volume at crank angle

    switch true
        case (ca(n) == 1)
            p(n) = p(1);
            T(n) = T(1);
            V(n) = V(1);

        %% Intake Stroke
        case (ca(n) <= 180)
            p(n) = p(1);                                                                                 
            T(n) = T(1);    
            V(n) = V(n);

        %% Compression Stroke
        case (ca(n) > 180 && ca(n) <= 359)
              %Elements 
            iSp = myfind({Sp.Name},{'O2','N2','CO2','H2O','Diesel'});
            SpS = Sp(iSp);
            NSp = length(SpS);
            Mi = [SpS.Mass];

            %Universal gas constant int
            R_O2 = Runiv/(Mi(1));
            R_N2 = Runiv/(Mi(2));
            R_Air = (0.79*R_O2)+(0.21*R_N2);
            
            %Cp int
            CpO2 = CpNasa(T(n-1),SpS(1));
            CpN2 = CpNasa(T(n-1),SpS(2));
            CpAir = (0.79*CpO2)+(0.21*CpN2);
            
            %Gamma int
            gamma = CpAir/(CpAir-R_Air);
            
            %Adiabatic relations int
            T(n) = T(n-1) * (V(n-1) / V(n))^(gamma - 1);                                                 
            p(n) = p(n-1) * (V(n-1) / V(n))^gamma;

          

        %% Combustion Phase
        case (ca(n) > 359 && ca(n) <= 390)
            
            %Cp int air
            CpO2 = CpNasa(T(n-1),SpS(1));
            CpN2 = CpNasa(T(n-1),SpS(2));
            CpAir = (0.79*CpO2)+(0.21*CpN2);
            
            %Cp int fuel, stays the same for HVO50
            CpDiesel = CpNasa(T(n-1),SpS(5));
            
            %M int
            M_air = (0.79*Mi(1))+(0.21*Mi(2));
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
            p(n) = p(n-1);
            T(n) = T(359) + (V(n)*(p(n-1)-p(n)) + dQ_comb./combLength - dQ(n)) / (mtot * CpComb);        % Based on the first law of thermodynamics
                %Note, dp = 0, dq = 0, so ONLY CpComb, dQ_comb as variables
                %that change.
            
        %% Expansion Stroke
        case (ca(n) > 390 && ca(n) <= 540)
            
            % Defining moles of each fraction
            w_HVO = displacement * rho_HVO;   
            w_Diesel = displacement * rho_diesel; 
        
            
            M_HVO = 226;      
            M_Diesel = 167;   
            
            
            m_HVO = w_HVO * md;        
            m_Diesel = w_Diesel * md;  
            
            n_HVO = m_HVO / M_HVO;          
            n_Diesel = m_Diesel / M_Diesel;  


            R_CO2 = Runiv/(Mi(3));
            R_H2O= Runiv/(Mi(4));
             
            
            
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
            
            CpO2 = CpNasa(T(n-1),SpS(1));
            CpN2 = CpNasa(T(n-1),SpS(2));
            CpCO2 = CpNasa(T(n-1),SpS(3));
            CpH2O = CpNasa(T(n-1),SpS(4));

            RExp = mfO2*R_O2 + mfH2O*R_H2O + mfCO2*R_CO2 + mfN2*R_N2;
            CpExp = mfO2*CpO2 + mfH2O*CpH2O + mfCO2*CpCO2 + mfN2*CpN2;
            gamma = CpExp/(CpExp-RExp);
            
            T(n) = T(n-1) * (V(n-1) / V(n))^(gamma - 1);                                                 
            p(n) = p(n-1) * (V(n-1) / V(n))^gamma;

        %% Exhaust Phase
        case (ca(n) > 540 && ca(n) <= 720)
            p(n) = p(1);
            T(n) = T(1);

        otherwise
            disp('Crank angle out of range (0 to 720 degrees).');
    end

    %% Instantaneous Work Done
    if n > 2
        W_inst(n) = trapz(V(n-1:n), p(n-1:n));                      % Work done between crank angles n-1 and n
    end
end

%% Average Work, Power, and BSFC over the cycle
m_fuel_dot = md  /RPMn *60*2; % Fuel mass flow rate
W_total = sum(W_inst); % Total work done over the cycle
eta = W_total/(sum(lhv.*[md*0.5;md*0.5]));

plot_results(ca, W_inst, m_fuel_dot, W_total, eta, V, p, T);

%% Volume Function
function V = Volume(ca, cRatio, conrodLength, bore, displacement, r)
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
    V = Vc + (pi / 4) * bore^2 * (r + conrodLength - (r * cos(ca) + sqrt(conrodLength^2 - r^2 * sin(ca)^2)));
end
function plot_results(ca, W_inst, m_fuel_dot, W_total, eta, V, p, T)
    % This function plots various results from the engine cycle analysis.
    % Inputs:
    % - ca: Crank angle (degrees)
    % - W_inst: Instantaneous work (J)
    % - m_fuel_dot: Mass fuel flow rate
    % - W_total: Total work done (J)
    % - eta: Efficiency of the cycle
    % - V: Volume (m^3)
    % - p: Pressure (Pa)
    % - T: Temperature (K)
    
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
    figure(1);
    subplot(1, 2, 1);
    grid on;
    plot(V(2:end), p(2:end), 'LineWidth', 2);
    hold on;

    % Annotate stroke phases on P-V diagram
    text(V(1), p(1), 'Intake', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(V(180), p(180), 'Compression', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
    text(V(359), p(359), 'Combustion', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);

    % Mark transition points on P-V diagram (ca. 180, 360, 540)
    plot(V(180), p(180), 'ro', 'MarkerFaceColor', 'r');
    plot(V(360), p(360), 'go', 'MarkerFaceColor', 'g');
    plot(V(540), p(540), 'bo', 'MarkerFaceColor', 'b');

    grid off;
    xlabel('Volume (m^3)');
    ylabel('Pressure (Pa)');
    title('P-V Diagram');
    legend('P-V curve', 'Location', 'best');

    %% Plotting Temperature vs Crank Angle
    subplot(1, 2, 2);
    plot(ca, T, 'LineWidth', 2);
    hold on;


    xlabel('Crank Angle (degrees)');
    ylabel('Temperature (K)');
    title('Temperature vs Crank Angle');
    grid on;

    %% Log-Log Plot of Volume vs Pressure
    figure(2);
    plot(log(V(2:end)), log(p(2:end)), 'LineWidth', 2);
    xlabel('log Volume (m^3)');
    ylabel('log Pressure (Pa)');
    title('Log-Log P-V Curve');
    grid on;
end
