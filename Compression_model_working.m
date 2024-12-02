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

% Diesel properties : For now we will use B7 diesel, which is the
% commercially available in Europe

rho = 836.1;                      % Density [kg/m^3]
cn = 52.2;                        % Cetane ratio, how ignitable it is [-]
eta = 2.7638e-6;                  % Viscosity [m2/kg]

%% Assumed data Ingore these for stoichiometric analyis except lhv
lhv = [43e6;37e6];                % Lower heating value, energy content per kg of fuel [J/kg] left pure diesel, right FAME
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

% Initialize pressure and temperature
p(1) = P0;
T(1) = Tref;


W_inst = zeros(1, 720);                                         % Instantaneous work done at each crank angle


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
        case (ca(n) > 180 && ca(n) <= 360)
              %Elements 
            iSp = myfind({Sp.Name},{'O2','N2','CO2','H2O','Diesel'});
            SpS = Sp(iSp);
            NSp = length(SpS);
            Mi = [SpS.Mass];

            R_O2 = Runiv/(Mi(1));
            R_N2 = Runiv/(Mi(2));
            R_Air = (0.79*R_O2)+(0.21*R_N2);
            
            CpO2 = CpNasa(T(n-1),SpS(1));
            CpN2 = CpNasa(T(n-1),SpS(2));

            CpAir = (0.79*CpO2)+(0.21*CpN2);
            gamma = CpAir/(CpAir-R_Air);
            
            
            T(n) = T(n-1) * (V(n-1) / V(n))^(gamma - 1);                                                 
            p(n) = p(n-1) * (V(n-1) / V(n))^gamma;

          

        %% Combustion Phase
        case (ca(n) > 360 && ca(n) <= 390)
            CpO2 = CpNasa(T(n-1),SpS(1));
            CpN2 = CpNasa(T(n-1),SpS(2));
            CpAir = (0.79*CpO2)+(0.21*CpN2);
            M_air = (0.79*Mi(1))+(0.21*Mi(2));

      

            M_Diesel = Mi(5);
            
            Mf_Air = M_air/(M_Diesel+M_air);
            Mf_Diesel = M_Diesel/(M_Diesel+M_air);

            CpDiesel = CpNasa(T(n-1),SpS(5));
            
            CpComb = (Mf_Air*CpAir) + (Mf_Diesel*CpDiesel); 

            dQ_comb = sum(lhv .* [md * 0.93; md * 0.07]);           % Heat released from 93% pure diesel and 7% FAME
            dQ(n) = 0;                                              % For easier calculations, assume that there is no heat exchange in this small time instance
            p(n) = p(n-1);
            T(n) = T(359) + (V(n)*(p(n-1)-p(n)) + dQ_comb - dQ(n)) / (mtot * CpComb);        % Based on the first law of thermodynamics
            
            
        %% Expansion Stroke
        case (ca(n) > 390 && ca(n) <= 540)
            
    
            M_diesel = Mi(5);

            R_CO2 = Runiv/(Mi(3));
            R_H2O= Runiv/(Mi(4));
             


            O2toDiesel = 71/4;
            CO2toDiesel = 48/4;
            H2OtoDiesel = 46/4;

            %reactants
            molDiesel = md/M_diesel;
            molO2 = molDiesel * O2toDiesel;
            %products
            molCO2 = molDiesel * CO2toDiesel;
            molH2O = molDiesel * H2OtoDiesel;

            molAir = ma/M_air;

            molO2air = 0.79*molAir;
            molN2 = 0.21 * molAir;

            molO2Final = molO2air - molO2;

            moltot = molCO2 + molH2O + molO2Final + molN2;

            mfO2 = molO2Final/moltot;
            mfH2O = molH2O/moltot;
            mfCO2 = molCO2/moltot;
            mfN2 = molN2/moltot;



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
eta = W_total/(sum(lhv.*[md*0.93;md*0.07]));

%% Plot Results


% Display results
fprintf('The mass fuel flow rate %.10f \n', m_fuel_dot);


fprintf('The total work done during the cycle is %.2f J\n', W_total);
fprintf('The efficiency of the cycle is %.2f \n', eta);
figure(3)
subplot(1, 2, 1)
plot(ca, W_inst, 'LineWidth', 2);
xlabel('Crank Angle (degrees)');
ylabel('Instantaneous Work (J)');
title('Instantaneous Work vs. Crank Angle');
grid on;


%% Plotting Results

% Create the P-V diagram
figure(1)
subplot(1, 2, 1)
grid on
plot(V(2:end), p(2:end), 'LineWidth', 2)
hold on;

% Annotate stroke phases on P-V diagram
text(V(1), p(1), 'Intake', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(V(180), p(180), 'Compression', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(V(359), p(359), 'Combustion', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);

% Mark transition points on P-V diagram (ca. 180, 360, 540)
plot(V(180), p(180), 'ro', 'MarkerFaceColor', 'r');
plot(V(360), p(360), 'go', 'MarkerFaceColor', 'g');
plot(V(540), p(540), 'bo', 'MarkerFaceColor', 'b');

grid off
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('P-V Diagram');
legend('P-V curve', 'Location', 'best');

% Create the Temperature vs Crank Angle plot
subplot(1, 2, 2)
plot(ca, T, 'LineWidth', 2)
hold on;

% Annotate phases with text labels
text(90, T(90), 'Intake', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 12);
text(270, T(270), 'Compression', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(360, T(360), 'Combustion', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);
text(450, T(450), 'Expansion', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);
text(630, T(630), 'Exhaust', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);

% Mark transition points
plot([180 180], [min(T) max(T)], 'r--'); % Intake to Compression
plot([360 360], [min(T) max(T)], 'g--'); % Compression to Combustion
plot([540 540], [min(T) max(T)], 'b--'); % Combustion to Expansion

xlabel('Crank Angle (degrees)');
ylabel('Temperature (K)');
title('Temperature vs Crank Angle');
grid on;

%% Log-Log Plot of Volume vs Pressure
figure(2)
plot(log(V(2:end)), log(p(2:end)), 'LineWidth', 2);
xlabel('log Volume (m^3)');
ylabel('log Pressure (Pa)');
title('Log-Log P-V Curve');
grid on;

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