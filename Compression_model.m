clc, clear, close all
%% Initial Parameters
bore = 0.104;                     % Diameter of piston [m]
stroke = 0.085;                   % Full travel of piston [m]
conrodLength = 0.1365;            % Length of rod which connects the piston head and crankshaft [m]
cRatio = 21.5;                    % Ratio between min and max pressure[-]
displacement = 722e-6;            % The volume of the cylinder [m^3]
P0 = 101325;                      % Ambient Pressure [Pa]
 
global Runiv Pref Tref
Runiv= 8.314472;                   % Universal gas constant
Pref= 1.01235e5;                   % Reference pressure, 1 atm!
Tref= 298.15;                      % Reference Temperature

% Diesel properties : For now we will use B7 diesel, which is the
% commercially available in Europe

rho = 836.1;                      % Density [kg/m^3]
cn = 52.2;                        % Cetane ratio, how ignitable it is [-]
lhv = 43e6;                         % Lower heating value, energy content per kg of fuel [J/kg]
eta = 2.7638e-6;                  % Viscosity [m2/kg]
gamma = 1.4;                      % At this point we assume it to be 1.4, for ease of use

RPMn = 1000;                      % Unloaded RPM of a typical diesel engine
rotTime = 60/RPMn;                % Time of a full rotation of 360 degrees [s]
R_specific = 37;                  % This is an arbitrary value which is relative to the specific constant of B7 diesel, this will need to be found through stoichiometric relations
m = 0.00009;                      % Also and arbitrary value for now
Cp = 2.1e3;
% Crank angle and time setup
for i=1:1:720
    t(i) = i*(rotTime)/720;
    ca(i) = i;  % Crank angle in degrees, from 1 to 720
end

% Initialize pressure and temperature
p(1) = P0;
T(1) = Tref;

% Loop over each crank angle to compute volume, pressure, and temperature
for n=2:1:720
    % Pass a single crank angle ca(n) to the Volume function
    V(n) = Volume(ca(n), cRatio, conrodLength, bore, displacement);  % Corrected here: Pass scalar ca(n)

    switch true
        %% Intake Stroke
        case (ca(n) <= 180)
            p(n) = p(1);                                    % Pressure [Pa]
            T(n) = T(1);                                    % Temperature [K]
        %% Compression Stroke
        case (ca(n) > 180 && ca(n) <= 359)
            dQ_comb(n)=0;                                   % Heat released from combustion [J]
            dQ(n)=0;                                        % Heat extracted from cycle [J]
            T(n)=T(n-1)*(V(n-1)/V(n))^(gamma-1);            % Temperature after compression following first law [K]
            p(n)=p(n-1)*(V(n-1)/V(n))^gamma;                % Pressure during compression following Poisson relations [Pa]

        %% Combustion Phase
        case (ca(n) > 359 && ca(n) <= 360)
            dQ_comb(n) = lhv * m;    % Heat released from combustion [J] based on stoichiometric content
            dQ(n)=0;                                       % Heat extracted from cycle [J]
            T(n)=T(n-1)+(p(n-1)*(V(n-1)-V(n))+dQ_comb(n)-dQ(n))/(m*Cp);      % Temperature after combustion following first law [K]
            p(n)=(T(n)*R_specific*m)/(V(n));               % Pressure after combustion [Pa]

        %% Expansion Stroke
        case (ca(n) > 360 && ca(n) <= 540)
            T(n)=T(n-1)+(p(n-1)*(V(n-1)-V(n)))/(m*Cp);  % Temperature during expansion following first law [K]
            p(n)=p(n-1)*(V(n-1)/V(n))^gamma;               % Pressure during expansion following Poisson relations [Pa]

        %% Exhaust Phase
        case (ca(n) > 540 && ca(n) <= 720)
            p(n)=p(1);
            T(n)=T(1);                                    % Temperature of the mixture (K)

        otherwise
            % Handle unexpected values of ca
            disp('Crank angle out of range (0 to 720 degrees).');
    end
end

%% Plotting Results

% Create the P-V diagram
figure(1)
subplot(1,2,1)
grid on
plot(V(1,2:end), p(1,2:end), 'LineWidth', 2)
hold on;

% Annotate stroke phases on P-V diagram
text(V(1,180), p(1,180), 'Intake', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(V(1,360), p(1,360), 'Compression', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(V(1,540), p(1,540), 'Combustion', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 12);

% Mark transition points on P-V diagram (ca. 180, 360, 540)
plot(V(1,180), p(1,180), 'ro', 'MarkerFaceColor', 'r');
plot(V(1,360), p(1,360), 'go', 'MarkerFaceColor', 'g');
plot(V(1,540), p(1,540), 'bo', 'MarkerFaceColor', 'b');

grid off
xlabel('Volume (m^3)');
ylabel('Pressure (Pa)');
title('P-V Diagram');
legend('P-V curve', 'Location', 'best');

% Create the Temperature vs Crank Angle plot
subplot(1,2,2)
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
plot(log(V(1,2:end)), log(p(1,2:end)), 'LineWidth', 2);
xlabel('log Volume (m^3)');
ylabel('log Pressure (Pa)');
title('Log-Log P-V Curve');
grid on;

%% Volume Function
function V = Volume(ca, cRatio, conrodLength, bore, displacement)
    % Input parameters:
    % ca: Crank angle (in degrees)
    % cRatio: Compression ratio
    % conrodLength: Length of the connecting rod (m)
    % bore: Bore diameter of the engine (m)
    % displacement: Total displacement volume (m^3)

    % Convert crank angle to radians
    ca = deg2rad(ca);

    % Calculate stroke length based on the crank angle
    stroke = 2 * conrodLength * (1 - cos(ca));  % Stroke distance based on crank angle (m)
    r = conrodLength + stroke / 2;  % Piston's radius of travel

    % Clearance volume based on compression ratio
    Vc = displacement / (cRatio - 1);

    % Volume calculation at a given crank angle
    %V = Vc + (pi / 4) * bore^2 * (r - (r * cos(ca) + sqrt(conrodLength^2 - r^2 * sin(ca)^2)));
    V = Vc + (pi / 4) * bore^2 * (r *(1- cos(ca) + conrodLength*sqrt(1 - (r/conrodLength)^2 * sin(ca)^2)));
end