function [V] = CylinderVolume(Ca,Cyl)
% This function provides the cylinder volume as function of 
% Ca : Crankangle [degrees]
% Cyl :  a struct containing
%   Cyl.S : Stroke
%   Cyl.B                   : Bore
%   Cyl.ConRod              : Connecting Rod length
%   Cyl.CompressionRatio    : Compession Ratio
%   Cyl.TDCangle            : Angle associated with the Top Dead Center
%----------------------------------------------------------------------
fprintf('WARNING------------------------------------------------------------------\n');
fprintf(' Modify this function to yours. Now it is just a sinusoidal expression\n');
fprintf(' This function is %s\n',mfilename('fullpath'));
fprintf('END OF WARNING ----------------------------------------------------------\n');
B   = Cyl.Bore;
S   = Cyl.Stroke;
cr  = Cyl.CompressionRatio;
r   = S/2;
l   = Cyl.ConRod;
d   = Cyl.Displacement;
%-------------------------------------------------------------------------------------------------------
Ca = deg2rad(Ca);

%stroke = 2 * l * (1 - cos(Ca)); 
%r = l + S / 2;  % Piston's radius of travel

% Clearance volume based on compression ratio
Vc = d / (cr - 1);
% Volume calculation at a given crank angle
%V = Vc + (pi / 4) * bore^2 * (r - (r * cos(ca) + sqrt(conrodLength^2 - r^2 * sin(ca)^2)));
V = Vc + (pi / 4) * B.^2 .* (r *(1- cos(Ca) + l*sqrt(1 - (r/l).^2 .* sin(Ca).^2)));





