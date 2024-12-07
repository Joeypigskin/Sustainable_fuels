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
B   = Cyl.Bore;
S   = Cyl.Stroke;
cr  = Cyl.CompressionRatio;
r   = S/2;
l   = Cyl.ConRod;
%-------------------------------------------------------------------------------------------------------
CAl     = Ca-Cyl.TDCangle;
Apiston = pi*(B/2)^2;
Vd      = pi*(B/2)^2*S;
Vc      = Vd/(cr-1);
V       = Vc + Apiston*(r+l-(r*cosd(Ca) + sqrt(l^2-r^2 * (sind(Ca).^2)))); % 'sind' is the sine function taking arguments in degrees instead of radians





