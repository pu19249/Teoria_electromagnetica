# Teoria electromagnetica
# Seccion 30
# Integrantes: Alejandro Duarte 19446, Marco Duarte 19022, Jonathan Pu 19249
# Descripcion: proyecto de representacion grafica de campo electrico y potencial
# de dipolos y medios polarizados

% DIPOLOS ELECTRICOS Y EXPANSION MULTIPOLAR
% Primer inciso
clear; close all; clc;
N = 20;
minTheta = -2*pi; maxTheta = 2*pi;
minPhi = -2*pi; maxPhi = 2*pi;
minR = 2; maxR = 2;

#Vectores para dar un rango a cada variable
Theta = [-2*pi -pi 0 pi 2*pi];
R = [-2 -1 0 1 2];
Phi = [-2*pi -pi 0 pi 2*pi];

#Definicion de ejes
# Theta_axis = linspace(minTheta, maxTheta, N);
# Phi_axis = linspace(minPhi, maxPhi, N);
# R_axis = linspace(minR, maxR, N);

R_axis = R.*sin(Theta).*cos(Phi);
Phi_axis = R.*sin(Theta).*sin(Phi);
Theta_axis = R.*cos(Theta);

[RG, PhiG, ThetaG] = meshgrid(1:4,1:4,1:4);



#Componentes del campo electrico
eps0 = 8.85E-12;
po = 100; #Un valor del momento dipolar para simular las graficas
k = po/(4*pi*eps0);
Er = k.*(-2.*sin(ThetaG).*cos(PhiG)-2.*sin(ThetaG).*sin(PhiG)-2.*cos(ThetaG))./(RG.^5);
Ephi = k.*(sin(PhiG)-cos(PhiG))/(RG.^3); 
Etheta = k.*(-cos(ThetaG).*cos(PhiG)-cos(ThetaG).*sin(PhiG)+sin(ThetaG))./(RG.^5);

E = sqrt(Er.^2 + Ephi.^2 + Etheta.^2);

u = Er./E;
v = Ephi./E;
w = Etheta./E;

figure();
h = quiver3(RG, PhiG, ThetaG, u, v, w);
