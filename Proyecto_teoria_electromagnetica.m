# Teoria electromagnetica
# Seccion 30
# Integrantes: Alejandro Duarte 19446, Marco Duarte 19022, Jonathan Pu 19249
# Descripcion: proyecto de representacion grafica de campo electrico y potencial
# de dipolos y medios polarizados

##% DIPOLOS ELECTRICOS Y EXPANSION MULTIPOLAR
##% Primer inciso
##clear; close all; clc;
##N = 20;
##minTheta = -2*pi; maxTheta = 2*pi;
##minPhi = -2*pi; maxPhi = 2*pi;
##minR = 2; maxR = 2;
##
###Vectores para dar un rango a cada variable
##Theta = [-2*pi -pi 0 pi 2*pi];
##R = [-2 -1 0 1 2];
##Phi = [-2*pi -pi 0 pi 2*pi];
##
###Definicion de ejes
### Theta_axis = linspace(minTheta, maxTheta, N);
### Phi_axis = linspace(minPhi, maxPhi, N);
### R_axis = linspace(minR, maxR, N);
##
##R_axis = R.*sin(Theta).*cos(Phi);
##Phi_axis = R.*sin(Theta).*sin(Phi);
##Theta_axis = R.*cos(Theta);
##
##[RG, PhiG, ThetaG] = meshgrid(1:4,1:4,1:4);
##
##
##
###Componentes del campo electrico
##eps0 = 8.85E-12;
##po = 100; #Un valor del momento dipolar para simular las graficas
##k = po/(4*pi*eps0);
##Er = k.*(-2.*sin(ThetaG).*cos(PhiG)-2.*sin(ThetaG).*sin(PhiG)-2.*cos(ThetaG))./(RG.^5);
##Ephi = k.*(sin(PhiG)-cos(PhiG))/(RG.^3); 
##Etheta = k.*(-cos(ThetaG).*cos(PhiG)-cos(ThetaG).*sin(PhiG)+sin(ThetaG))./(RG.^5);
##
##E = sqrt(Er.^2 + Ephi.^2 + Etheta.^2);
##
##u = Er./E;
##v = Ephi./E;
##w = Etheta./E;
##
##figure();
##h = quiver3(RG, PhiG, ThetaG, u, v, w);
##
###Potencial
##
##[pot1, pot2, pot3] = meshgrid(1:4,1:4,1:4);
##Pot = po.*((sin(pot3)).*cos(pot2) + sin(pot3).*sin(pot2) + cos(pot3))./(4*pi.*pot1.^2);
##figure();
##i = quiver3(pot1, pot2, pot3, 1:64, 1:64, 1:64);
##
##
##
##%POLARIZACION
###Problema 1
###Densidades de carga de polarizacion
##
###Definicion de variables
##
##k = 5;
##alpha0 = 10;
##
###Graficar polarizacion
##x=linspace(-1,1,10);
##y=linspace(-2*pi, 2*pi, 10);
##[ss, pp] = meshgrid(x, y);
##u = k.*ss.*e.^((ss.^2)/alpha0).*cos(pp);
##v = k.*ss.*e.^((ss.^2)/alpha0).*sin(pp);
##quiver(ss, pp, u, v)
##axis equal
##j = figure();
##
##xlabel('Eje x')
##ylabel('Eje y')
####s = linspace(1,10,1);
####Pphi = 0:2*pi;
####[ss, pp] = meshgrid(s, Pphi);
######xx = ss.*cos(tt);
######yy = ss.*sin(tt);
####componente_i = k.*e.^((-ss.^2)/alpha0).*cos(pp);
####componente_j = k.*e.^((-ss.^2)/alpha0).*sin(pp);
####
####normalizar = sqrt(componente_i.^2 + componente_j.^2);
####comp_i = componente_i./normalizar;
####comp_j = componente_j./normalizar;
####quiver(ss, pp, comp_i, comp_j)
####axis equal
##
####Pol = k.*e.^((s.^2)/alpha0);
####ThePol = 2*pi;
####mesh(Pol, ThePol, 0*Pol)
####figure();
###volumetrica = (-k.*(e.^(-x.^2)/alpha0))*((1./x-(2.*x/alpha0)));


clear; close all; clc;
% POTENCIAL APROXIMADO DE LA DENSIDAD DE CARGA
r = [-2*pi:pi/10:2*pi];
theta = [-2*pi:pi/10:2*pi];
[Mx, My] = meshgrid(r, theta);
Mz = 89991804694*(0.25-((Mx.^2+My.^2)/120).*((3*(cos(atan(My./Mx)).^2)-1)/2));
mesh(Mx, My, Mz)
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('POTENCIAL (EXPANSION MULTIPOLAR)')
