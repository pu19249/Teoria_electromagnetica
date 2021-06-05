# Teoria electromagnetica
# Seccion 30
# Integrantes: Alejandro Duarte 19446, Marco Duarte 19022, Jonathan Pu 19249
# Descripcion: proyecto de representacion grafica de campo electrico y potencial
# de dipolos y medios polarizados

% ------------------ DIPOLOS ELECTRICOS Y EXPANSION MULTIPOLAR -----------------
% Primer inciso

clear; close all; clc;
#Potencial de un dipolo orientado en los tres ejes
figure
hold on
theta_0 = pi;
phi_0 = 2.*pi;
r_0= 1:1:10;
#Asumiendo que Po = 0.001
pot_0 = 8991804.694*((sin(theta_0)*cos(phi_0))+(sin(theta_0)*sin(phi_0)+cos(theta_0)./(r_0.^2)));
plot(r_0, pot_0, '-.og')
xlabel('r')
ylabel('V(r)')

title('POTENCIAL ORIENTADO EN LOS TRES EJES')

#Potencial de dipolo orientado en x
figure
hold on
theta2 = pi;
phi2 = 2.*pi;
r2 = 1:1:10;
#Asumiendo que Po = 0.001
pot_x = 8991804.694*((sin(theta2))*cos(phi2))./(r2.^2);
plot(r2, pot_x, '-.og')
xlabel('r')
ylabel('V(r)')

title('POTENCIAL ORIENTADO EN X')

#Potencial de dipolo orientado en y
figure
theta3 = pi;
phi3 = 2.*pi;
r3 = 1:1:10;
#Asumiendo que Po = 0.001
pot_y = 8991804.694*((sin(theta3))*sin(phi3))./(r3.^2);
plot(r3, pot_y, '-.ob')
xlabel('r')
ylabel('V(r)')

title('POTENCIAL ORIENTADO EN Y')

#Potencial de dipolo orientado en z
figure
theta4 = pi;
phi4 = 2.*pi;
r4 = 1:1:10;
#Asumiendo que Po = 0.001
pot_z = 8991804.694*((cos(theta4))./(r4.^2));
plot(r4, pot_z, '-.om')
xlabel('r')
ylabel('V(r)')

title('POTENCIAL ORIENTADO EN Z')

#Potencial de un cuadripolo paralelo al eje z
figure
a = 5;        #La separacion entre sus centros
q = 1e-9;     #Se asume una carga de 1nC
d = 6;        #La separcion entre las cargas
R = 0:1:15;   #Si se mide desde el centro hasta una distancia R de 15
alpha1 = asin(((d/2)^2-R.^2-(a^2)/4)./(a*R.^2)); #Varian junto con el radio de
alpha2 = asin((R+R.^2+(a^2)/4)./a.*R);           #donde se miden



pot_qua = (q*d/(1.11e-10))*(((cos(alpha1))./(R.^2+(a^2)/4+a*sin(alpha1).*R))-((cos(alpha2))./(R.^2+(a^2)/4-a*sin(alpha2).*R)));
plot(R, pot_qua, '-.ob')
xlabel('R')
ylabel('V(R)')

title('POTENCIAL DOS DIPOLOS PARALELOS')

% POTENCIAL APROXIMADO DE LA DENSIDAD DE CARGA
figure
r = [-2*pi:pi/10:2*pi];
theta = [-2*pi:pi/10:2*pi];
[Mx, My] = meshgrid(r, theta);
Mz = 89991804694*(0.25-((Mx.^2+My.^2)/120).*((3*(cos(atan(My./Mx)).^2)-1)/2));
mesh(Mx, My, Mz)
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('POTENCIAL (EXPANSION MULTIPOLAR)')

figure
r1 = [-2*pi:pi/10:2*pi];
theta1 = [-2*pi:pi/10:2*pi];
[Mx1, My1] = meshgrid(r1, theta1);
Mz = 89991804694*(sin(My1)./(r1.^2));
mesh(Mx, My, Mz)
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('POTENCIAL ORIENTADO EN X')


% ------------------------------ POLARIZACION ----------------------------------

#Densidades de carga
figure
a = 5;  #radio interno a
b = 10; #radio externo b
k = 2;
alpha0 = 0.01;
d_radio = linspace(a, b, 11);
theta5 = linspace(0, 2*pi, 11);
[Mr, Mt] = meshgrid(d_radio, theta5);
xr = Mr.*cos(Mt); #parametrizo mis valores
yr = Mr.*sin(Mt);
densidad = k*exp(-xr.^(2/alpha0))*(1./xr-2*xr/alpha0);
surf(xr, yr, densidad, Mr)
colorbar
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('DENSIDAD DE CARGA VOLUMETRICA')
densidad_a = -k*exp(-xr.^(2/alpha0));
densidad_b = k*exp(-xr.^(2/alpha0));
figure
hold on
surf(xr, yr, densidad_a, Mr)
surf(xr, yr, densidad_b, Mr)
colorbar
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('DENSIDAD DE CARGA SUPERFICIAL')
hold off

#Campo electrico
figure
for x = 0:b #Desde el centro hasta mi radio externo
  x = linspace(-10, b, 11);
  y = x';   #La derivada de x
  u = -2.259*10^10.*exp(-(x.^2+y.^2)/0.01);
  [ex, ey] = gradient(u);
  
endfor
hold on
quiver(x, y, ex, ey)
contour(x, y, u)
hold off
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('CAMPO ELECTRICO')

#Densidades de las esferas
#a<r<b
k_esf = 2;
Q = 1e-9;
a_esf = 5;
b_esf = 10;
e0 = 8.85e-12;
er = 2.24; #dielectrico
xe = 1.24; #er=xe-1
radio_esf = linspace(a_esf, b_esf, 11);
theta_esf = linspace(0, 2*pi, 11);
[Mesf1, Mesft1] = meshgrid(radio_esf, theta5);
xesf1 = Mesf1.*cos(Mesft1); #parametrizo mis valores
yesf1 = Mesf1.*sin(Mesft1);
densidad_esf1 = -((e0*xe)+Q)./(4*pi*er*e0*xesf1.^2); #superficial en a
densidad_esf1_1 = ((e0*xe)+Q)./(4*pi*er*e0*xesf1.^2); #superficial en b
figure
hold on
surf(xesf1, yesf1, densidad_esf1, Mesf1)
surf(xesf1, yesf1, densidad_esf1_1, Mesf1)
hold off
colorbar
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('DENSIDAD DE CARGA VOLUMETRICA ESFERA a<r<b')
figure
hold on
plot(xesf1, densidad_esf1)
plot(xesf1, densidad_esf1_1)
hold off
colorbar
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('DENSIDAD DE CARGA SUPERFICIAL ESFERA a<r<b')

#a<r<b
radio_esf2 = linspace(b_esf, 2*b_esf, 11);
theta_esf2 = linspace(0, 2*pi, 11);
[Mesf2, Mesft2] = meshgrid(radio_esf2, theta_esf2);
xesf2 = Mesf2.*cos(Mesft2); #parametrizo mis valores
yesf2 = Mesf2.*sin(Mesft2);
densidad_esf2 = ((xesf2.^2)*k_esf); #superficial en b
densidad_esf2_1 = -4*(((xesf2.^2)*k)); #superficial en 2b
densidad_esf2_2 = 4*k_esf*xesf2.^2; #volumetrica
figure
surf(xesf1, yesf1, densidad_esf2_2, Mesf1)
colorbar
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('DENSIDAD DE CARGA VOLUMETRICA ESFERA b<r<2b')
figure
hold on
plot(xesf1, densidad_esf2)
plot(xesf1, densidad_esf2_1)
hold off
colorbar
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('DENSIDAD DE CARGA SUPERFICIAL ESFERA b<r<2b')

#Grafica de Desplazamiento vs distancia
distancia1 = 0:1:a_esf;
D1 = 0./4*pi*distancia1.^2;
distancia2 = a_esf:1:2*b_esf;
D2 = Q./4*pi*distancia2.^2;
figure
hold on
plot(distancia1, D1, '-.or')
plot(distancia2, D2, '-.ob')
hold off
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('DESPLAZAMIENTO ELECTRICO VS DISTANCIA')

#Grafica de Campo electrico vs distancia
distancia3 = linspace(0, a_esf, 5);
distancia4 = linspace(a_esf, b_esf, 5);
distancia5 = linspace(b_esf, 2*b_esf, 5);
distancia6 = linspace(2*b_esf, 30, 5);  #Al infinito pero se acorta para la grafica y ver su comporatmiento

E_1 = 0./4*pi*er*e0*distancia3.^2; #r<a
E_2 = Q./4*pi*er*e0*distancia4.^2; #a<r<b
E_3 = (Q./4*pi*er*e0*distancia5.^2)+(k_esf/e0)*distancia5.^2; #b<r<2b
E_4 = Q./4*pi*e0*distancia6.^2; #r>2b

figure

plot(distancia3, E_1, '-.or')
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('CAMPO ELECTRICO VS DISTANCIA r<a')

figure
plot(distancia4, E_2, '-.ob')
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('CAMPO ELECTRICO VS DISTANCIA a<r<b')

figure
plot(distancia5, E_3, '-.og')
xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('CAMPO ELECTRICO VS DISTANCIA b<r<2b')

figure
plot(distancia6, E_4, '-.om')

xlabel('EJE X')
ylabel('EJE Y')
zlabel('EJE Z')
title('CAMPO ELECTRICO VS DISTANCIA r>2b')



