clc
clear
close all
A = [-0.005509 0 0 -0.1588;
    0 -0.2062 0 0;
    -0.01216 0 0 -0.5672;
    0 0 0 -0.040]
B = [0.28 0 -0.01348 0;
    -9.375 7.658 0 0;
    0 0 0.7317 0;
    0.02999 0 0 0.04]
C = [14.21 0 0 0;
    0 1 0 0;
    0.3221 0 0.1434 11.16;
    0.4133 0 0 19.28]
D = [0 0 0 0;
    0 0 0 0;
    1.272 0 -0.208 0;
    0 0 0 0]
sys = ss(A,B,C,D)
figure(1)
step(sys)
grid
Gf = tf(sys)
%% Analisis RGA
s = tf('s')
Kst = dcgain(Gf)
Kv = dcgain(s*Gf)
K = [Kst(1:2,:); Kv(3,:); Kst(4,:)]
RGA = K.*(inv(K))'
%% Diseño PID
G14 = Gf(1,4)
C14 = -0.0045*(s+0.09)*(s+0.0055)/(s*(s+0.1))
G22 = Gf(2,2)
C22 = 0.00016*(s+7)*(s+0.7)/(s*(s+2.32))
G33 = Gf(3,3)
C33 = 0.1/(s+1)
G41 = Gf(4,1)
C41 = 0.021615*(s+0.05502)/s
%% Respuesta transitoria de la variable
T = feedback(G33*C33,1)
U = feedback(C33,G33)
figure(1)
step(T)
grid
figure(2)
step(U)
grid
%% Análisis de acoplamiento
r = 4 % Número de salidas
Gc = [0 0 0 C41;
      0 C22 0 0;
      0 0 C33 0;
      C14 0 0 0]
T = feedback(sys*Gc,eye(r))
U = feedback(Gc,sys)
figure(1)
step(T)
grid
figure(2)
step(U)
grid
%% Pto Equilibrio
xini = [0.2; 0.3; -0.4; -0.1]
xo = [xini; zeros(6,1)]
figure(1)
initial(T,xo)
grid
figure(2)
initial(U,xo)
grid
%% Diseño PID 2 DOF
G14 = Gf(1,4)
C14a = -0.0075*(s+0.00519)/s
C14b = -0.01*s/(s+0.1)
G22 = Gf(2,2)
C22a = 0.0125/s 
C22b = 5*(s+1)/(s+5)
G33 = Gf(3,3)
C33a = 0.0045/s
C33b = 0.23*(s+1.2)/(s+0.6)
G41 = Gf(4,1)
C41a = 0.001/s
C41b = 0.005*(s+1)/(s+5)
%% Respuesta transitoria de la variable
T = feedback(G14*C14a,1+C14b/C14a)
U = feedback(C14a,G14*(1+C14b/C14a))
figure(1)
step(T)
grid
figure(2)
step(U)
grid