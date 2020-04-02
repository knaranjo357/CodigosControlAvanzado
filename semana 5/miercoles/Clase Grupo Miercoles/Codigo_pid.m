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
Gf = tf(sys)
eig(Gf)
s = tf('s')
Kst = dcgain(Gf)
Kv = dcgain(s*Gf)
K = [Kst(1:2,:);
    Kv(3,:)
    Kst(4,:)]
RGA = K.*(inv(K))'
figure(1)
step(sys)
grid
%% Diseño PID
G14 = Gf(1,4)
C14 = -0.0066*(s+0.0056)/s
G22 = Gf(2,2)
C22 = 0.000238*(s+4)/s
G33 = Gf(3,3)
C33 = 0.025/(s+0.125)  %  Mp = 0, tr = 82.3 sg, ts = 140 sg
G41 = Gf(4,1)
C41 = 0.06*(s+0.06159)/s
%%
T = feedback(G33*C33,1)
U = feedback(C33,G33)
figure(1)
step(T)
grid
figure(2)
step(U)
grid
%% Análisis de acoplamiento
r = 4 % Número de salidas del sistema
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
%% Análisis pto equilibrio
xini = [0.2; -0.6; -0.1; 0.6]
xo = [xini; zeros(4,1)]
figure(1)
initial(T,xo)
grid
figure(2)
initial(U,xo)
grid
%% PID 2 DOF
G14 = Gf(1,4)
C14a = -0.125*(s+0.0052)/s
C14b = -521.35*(s+0.004076)/(s+40)
G22 = Gf(2,2)
C22a = 0.0014*(s+2.28)/s
C22b = 0.063*(s+3.31)/(s+1.47)
G33 = Gf(3,3)
C33a = 0.01/s
C33b = 0.18*(s+0.7)/(s+0.18)
G41 = Gf(4,1)
C41a = 0.09*(s+2.211)/s
C41b = 1.81*(s+2.727)/(s+0.5)
%%
xini = [0.2; -0.6; -0.1; 0.6]
T = ss(feedback(G33*C33a,1+C33b/C33a))
U = ss(feedback(C33a,G33*(1+C33b/C33a)))

figure(1)
initial(T,xini)
grid
figure(2)
initial(U,xini)
grid