clc
clear
close all
syms x1 x2 x3 u1 u2 u3 a b c
f1 = -a*x1^2 + x2 + u2*u3;
f2 = -x2 + b*x3 - u1*u3;
f3 = -c*x3^2 + exp(-u1);
F = [f1; f2; f3];
X = [x1; x2; x3];
U = [u1; u2; u3];
A = jacobian(F,X)
B = jacobian(F,U)
u10 = 1;
u20 = 2;
u30 = 3;
x30 = sqrt((exp(-u10))/c)
x20 = b*x30 - u10*u30
x10 = sqrt((x20 + u20*u30)/a)
Ap = subs(A,{x1,x2,x3,u1,u2,u3},{x10,x20,x30,u10,u20,u30})
Bp = subs(B,{x1,x2,x3,u1,u2,u3},{x10,x20,x30,u10,u20,u30})
CL = [2 1 0; 0 2 4]
DL = [0 0 0;
      0 0 0]
%% Modelo incertidumbres
a = 1;
b = 1;
c = 1;
x30n = sqrt((exp(-u10))/c)
x20n = b*x30n - u10*u30
x10n = sqrt((x20n + u20*u30)/a)

a11 = -2*a*((b*(828390857088487/(2251799813685248*c))^(1/2) + 3)/a)^(1/2)
a33 = -2*c*(828390857088487/(2251799813685248*c))^(1/2)
Tol = 10;
a11n = ureal('a11n',a11,'Percentage',[-Tol Tol])
a33n = ureal('a33n',a33,'Percentage',[-Tol Tol])
bn = ureal('bn',b,'Percentage',[-Tol Tol])
Ac = [a11n,1, 0;
      0, -1, bn;
      0, 0, a33n]
BL = double(Bp)
sys_inc = ss(Ac,BL,CL,DL)
sys = sys_inc.NominalValue
AL = sys.a
%% Diseño PID
s = tf('s')
Gf = tf(sys)
K = dcgain(Gf)
RGA = K.*(pinv(K))'
%% 
G12 = Gf(1,2)
C12 = 0.327*(s+3.71)/s % Mp = 0%, tr = 2 sg , ts = 2 sg
G21 = Gf(2,1)
C21 = -0.44947*(s+1.9)/s  % Mp = 7.4 %, tr = 0.442 sg , ts = 1.94 sg

Gc = [0 C21;
      C12 0;
      0 0]
T = feedback(sys_inc*Gc,eye(2))
U = feedback(Gc,sys_inc)
figure(1)
step(T)
grid
figure(2)
step(U)
grid
%% Condiciones iniciales
xini = [0.1; 0.2; 0.6]
xo = [xini; zeros(2,1)]
figure(1)
initial(T,xo)
grid
figure(2)
initial(U,xo)
grid
%% PID 2 DOF
G12 = Gf(1,2)
C12a = 0.3177*(s+73.9)/s % Mp = 1.44%, tr = 1.58 sg , ts = 2.37 sg
C12b = 344.6*(s+2.404)/(s+36.9)
G21 = Gf(2,1)
C21a = -0.43366/s  % Mp = 1.24 %, tr = 1.29 sg , ts = 1.98 sg
C21b = -0.20261*(s+4.127)/(s+3.428)
%% Analisis de incertidumbres PID 2 DOF
G12inc = ss(Ac,BL(:,2),CL(1,:),DL(1,2))
G21inc = ss(Ac,BL(:,1),CL(2,:),DL(2,1))
T1 = feedback(C12a*G12inc,1+C12b/C12a)
U1 = feedback(C12a,G12inc*(1+C12b/C12a))
figure(1)
step(T1)
grid
figure(2)
step(U1)
grid

T2 = feedback(C21a*G21inc,1+C21b/C21a)
U2 = feedback(C21a,G21inc*(1+C21b/C21a))
figure(2)
step(T2)
grid
figure(2)
step(U2)
grid
