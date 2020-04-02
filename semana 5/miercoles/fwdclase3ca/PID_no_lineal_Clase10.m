clc
clear
close all
syms x1 x2 x3 u1 u2 u3 a b c
f1 = -a*x1^2 + x2 +u2*u3
f2 = -x2 + b*x3 - u1*u3
f3 = -c*x3^2 + exp(-u1)
F = [f1; f2; f3]
X = [x1; x2; x3]
U = [u1;u2;u3]
A = jacobian(F,X)
B = jacobian(F,U)
u10 = 1
u20 = 2
u30 = 3
x30 = sqrt(exp(-u10)/c)
x20 = b*x30 - u10*u30
x10 = sqrt((x20+u20*u30)/a)
Ap = subs(A,{x1,x2,x3,u1,u2,u3},{x10,x20,x30,u10,u20,u30})
Bp = subs(B,{x1,x2,x3,u1,u2,u3},{x10,x20,x30,u10,u20,u30})
%% Modelo incetidumbres
a = 1
b = 1
c = 1
x30n = sqrt(exp(-u10)/c)
x20n = b*x30n - u10*u30
x10n = sqrt((x20n+u20*u30)/a)

a11 = -2*a*((b*(828390857088487/(2251799813685248*c))^(1/2) + 3)/a)^(1/2)
a33 = -2*c*(828390857088487/(2251799813685248*c))^(1/2)
Tol = 5
a11n = ureal('a11n',a11,'Percentage',[-Tol,Tol])
a33n = ureal('a33n',a33,'Percentage',[-Tol,Tol])
bn = ureal('bn',b,'Percentage',[-Tol,Tol])
Ainc = [a11n, 1, 0;
        0, -1, bn;
        0, 0 , a33n]
AL = Ainc.NominalValue   
BL = [0 3 2;
    -3 0 -1;
    -exp(-1) 0 0]
CL = [2 1 0;
      0 2 4]
DL = zeros(2,3)
sys_inc = ss(Ainc,BL,CL,DL)
sys = ss(AL,BL,CL,DL)
eig(AL)
figure(1)
step(sys)
grid
%% PID 1 DOF
Gp = tf(sys)
K = dcgain(Gp)
RGA = K.*(pinv(K))'

s = tf('s')
G12 = Gp(1,2)
C12 = 0.17*(s+5.8)/s % Mp = 0%, tr = 1.28 sg, ts = 2.14 sg 
G21 = Gp(2,1)
C21 = -0.34*(s+2.9)/s % Mp = 14.3%, tr = 0.436 sg, ts = 1.85 sg 
%% Analisis de acoplamiento
Gc = [0 C21;
      C12 0;
      0 0]
T = feedback(sys*Gc,eye(2))
U = feedback(Gc,sys)
figure(1)
step(T)
grid
figure(2)
step(U)
grid
%%
xini = [0.2; -0.3; 0.5]
xo = [xini; zeros(2,1)]
figure(1)
initial(T,xo)
grid

Gz = c2d (sys,Tm)
figure(2)
initial(U,xo)
grid