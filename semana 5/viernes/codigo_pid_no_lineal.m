clc
clear
close all
syms x1 x2 x3 u1 u2 u3 a b c
f1 = -a*x1^2+x2+u2*u3;
f2 = -x2+b*x3-u1*u3;
f3 = -c*x3^2+exp(-u1);
F = [f1; f2; f3]
X = [x1; x2; x3]
U = [u1; u2; u3]
A = jacobian(F,X)
B = jacobian(F,U)
u10 = 1
u20 = 2
u30 = 3
x30 = sqrt((exp(-u10))/c)
x20 = b*x30 - u10*u30
x10 = sqrt((x20+u20*u30)/a)
Ap = subs(A,{x1,x2,x3,u1,u2,u3},{x10,x20,x30,u10,u20,u30})
Bp = subs(B,{x1,x2,x3,u1,u2,u3},{x10,x20,x30,u10,u20,u30})
%% Modelo incertidumbres
a = 1
b = 1
c = 1
x30n = sqrt((exp(-u10))/c)
x20n = b*x30n - u10*u30
x10n = sqrt((x20n+u20*u30)/a)
a11 = -2*a*((b*(828390857088487/(2251799813685248*c))^(1/2) + 3)/a)^(1/2)
a33 = -2*c*(828390857088487/(2251799813685248*c))^(1/2)
Tol = 10
a11n = ureal('a11n',a11,'Percentage',[-Tol Tol])
a33n = ureal('a33n',a33,'Percentage',[-Tol Tol])
bn = ureal('bn',b,'Percentage',[-Tol Tol])
Ac = [a11n 1 0;
     0, -1, bn;
     0, 0, a33n]
BL = double(Bp)
CL = [2 1 0;
      0 2 4]
DL = [0 0 0;
      0 0 0]
sys_inc = ss(Ac,BL,CL,DL)
sys = sys_inc.NominalValue
AL = sys.a
%% An�lisis de ganancias relativas
Gf = tf(sys)
K = dcgain(Gf)
RGA = K.*(pinv(K))'
%% Dise�o de control PID
r = 2
s = tf('s')
G12 = Gf(1,2)
C12 = 0.26*(s+3.29)/s % Mp = 0; tr = 1.66 sg; ts = 3 sg
G21 = Gf(2,1)
C21 = -0.15*(s+1)/s % Mp = 0; tr = 1.86 sg; ts = 3.2 sg
Gc = [0 C21;
      C12 0;
      0 0]
T = feedback(sys_inc*Gc,eye(r)) 
U = feedback(Gc,sys_inc)
figure(1)
step(T)
grid
figure(2)
step(U)
grid
%% Condiciones iniciales PID
xini = [0.2; 0.3; 0.4]
xo = [xini; zeros(2,1)]
figure(1)
initial(T,xo)
grid
figure(2)
initial(U,xo)
grid
%% PID 2 DOF
G12 = Gf(1,2)
C12a = 0.16*(s+5.44)/s % Mp = 0; tr = 1.56 sg; ts = 3.06 sg
C12b = 0.0048*(s+5.23)/(s+0.5)
G21 = Gf(2,1)
C21a = -1*(s+1)/s % Mp = 8; tr = 0.23 sg; ts = 3 sg
C21b = 0.09*(s+2)/(s+1)

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
figure(3)
step(T2)
grid
figure(4)
step(U2)
grid