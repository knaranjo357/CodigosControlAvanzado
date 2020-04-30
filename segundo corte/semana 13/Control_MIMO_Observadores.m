clc
clear
close all
syms x1 x2 x3 u1 u2 u3 a b c d
f1 = -a*x2 - b*x1^2 + u1*u2;
f2 = -c*x3*x2 + u2^2 ;
f3 = -d*x3 + u3*x1 + u1;

F = [f1; f2; f3]
X = [x1; x2; x3]
U = [u1; u2; u3]
A = jacobian(F,X)
B = jacobian(F,U)
x10 = 2;
x20 = 4;
x30 = 4;
u20 = sqrt(c*x30*x20)
u10 = (a*x20 + b*x10^2)/u20
u30 = (d*x30 - u10)/x10

AL = subs(A,{x1,x2,x3,u1,u2,u3},{x10,x20,x30,u10,u20,u30})
BL = subs(B,{x1,x2,x3,u1,u2,u3},{x10,x20,x30,u10,u20,u30})
%%
a1 = 2
b1 = 4
c1 = 1
d1 = 3
x10 = 2;
x20 = 4;
x30 = 4;
u20 = double(sqrt(c1*x30*x20))
u10 = double((a1*x20 + b1*x10^2)/u20);
u30 = double((d1*x30 - u10)/x10)
ALn = double(subs(AL,{a,b,c,d},{a1,b1,c1,d1}))
BLn = double(subs(BL,{a,b,c,d},{a1,b1,c1,d1}))
CLn = [1 0 0;
       0 0 1]
DLn = [0 0 0;
       0 0 0]
yo = CLn*[x10;x20;x30] + DLn*[u10;u20;u30] % Pto Equilibrio vector salida
sys =ss(ALn,BLn,CLn,DLn)
eig(ALn)
Gs = zpk(sys)
K = dcgain(Gs)
RGA = K.*pinv(K)'
%% Matriz de controlabilidad
n = 3 % Numero de estados
p = 3 % N�mero de actuadores
r = 2 % N�mero de salidas
Ahat = [ALn zeros(n,r);
        -CLn zeros(r,r)];
Bhat = [BLn;
        -DLn]    
M = ctrb(Ahat,Bhat)
rank(M)
%% Control LQR
Q = diag([1 1 1 30 10]);
R = diag([0.1 0.1 0.1])
[Kt,S,E] = lqi(sys,Q,R)
Kest = Kt(:,1:n)
Ki = -Kt(:,n+1:end);
%% Modelo lineal en lazo cerrado para seguimiento
AA = [ALn - BLn*Kest, BLn*Ki;
     -CLn + DLn*Kest, -DLn*Ki]
eig(AA)
BB = [zeros(n,r);
      eye(r)];
CC = [CLn - DLn*Kest, DLn*Ki]
DD = zeros(r,r);
sys_new = ss(AA,BB,CC,DD)
figure(1)
step(sys_new)
grid

% Acci�n de control
CC1 = -Kt;
DD1 = zeros(p,r);
sys_u = ss(AA,BB,CC1,DD1)
figure(2)
step(sys_u)
grid
%% Obsevabor de Luenberger polos reales
Vob = obsv(ALn,CLn)
rank(Vob) % rango debe ser igual al numero de estados
ts = 1
s1 = -4/ts
delta = diag([s1 3*s1 4*s1]) % Matriz cuadrada igual al numero de estados n*n
% Grand = 2*rand(r,n)-1;  % Matriz de r*n
Grand = [-0.443003562265903 0.915013670868595 -0.684773836644903;0.0937630384099680 0.929777070398553 0.941185563521231];
X = lyap(ALn',-delta,-CLn'*Grand)
det(X)
L = (Grand*inv(X))'
eig(ALn-L*CLn) 
%% Obsevabor de Luenberger polos complejos
ts = 0.3
Mp = 0.1
zita = abs(log(Mp))/sqrt(pi^2+(log(Mp))^2)
wn = 4/(ts*zita)
wd = wn*sqrt(1-zita^2)
s1 = -zita*wn + wd*j
delta = [real(s1) imag(s1) 0;
         -imag(s1) real(s1) 0;
          0 0 2*real(s1)]
% Grand = 2*rand(r,n)-1;  % Matriz de r*n
Grand = [0.389657245951634 0.900444097676710 -0.122511280687204;-0.365801039878279 -0.931107838994183 -0.236883085813983];
X = lyap(ALn',-delta,-CLn'*Grand)
L = (Grand*inv(X))'
eig(ALn-L*CLn) 
%% Observador de Kalman
sys_new = ss (ALn,[BLn BLn],CLn,[DLn DLn])
Qn = diag([0.04^2 0.01^2 0.025^2])  % ruido del actuador
Rn = diag([0.04^2 0.01^2]) % ruido del sensor
[sys_obv,L,S] = kalman(sys_new,Qn,Rn)
eig(sys_obv.a)
eig(ALn-L*CLn) 
%% Control 2 DOF Continua: Controlador & Observador
Ac_2 = [ALn-BLn*Kest-L*CLn+L*DLn*Kest, BLn*Ki - L*DLn*Ki;
        zeros(r,n), zeros(r,r)];
Bc_2 = [zeros(n,r),L;
        eye(r), -eye(r)];  
Cc_2 = [-Kest, Ki];    
Dc_2 = zeros(p,2*r);
sys_c2 = ss(Ac_2,Bc_2,Cc_2,Dc_2)
Bsp = Bc_2(:,1:r);
Bpv = Bc_2(:,r+1:2*r);
ALC = [ALn, BLn*Cc_2;
       Bpv*CLn, Ac_2+Bpv*DLn*Cc_2]
BLC = [zeros(n,r);
       Bsp]
CLC = [CLn, DLn*Cc_2]
DLC = zeros(r,r)
sysLC = ss(ALC,BLC,CLC,DLC)
CuLC = [zeros(p,n), Cc_2]
DuLC = zeros(p,r)
sysuLC = ss(ALC,BLC,CuLC,DuLC)
figure(3)
step(sysLC,5)
grid
figure(4)
step(sysuLC)
grid
%% Control 2 DOF discreta: Controlador & Observador
Tm = 0.01
sys_d2 = c2d(sys_c2,Tm)
sysd = c2d(sys,Tm)

Gp = sysd.a
Hp = sysd.b
Cp = sysd.c
Dp = sysd.d
Gc = sys_d2.a
Hc = sys_d2.b
Hsp = sys_d2.b(:,1:r)
Hpv = sys_d2.b(:,r+1:end)
Cc = sys_d2.c
Dc = sys_d2.d
GLc = [Gp Hp*Cc;
       Hpv*Cp, Gc+Hpv*Dp*Cc]
HLc = [zeros(n,r);
       Hsp]
CLc = [Cp Dp*Cc]
DLc = zeros(r,r)
CULc = [zeros(p,n), Cc]
DULc = zeros(p,r) 
sysdLc = ss(GLc, HLc, CLc, DLc, Tm)
sysULc = ss(GLc, HLc, CULc, DULc, Tm)
figure(5)
step(sysdLc,5)
grid
figure(6)
step(sysULc)
grid
%% Ganancia antiwindup
Kaw = [0 0 0;
       0 0 0;
       0 0 0;
       0 0 0;
       0 0 0]
%% comando LQG
Q = diag([1 1 1])  % Penalidad de los estados
R = diag([0.1 0.1 0.1]) % Penalidad de la se�al de control
QI = diag([30 10]) % Penalidad de la componente integral
Qn = diag([0.01^2 0.01^2 0.01^2])   % ruido de los estados
Rn = diag([0.01^2 0.01^2]) % ruido del sensor
QXU = blkdiag(Q,R);
QWV = blkdiag(Qn,Rn);
Reg_1dof = lqg(sys,QXU,QWV,QI,'1dof')
Reg_2dof = lqg(sys,QXU,QWV,QI,'2dof')