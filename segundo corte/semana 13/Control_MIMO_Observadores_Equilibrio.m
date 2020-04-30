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
p = 3 % Número de actuadores
r = 2 % Número de salidas
M = ctrb(ALn,BLn)
rank(M)
%% Control LQR
Q = diag([1 1 10]);
R = diag([3 3 3]);
[Kest,S,E] = lqr(ALn,BLn,Q,R)
%% Modelo lineal en lazo cerrado para pto equilibrio
xini = [0.5; -0.3; -0.7]
Anew = ALn-BLn*Kest
eig(Anew)
Bnew = zeros(n,p)
Cnew = CLn
Dnew = DLn
sys_new = ss(Anew,Bnew,Cnew,Dnew)
figure(1)
initial(sys_new,xini)
grid
Cu = -Kest
Du = zeros(p,p)
sys_u = ss(Anew,Bnew,Cu,Du)
figure(2)
initial(sys_u,xini)
grid
%% Obsevabor de Luenberger polos reales
Vob = obsv(ALn,CLn)
rank(Vob) % rango debe ser igual al numero de estados
ts = 0.2
s1 = -4/ts
delta = diag([s1 2*s1 4*s1]) % Matriz cuadrada igual al numero de estados n*n
% Grand = 2*rand(r,n)-1;  % Matriz de r*n
Grand = [-0.544671404366893 -0.377795426699174 -0.139585217340832;-0.128602631792202 0.846759284206488 -0.630367359751728]
X = lyap(ALn',-delta,-CLn'*Grand)
det(X)
L = (Grand*inv(X))'
eig(ALn-L*CLn) 
%% Obsevabor de Luenberger polos complejos
ts = 0.2
Mp = 0.1
zita = abs(log(Mp))/sqrt(pi^2+(log(Mp))^2)
wn = 4/(ts*zita)
wd = wn*sqrt(1-zita^2)
s1 = -zita*wn + wd*j
delta = [real(s1) imag(s1) 0;
         -imag(s1) real(s1) 0;
          0 0 5*real(s1)]
Grand = 2*rand(r,n)-1;  % Matriz de r*n
Grand = [0.602029245539478 0.857708278956089 -0.0227820523928417;-0.941559444875707 0.460661725710906 0.157050122046878] 
X = lyap(ALn',-delta,-CLn'*Grand)
L = (Grand*inv(X))'
eig(ALn-L*CLn) 
%% Observador de Kalman
sys_new = ss (ALn,[BLn BLn],CLn,[DLn DLn])
Qn = diag([0.07^2 0.0004^2 0.005^2])  % ruido del actuador
Rn = diag([0.001^2 0.001^2]) % ruido del sensor
[sys_obv,L,S] = kalman(sys_new,Qn,Rn)
eig(sys_obv.a)
eig(ALn-L*CLn) 
%% Control unificado continua: controlador & Observador
Ac = ALn - L*CLn - (BLn-L*DLn)*Kest
Bc = L
Cc = -Kest
Dc = zeros(p,r)
sys_c = ss(Ac,Bc,Cc,Dc)
ALc = [ALn, BLn*Cc;
       Bc*CLn, Ac+Bc*DLn*Cc]
BLc = zeros(2*n,1)
CLc = [CLn DLn*Cc]
DLc = zeros(r,1)

CuLc = [zeros(p,n) Cc]
DuLc = zeros(p,1)
xininew = [xini; zeros(n,1)]
sys_new_Lc = ss(ALc,BLc,CLc,DLc)
figure(3)
initial(sys_new_Lc,xininew)
grid

sys_u_Lc = ss(ALc,BLc,CuLc,DuLc)
figure(4)
initial(sys_u_Lc,xininew)
grid
%% Control unificado discreta: controlador & Observador
Tm = 0.01;
sys_c_dig = c2d(sys_c,Tm)
Gc = sys_c_dig.a
Hc = sys_c_dig.b
Cc = sys_c_dig.c
sys_d = c2d(sys,Tm)
Gp = sys_d.a
Hp = sys_d.b
Cp = sys_d.c
Dp = sys_d.d
GLc = [Gp, Hp*Cc;
       Hc*Cp, Gc+Hc*Dp*Cc]
HLc = zeros(2*n,1)
CLc = [Cp Dp*Cc]
DLc = zeros(r,1)

CuLc = [zeros(p,n) Cc]
DuLc = zeros(p,1)
xininew = [xini; zeros(n,1)]
sysd_Lc = ss(GLc,HLc,CLc,DLc,Tm)
figure(5)
initial(sysd_Lc,xininew)
grid

sysdu_Lc = ss(GLc,HLc,CuLc,DuLc,Tm)
figure(6)
initial(sysdu_Lc,xininew)
grid
%% comando LQG
Q = diag([1 1 1])  % Penalidad de los estados
R = diag([0.1 0.1 0.1]) % Penalidad de la señal de control
Qn = diag([0.01^2 0.01^2 0.01^2])   % ruido de los estados
Rn = diag([0.01^2 0.01^2]) % ruido del sensor
QXU = blkdiag(Q,R);
QWV = blkdiag(Qn,Rn);
Reg = lqg(sys,QXU,QWV)
%%
Ac = Reg.a
Bc = Reg.b
Cc = Reg.c
Dc = Reg.d
ALc = [ALn, BLn*Cc;
       Bc*CLn, Ac+Bc*DLn*Cc]
BLc = zeros(2*n,1)
CLc = [CLn DLn*Cc]
DLc = zeros(r,1)

CuLc = [zeros(p,n) Cc]
DuLc = zeros(p,1)
xininew = [xini; zeros(n,1)]
sys_new_Lc = ss(ALc,BLc,CLc,DLc)
figure(3)
initial(sys_new_Lc,xininew)
grid

sys_u_Lc = ss(ALc,BLc,CuLc,DuLc)
figure(4)
initial(sys_u_Lc,xininew)
grid