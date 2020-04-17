%% PARCIAL
clc
clear
close all
%paramtros
syms a b c d e f g 
an=3
bn=4
cn=2
dn=9
en=4
fn=5
gn=2
%se asume
u10=2;
u20=4;
u30=1;
%punto de equilibrio en parametros 
x20=0
x30=fn*u30/(gn*u20)
x10=-bn/cn*x30*x20+dn/cn*u10*u20

%% sistema lineal en parametros
syms x1 x2 x3 u1 u2 u3 
f1 = x2;
f2 = -b/a*x3*x2-c/a*x1+d/a*u1*u2;
f3 = f/e*u3-g/e*x3*u2;
F=[f1;f2;f3];
X=[x1;x2;x3];
U=[u1;u2;u3];
y1=x1;
y2=x3;
Y=[y1;y2];
A=jacobian(F,X)
B=jacobian(F,U)
C = jacobian(Y,X)
D = jacobian(Y,U)

AL = double(subs(A,{x1,x2,x3,u1,u2,u3,a,b,c,d,e,f,g},{x10,x20,x30,u10,u20,u30,an,bn,cn,dn,en,fn,gn}))
BL = double(subs(B,{x1,x2,x3,u1,u2,u3,a,b,c,d,e,f,g},{x10,x20,x30,u10,u20,u30,an,bn,cn,dn,en,fn,gn}))
CL = double(subs(C,{x1,x2,x3,u1,u2,u3,a,b,c,d,e,f,g},{x10,x20,x30,u10,u20,u30,an,bn,cn,dn,en,fn,gn}))
DL = double(subs(D,{x1,x2,x3,u1,u2,u3,a,b,c,d,e,f,g},{x10,x20,x30,u10,u20,u30,an,bn,cn,dn,en,fn,gn}))
eig(AL)

%% puntos de equilibrio salida
x0=[x10;x20;x30]
y0=CL*x0;
y10=y0(1,1)
y20=y0(2,1)
%% sys
sys = ss(AL,BL,CL,DL);
Gf = tf(sys);
eig(AL)
figure(1)
step(sys)
grid

%% modelo discreta 
%periodo muestre0
w3db=2.41;
ws=20*w3db;
Tm=2*pi/ws

sysd = c2d(sys,Tm)
GL = sysd.a
HL = sysd.b

figure(1)
step(sys,sysd)
grid
figure(2)
bodemag(sys,sysd)
grid
%% Respuesta transitoria
syms s
mat_inv=inv(s*eye(3)-AL)

x_ini=[0.01;0.1;-0.02]
us=[0;0;0] %entrada

%salida
ys=C*mat_inv*(x_ini+BL*us)+DL*us
syms t
yt=collect(ilaplace(ys, t))

kn = [0,1,2,3,4]
tn=kn*Tm
puntos_dc=double(subs(yt,{t},{tn}))

ut = us;
yo = CL*x_ini + DL*ut;
x1d = GL*x_ini +  HL*ut;
y1 = CL*x1d  + DL*ut;
x2d = GL*x1d + HL*ut;;
y2 = CL*x2d  + DL*ut;
x3d = GL*x2d + HL*ut;;
y3 = CL*x3d  + DL*ut;
x4d = GL*x3d + HL*ut;;
y4 = CL*x4d + DL*ut;;

punto_dis = double([yo y1 y2 y3 y4])

%% Matriz de ganancias relatias RGA
s=tf('s');
Gp = tf(sys)
K = dcgain(Gp)
RGA = K.*(pinv(K))'

%% ---------------------------- CONTROL PID 1DOF --------------------------------------------

%%
G11=Gp(1,1)
C11=27.421*(s+0.856)*(s+0.9955)/(s*(s+111.3))

G23=Gp(2,3)
C23=1.8527*(s+9.8)*(s+3.56)/(s*(s+7.56))

%% Análisis de acoplamiento controlador 1DOF
r = 2 % Número de salidas
Gc = [C11 0;
    0 0;
    0 C23];
T = feedback(sys*Gc,eye(r))
U = feedback(Gc,sys)
figure(2)
step(T)
grid
figure(3)
step(U)
grid

%% Análisis de condiciones iniciales para modelo lineal 
xini = 1.2*[x10; 0.1; x30]

xo = [xini; zeros(4,1)]
figure(4)
initial(T,xo)
grid
figure(5)
initial(U,xo)
grid

%% ---------------------------- CONTROL MULTIVARIABLE SEGUIMIENTO ----------------------------------------

%% Matriz controlabilidad
n = 3; % Número de estados
p = 3; % Número de actuadores
r = 2; % Número de salidas
Ahat = [AL, zeros(n,r);
        -CL, zeros(r,r)]
Bhat = [BL;
        -DL]
Mc = ctrb(Ahat,Bhat)
rank(Mc) % El rango debe ser igual a n+r

%% Diseño polos complejos
ts = 3;
Mp = 0.05
zita = abs(log(Mp))/sqrt(pi^2+(log(Mp))^2);
wn = 4/(zita*ts)
wd = wn*sqrt(1-zita^2)
s1 = -zita*wn + j*wd
delta = [real(s1) imag(s1) 0 0 0;  % Matriz cuadrada de n+r
        -imag(s1) real(s1) 0 0 0;
        0 0 1.2*real(s1) 0 0 ;
        0 0 0 20*real(s1) 0 ;
        0 0 0 0 25*real(s1)]
    
% Grand = 2*rand(p,n+r)-1;
Grand = [ -0.8483    0.5583    0.1376   -0.3258   -0.3776;
   -0.8921    0.8680   -0.0612   -0.6756    0.0571;
    0.0616   -0.7402   -0.9762    0.5886   -0.6687];
    
X = lyap(Ahat, -delta, -Bhat*Grand);
det(X)
K = Grand*inv(X)
Kes_PC = K(:,1:n)
Ki_PC = -K(:,n+1:end)

%%Modelo en lazo cerrado
Ki=Ki_PC;
Kes=Kes_PC;

AA = [AL - BL*Kes, BL*Ki;
     -CL + DL*Kes , -DL*Ki]
eig(AA)
BB = [zeros(n,r);
      eye(r)]
CC = [CL - DL*Kes, DL*Ki]
DD = zeros(r,r);
sys_new = ss(AA,BB,CC,DD)
figure(1)
step(sys_new)
grid

%%Acción de control
CC1 = [-Kes, Ki];
DD1 = zeros(p,r);
sys_u = ss(AA,BB,CC1,DD1)
figure(2)
step(sys_u)
grid

%% ---------------------------- CONTROL MULTIVARIABLE PTO DE EQUILIBRIO ----------------------------------------

%% Matriz controlabilidad
n = 3; % Número de estados
p = 3; % Número de actuadores
r = 2; % Número de salidas
Mc = ctrb(AL,BL)
rank(Mc) % El rango debe ser igual a n
xini = [x10; 0.1; x30]

%% Control óptimo
Q = diag([20 1 1])  %Matriz diagonal de tamaño n
R = diag([1 1 1])      %Matriz diagonal de tamaño p
[Kest_LQR,S,E] = lqr(AL,BL,Q,R)

%%Modelo en lazo cerrado
Kest=Kest_LQR;
Anew = AL-BL*Kest
eig(Anew)
Bnew = zeros(n,p)
Cnew = CL
Dnew = DL
sys_new = ss(Anew,Bnew,Cnew,Dnew)
figure(1)
initial(sys_new,xini)
grid

%accion de control
Cu = -Kest
Du = zeros(p,p)
sys_u = ss(Anew,Bnew,Cu,Du)
figure(2)
initial(sys_u,xini)
grid

%% Condicion inicial
xini1 = 0.05*xini
