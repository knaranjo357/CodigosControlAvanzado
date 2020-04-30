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
u10 = double((a1*x20 + b1*x10^2)/u20)
u30 = double((d1*x30 - u10)/x10)
ALn = double(subs(AL,{a,b,c,d},{a1,b1,c1,d1}))
BLn = double(subs(BL,{a,b,c,d},{a1,b1,c1,d1}))
CLn = [1 0 0;
       0 0 1]
DLn = [0 0 0;
       0 0 0]
sys =ss(ALn,BLn,CLn,DLn)
eig(ALn)
Gs = tf(sys)
K = dcgain(Gs)
RGA = K.*pinv(K)'
%% Seleccion periodo de muestreo
bodemag(sys,{0.01 10000})
grid
%%
Tm = 10e-3
%%
sysd = c2d(sys,Tm)
figure(1)
bodemag(sys,sysd)
grid
figure(2)
step(sys,sysd)
grid
GLn = sysd.a
HLn = sysd.b
%% Análisis RGA en discreta
Gz = tf(sysd)
Kz = dcgain(Gz)
RGAz = Kz.*pinv(Kz)'
G11 = zpk(Gz(1,1))
G23 = zpk(Gz(2,3))
%% PID 1 DOF
z = tf('z',Tm);
C11 = 2*(z-0.84)/(z-1)
C23 = 7*(z-0.8)/(z-1)

% Análisis acoplamiento
Gcz = [C11 0;
    0 0;
    0 C23];
r = 2 % Variables a controlar
Tz = feedback(sysd*Gcz,eye(r))
Uz = feedback(Gcz,sysd)

figure(1)
step(Tz)
grid
figure(2)
step(Uz)
grid
%% PID 2 DOF

C11a= 6.8*(z-0.85)/(z-1)
C11b= 0.07*(z-0.77)/(z-0.93)

C23a= 2.5*(z-0.6)/(z-1)
C23b= 4*(z-0.4)/(z-0.8)


%% Seguimiento por espacio de estados
% Matriz de controlabilidad
n = 3 % Número de estados
p = 3 % Número de actuadores
r = 2 % Variables a controlar
Ghat = [GLn, zeros(n,r);
        -CLn*GLn, eye(r)]
Hhat = [HLn;
        -CLn*HLn]
Mc = ctrb(Ghat,Hhat)
rank(Mc)  % El rango de la matirz es n+r
%% Diseño polos reales
ts = 1;
%muestras = ts/Tm
s1 = -4/ts
Polc = [s1, 1.2*s1, 5*s1, 5*s1, 6*s1]
    Pold = exp(Tm*Polc)
    delta = diag(Pold)  % Matriz cuadrada de n+r
% Grand = 2*rand(p,n+r)-1;
Grand = [0.412092176039218 -0.907657218737692 0.389657245951634 -0.931107838994183 0.531033576298005;-0.936334307245159 -0.805736437528305 -0.365801039878279 -0.122511280687204 0.590399802274126;-0.446154030078220 0.646915656654585 0.900444097676710 -0.236883085813983 -0.626254790891243];
    X = lyap(Ghat, -delta, -Hhat*Grand);
    det(X)
    K = Grand*inv(X)
    Kes = K(:,1:n)
    Ki = -K(:,n+1:end)
%% Diseño polos complejos
ts = 1;
% muestras = ts/Tm
Mp = 0.1
    zita = abs(log(Mp))/sqrt(pi^2+(log(Mp))^2);
    wn = 4/(zita*ts)
    wd = wn*sqrt(1-zita^2)
    s1 = -zita*wn + j*wd
Polc = [s1, conj(s1), 2*real(s1), 5*real(s1), 6*real(s1)]
    Pold = exp(Tm*Polc)
delta = [   real(Pold(1))   imag(Pold(1))       0               0           0; 
            imag(Pold(2))   real(Pold(2))       0               0           0;
            0               0                   Pold(3)         0           0;
            0               0                   0               Pold(4)     0;
            0               0                   0               0           Pold(5)]
%Grand = 2*rand(p,n+r)-1;
Grand = [-0.300032468030383 0.232089352293278 0.661657255792582 0.834387327659620 0.507458188556991;-0.606809499137584 -0.0534223021945415 0.170528182305449 -0.428321962359253 -0.239108306049287;-0.497832284047938 -0.296680985874007 0.0994472165822791 0.514400458221443 0.135643281450442];
    X = lyap(Ghat, -delta, -Hhat*Grand);
    det(X)
    K = Grand*inv(X)
    Kes = K(:,1:n)
    Ki = -K(:,n+1:end)
%% Control óptimo
Q = diag([1 1 1 0.1 0.1]) % Matriz diagonal de tamaño n+r si Q es mas grande es mas rapido
%%primero se cambia lo ultimo de Q luego se cambia R y los valores iniciales de Q
R = diag([10 10 10])           % Matriz diagonal de tamaño p
%si R es mas pequeño es mas rapido
[K,S,E] = dlqr(Ghat,Hhat,Q,R)
Kes = K(:,1:n)
Ki = -K(:,n+1:end)
%% Modelo lineal en lazo cerrado para seguimiento
    AA = [GLn - HLn*Kes, HLn*Ki;
         -CLn*GLn + CLn*HLn*Kes ,eye(r)-CLn*HLn*Ki]
    eig(AA)
    BB = [zeros(n,r);
          eye(r)]
    CC = [CLn, zeros(r,r)]
    DD = zeros(r,r);
    sys_new = ss(AA,BB,CC,DD,Tm)
    figure(1)
    step(sys_new)
    grid

    % Acción de control
    CC1 = [-Kes Ki]
    CC1 = -K
    DD1 = zeros(p,r);
    sys_u = ss(AA,BB,CC1,DD1,Tm)
    figure(2)
    step(sys_u)
    grid
%% Diseño punto de equilibrio
Mc = ctrb(GLn,HLn)
rank(Mc)  % El rango de la matirz es n
% Diseño polos reales
ts = 1.5;
muestras = ts/Tm
s1 = -4/ts
Polc = [s1, 2*s1, 4*s1]
Pold = exp(Tm*Polc)
delta = diag(Pold)  % Matriz cuadrada de n
%Grand = 2*rand(p,n)-1;
Grand=[0.929777070398553,0.9914333896485891,-0.716227322745569;-0.684773836644903,-0.0292487025543176,-0.156477434747450;0.941185563521231,0.600560937777600,0.831471050378134]
X = lyap(GLn, -delta, -HLn*Grand);
det(X)
K = Grand*inv(X)
%% Diseño polos complejos
ts = 1.5;
muestras = ts/Tm
Mp = 0.1
zita = abs(log(Mp))/sqrt(pi^2+(log(Mp))^2);
wn = 4/(zita*ts)
wd = wn*sqrt(1-zita^2)
s1 = -zita*wn + j*wd
Polc = [s1, conj(s1), 4*real(s1)]
Pold = exp(Tm*Polc)
delta = [real(Pold(1)) imag(Pold(1)) 0;  % Matriz cuadrada de n
        imag(Pold(2)) real(Pold(2)) 0;
        0 0 Pold(3)]
% Grand = 2*rand(p,n)-1;
Grand = [0.584414659119109 -0.928576642851621 0.357470309715547;0.918984852785806 0.698258611737554 0.515480261156667;0.311481398313174 0.867986495515101 0.486264936249832];
X = lyap(GLn, -delta, -HLn*Grand);
det(X)
K = Grand*inv(X)
%% Control óptimo
Q = diag([2 0.1 0.8]) % Matriz diagonal de tamaño n+r
R = diag([10 10 1])           % Matriz diagonal de tamaño p
[K,S,E] = dlqr(GLn,HLn,Q,R)
%% Modelo lineal en lazo cerrado para punto equilibrio
xini = [0.2; 0.6; 0.1]
AA = GLn-HLn*K
eig(AA)
BB = zeros(n,p)
CC = CLn-DLn*K
DD = zeros(r,p);
sys_new = ss(AA,BB,CC,DD,Tm)
figure(1)
initial(sys_new,xini)
grid

% Acción de control
CC1 = -K;
DD1 = zeros(p,p);
sys_u = ss(AA,BB,CC1,DD1,Tm)
figure(2)
initial(sys_u,xini)
grid
%% Pto equilibrio del PID 1 DOF
xo = [xini; zeros(2,1)]
figure(1)
initial(Tz,xo)
grid
figure(2)
initial(Uz,xo)
grid