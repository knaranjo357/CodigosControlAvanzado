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
%% Matriz de controlabilidad
n = 3 % Numero de estados
p = 3 % Número de actuadores
r = 2 % Número de salidas
Ahat = [ALn zeros(n,r);
        -CLn zeros(r,r)];
Bhat = [BLn;
        -DLn]    
M = ctrb(Ahat,Bhat)
rank(M)
%% Diseño del Controlador
Q = diag([1 1 1 50 20])
R = diag([0.1 0.1 0.1])
[K,S,E] = lqi(sys,Q,R)
Kest = K(:,1:n)
Ki = -K(:,n+1:end)
%%
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

% Acción de control
CC1 = [-Kest Ki];
DD1 = zeros(p,r);
sys_u = ss(AA,BB,CC1,DD1)
figure(2)
step(sys_u)
grid

%% Obsevabor de Luenberger polos reales
Vob = obsv(ALn,CLn)
rank(Vob) % rango del numero de estados
ts = 0.3 %por lo menos 5 veces mas rapido que el sistema
s1 = -4/ts
delta = diag([s1 2*s1 4*s1]) % Matriz cuadrada igual al numero de estados n*n
%%Grand = 2*rand(r,n)-1;  % Matriz de r*n
Grand = [0.914333896485891 0.600560937777600 -0.156477434747450;-0.0292487025543176 -0.716227322745569 0.831471050378134];
X = lyap(ALn',-delta,-CLn'*Grand)
det(X)
L = (Grand*inv(X))'
eig(ALn-L*CLn) 
%% Unificaciòn del controlador y observador
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
