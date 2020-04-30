clc
clear
A = [-0.005509 0 0 -0.1588;
      0 -0.2062 0 0;
      -0.01216 0 0 -0.5672;
      0 0 0 -0.04]
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
%% Matriz de controlabilidad
n = 4 % Numero de estados
p = 4 % Número de actuadores
r = 4 % Número de salidas
Ahat = [A zeros(n,r);
        -C zeros(r,r)];
Bhat = [B;
        -D]    
M = ctrb(Ahat,Bhat)
rank(M)
%% Matriz de observabilidad
Vob = obsv(A,C)
rank(Vob)
%%
Q = diag([1 1 1 1])  % Penalidad de los estados
R = diag([100 100 100 100]) % Penalidad de la señal de control
QI = diag([1 1 1 1]) % Penalidad de la componente integral
Qn = diag([0.01^2 0.01^2 0.01^2 0.01^2])   % ruido de los estados
Rn = diag([0.01^2 0.01^2 0.01^2 0.01^2]) % ruido del sensor
QXU = blkdiag(Q,R);
QWV = blkdiag(Qn,Rn);
Reg_1dof = lqg(sys,QXU,QWV,QI,'1dof')
Reg_2dof = lqg(sys,QXU,QWV,QI,'2dof')
%% Control 1 DOF Continua: Controlador & Observador
sysLC = feedback(sys*Reg_1dof,eye(r))
sysuLC = feedback(Reg_1dof,sys)
figure(1)
step(sysLC)
grid
figure(2)
step(sysuLC)
grid
%% Control 2 DOF Continua: Controlador & Observador
Ac_2 = Reg_2dof.a;
Bc_2 = Reg_2dof.b;  
Cc_2 = Reg_2dof.c;    
Dc_2 = Reg_2dof.d;
Bsp = Bc_2(:,1:r);
Bpv = Bc_2(:,r+1:2*r);
ALC = [A, B*Cc_2;
       Bpv*C, Ac_2+Bpv*D*Cc_2]
BLC = [zeros(n,r);
       Bsp]
CLC = [C, D*Cc_2]
DLC = zeros(r,r)
sysLC = ss(ALC,BLC,CLC,DLC)
CuLC = [zeros(p,n), Cc_2]
DuLC = zeros(p,r)
sysuLC = ss(ALC,BLC,CuLC,DuLC)
figure(3)
step(sysLC)
grid
figure(4)
step(sysuLC)
grid