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
%% Diseño Controlador polos reales
ts = 3;
s1 = -4/ts
delta = diag([s1, 1.2*s1, 3*s1, 6*s1, 10*s1]) % Matriz cuadrada de n+r
% Grand = 2*rand(p,n+r)-1;
Grand = [0.402197511801853,0.396211040360617,-0.743971200559655,-0.934798358938944,0.338350609068788;0.332677703168851,0.333055826805174,0.998160789522721,0.122399585419320,-0.619133465640092;0.0782529300857133,-0.643735091199325,-0.657757867287136,0.763733000903620,-0.262166907872210];
det(X)
K = Grand*inv(X)
Kes = K(:,1:n)
Ki = -K(:,n+1:end)
%% Diseño polos complejos
ts = 3;
Mp = 0.05;
zita = abs(log(Mp))/sqrt(pi^2+(log(Mp))^2);
wn = 4/(zita*ts)
wd = wn*sqrt(1-zita^2)
s1 = -zita*wn + j*wd
delta = [real(s1) imag(s1) 0 0 0;  % Matriz cuadrada de n+r
        -imag(s1) real(s1) 0 0 0;
        0 0 2.5*real(s1) 0 0;
        0 0 0 3.1*real(s1) 0;
        0 0 0 0 8*real(s1)]
% Grand = 2*rand(p,n+r)-1;
Grand = [0.308891415514133,0.436717886411767,-0.349708636358880,0.557604483648185,-0.467057018441855;-0.184761605917695,0.937298660462187,-0.788741593341956,-0.153094162074523,-0.692686564817387;0.639962445563881,0.0626678131313490,0.221917317492401,-0.818353428425121,-0.437989394932258];
X = lyap(Ahat, -delta, -Bhat*Grand);
det(X)
K = Grand*inv(X)
Kes = K(:,1:n)
Ki = -K(:,n+1:end)
%% Control óptimo
Q = diag([1 1 1 1 1])      % Matriz diagonal de tamaño n+r
R = diag([0.1 0.1 0.1])           % Matriz diagonal de tamaño p
[K,S,E] = lqr(Ahat,Bhat,Q,R)
[K1,S,E] = lqi(sys,Q,R)
Kes = K(:,1:n)
Ki = -K(:,n+1:end)
%% Modelo en lazo cerrado
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

% Acción de control
CC1 = [-Kes, Ki];
DD1 = zeros(p,r);
sys_u = ss(AA,BB,CC1,DD1)
figure(2)
step(sys_u)
grid
