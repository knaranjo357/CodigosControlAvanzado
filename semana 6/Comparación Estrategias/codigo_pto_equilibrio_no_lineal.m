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
Mc = ctrb(AL,BL)
rank(Mc) % El rango debe ser igual a n
xini = [0.2; 0.3; 0.01]
%% Diseño Controlador polos reales
ts = 3;
s1 = -4/ts
delta = diag([s1, 3*s1, 3*s1])  % Matriz cuadrada de n
% Grand = 2*rand(p,n)-1;
Grand = [-0.622089969934911,-0.263030807019327,-0.837748462268430;0.373550866730630,0.251237121459381,0.858771941937460;-0.632977688525461,0.560454870302754,0.551425357216805];
X = lyap(AL, -delta, -BL*Grand);
det(X)
Kest = Grand*inv(X)
eig(AL-BL*Kest)
%% Diseño polos complejos
ts = 3;
Mp = 0.05;
zita = abs(log(Mp))/sqrt(pi^2+(log(Mp))^2);
wn = 4/(zita*ts)
wd = wn*sqrt(1-zita^2)
s1 = -zita*wn + j*wd
delta = [real(s1) imag(s1) 0;
        -imag(s1) real(s1) 0;
        0 0 1.2*real(s1)]
% Grand = 2*rand(p,n)-1;
Grand = [4.48711804017510e-05,0.219733296845117,0.610978849059371;-0.0401557177078791,0.235332779176909,0.153443031229370;0.809444476134726,0.718884611292425,-0.634155061170172];
X = lyap(AL, -delta, -BL*Grand);
det(X)
Kest = Grand*inv(X)
%% Control óptimo
Q = diag([1 1 6])  %Matriz diagonal de tamaño n
R = diag([1 1 0.1])      %Matriz diagonal de tamaño p
[Kest,S,E] = lqr(AL,BL,Q,R)
%% Modelo en lazo cerrado
% Anew = Ac-BL*Kest
% eig(Anew.NominalValue)
Anew = AL-BL*Kest
eig(Anew)
Bnew = zeros(n,p)
Cnew = CL
Dnew = DL
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

