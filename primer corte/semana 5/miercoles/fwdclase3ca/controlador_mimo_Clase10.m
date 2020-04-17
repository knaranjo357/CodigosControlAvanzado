clc
clear
close all
A = [-0.005509 0 0 -0.1588;
    0 -0.2062 0 0;
    -0.01216 0 0 -0.5672;
    0 0 0 -0.040]
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
figure(1)
step(sys)
grid
Gf = tf(sys)
%%
s = tf('s')
Kst = dcgain(Gf)
Kv = dcgain(s*Gf)
K = [Kst(1:2,:);
     Kv(3,:);
     Kst(4,:)]
RGA = K.*(inv(K))'
%% Diseño PID 1 DOF
G14 = Gf(1,4)
C14 = -0.02*(s+0.038)*(s+0.01)/(s*(s+0.1)) % Mp = 8.1%, tr = 75.4 sg, ts = 355 sg
G22 = Gf(2,2)
C22 = 0.0005*(s+0.6)/s % Mp = 0%, tr = 190 sg, ts = 341 sg
G33 = Gf(3,3)
C33 = 0.02/(s+0.2) % Mp = 11.5%, tr = 76.4 sg, ts = 327 sg
G41 = Gf(4,1)
C41 = 0.0015*(s+0.67)/s % Mp = 0%, tr = 187 sg, ts = 312 sg
%% Análisis de acoplamiento
r = 4 % Número de salidas
Gc = [0 0 0 C41;
    0 C22 0 0;
    0 0 C33 0;
    C14 0 0 0]
T = feedback(sys*Gc,eye(r))
U = feedback(Gc,sys)
figure(2)
step(T)
grid
figure(3)
step(U)
grid
%% Análisis de condiciones iniciales
xini = 0*[0.2; 0.3; 0; -0.2]
xo = [xini; zeros(5,1)]
figure(2)
initial(T,xo)
grid
figure(3)
initial(U,xo)
grid

%% Diseño PID 2 DOF
G14 = Gf(1,4)
C14a = -0.7*(s+0.014)/s
C14b = -55*(s+0.005)/(s+0.4516)

G22 = Gf(2,2)
C22a = 1.5*(s+0.014)/s
C22b = 55*(s+0.005)/(s+0.452)

G33 = Gf(3,3)
C33a = 0.072/s
C33b = 3.8*(s+0.002)/(s+0.012)

G41 = Gf(4,1)
C41a = 0.072/s
C41b = 3.8*(s+0.002)/(s+0.012)

T1 = feedback(C41a*G41,1+C41b/C41a)
U1 = feedback(C41a,G41*(1+C41b/C41a))
figure(1)
step(T1)
grid
figure(2)
step(U1)
grid