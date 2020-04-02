clc
clear all
Ra = 50;
Ki = 1.5;
Kv = 1.5;
Bm = 0.01;
JL = 1;
n = 10;
Jm = 0.1;

%Con w1 como estado
A = [-((Ki*Kv)/(Jm*Ra+(Ra*JL/n^2)))-(Bm/(Jm+(JL/n^2)))]

B = [Ki/(Jm*Ra+(Ra*JL/n^2))   -1/(Jm*n+JL/n)]

C = [1/n;
    -Kv/Ra]

D = [0 0;
    1/Ra 0]

sys = ss(A,B,C,D)
Gf = tf(sys)
Gf = zpk(sys)
G11 = Gf(1,1)
G12 = Gf(1,2)
G21 = Gf(2,1)
G22 = Gf(2,2)

step(Gf)
grid on
step(G11)

%% Diseño PID 1 DOF
s = tf('s')
Con = 5.6*(s+1.5)*(s+0.7)/(s*(s+1))

%% Diseño 2grados de libertad 
Ca = (17*(s+0.5))/s
Cb = 2*s/(s+0.5)
%% Cascada
Ga = Gf(2,1)
Gb = 0.02727/(0.02*s+0.001818)
C1 = 7.1/(s*(s+0.11))
C2 = 0.08*(s+0.094)/s
T = feedback(Ga*C1,1)*Gb
%Si una funcion de transferencia tiene igual numero de polos y ceros el
%valor inicial va a ser la division de las partes que multiplican a las "s"
%de mayor orden 
%Con w2 como estado
%{
Jt = 1 + n^2*Jm/JL

A = [(-((Ki*Kv*n^2)/(JLRa))-(Bm*n^2/(JL)))/Jt]

B = [n*Ki/(JL*R*Jt),   -1/(JL*Jt)]

C = [1;
    -Kv*n/Ra]

D = [0 0;
    1/Ra 0]

figure(1)
step(Gf)
grid

figure(2)
bodemag(Gf)
grid

roots([1 2.3 4.1 3.6 2.4])
%%
xini = [ 0.2; 0; 0.3; 0]
%}

