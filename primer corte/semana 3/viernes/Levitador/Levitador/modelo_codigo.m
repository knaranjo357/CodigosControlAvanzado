clc
clear
g = 9.81;
M = 0.05;
R = 10;
L = 0.5;
x10 = 0.16
x20 = 0;
x30 = sqrt(g*M*x10)
eo = R*x30
%%
syms x1 x2 x3 e g1 M1 R1 L1
f1 = x2;
f2 = g1 - x3^2/(M1*x1);
f3 = -R1*x3/L1 + e/L1;
F = [f1;f2;f3];
X = [x1;x2;x3];
U = [e];
% Y = [x1];
Y = [e*x3;
     R1*x3^2]
A = jacobian(F,X)
B = jacobian(F,U)
C = jacobian(Y,X)
D = jacobian(Y,U)
AL = double(subs(A,{x1,x2,x3,e,g1,M1,R1,L1},{x10,x20,x30,eo,g,M,R,L}))
BL = double(subs(B,{x1,x2,x3,e,g1,M1,R1,L1},{x10,x20,x30,eo,g,M,R,L}))
CL = double(subs(C,{x1,x2,x3,e,g1,M1,R1,L1},{x10,x20,x30,eo,g,M,R,L}))
DL = double(subs(D,{x1,x2,x3,e,g1,M1,R1,L1},{x10,x20,x30,eo,g,M,R,L}))
sys = ss(AL,BL,CL,DL)
eig(AL)
G = tf(sys)