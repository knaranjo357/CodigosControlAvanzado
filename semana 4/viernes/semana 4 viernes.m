clc
clear
A = [0 1;
    -6 -5]
B = [1 0;
    0 1]
C = [1 2]
D = [0 1]
sys = ss(A,B,C,D)
Tm = 0.1
sysd = c2d(sys,Tm)
GL = sysd.a
HL = sysd.b
den = poly(eig(GL))
[r,p,k] = residue([1 -0.9746],den)
xini = [0.3; 0.1]
%% 0<t<1
syms z n 
inv_mat = inv(z*eye(2)-GL)
phik = iztrans(inv_mat*z)
xk = phik*xini;
xknum = double(subs(xk,n,10))
yk = C*xk;
yknum = double(subs(yk,n,9))
%%
k = 10;
xk1 = [(0.8187)^k - 0.7*(0.7408)^k;
      -2*(0.8187)^k + 2.1*(0.7408)^k]
yk = -3*(0.8187)^k + 3.5*(0.7408)^k
%% 1<t
uz = [z/(z-1);
      -2*z/(z-1)]
  
yz2 = C*inv_mat*(z*xknum + HL*uz) + D*uz;
yz2 = collect(yz2,z)
[a1,b1] = numden(yz2)
n1 = sym2poly(a1)
d1 = sym2poly(b1)
[r,p,k] = residue(n1(1,end-1),d1)
yk2 = iztrans(yz2)
yk2num = double(subs(yk2,n,10))





%%%incertidumbres
clc
clear
Tol = 5
a = ureal('a',-6,'Percentage',[-Tol,Tol])
b = ureal('b',-5,'Percentage',[-Tol,Tol])
get(a)
c = a+b
c = a*b
c = a/b
c = 2*a
c = a^2
c = ureal('c',exp(-5),'Percentage',[-Tol,Tol])
A = [0 1;
    a -5]
B = [1 0;
    0 1]
C = [1 2]
D = [0 1]
sys = ss(A,B,C,D)
sysnom = sys.NominalValue
figure(1)
step(sys)
grid
figure(2)
bodemag(sys)
grid
