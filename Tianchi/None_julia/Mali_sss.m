clear all;
close all;
clc
format long;


M1=1; 
M2=1;
M3=1;

% X=OCT
%Y=TET
%Z=NANOG
%1% Dx00+NT2 -> Dx01
%2% Dx10+NT2 -> Dx11
%3% Dx00+O -> Dx10
%4% Dx01+O -> Dx11
%5% 

G=[0,   0,  -1,   -1,    0,   0,  -1,  -1,   0,   0,   1,   1,   0,   0,  1,  1, -1,  1, 1,  1,  1,  1,  0,  0, 0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0;%X O%1
0,   0,   0,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0, 0,  0,  0,  0,  1, 1,  1,  1,  0,  0,  0, -1,  0, -1,  1,  0,  0;%Y  T%2
0,   0,   0,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  1,  1,  0,  0, -1, -1,  1,  0,  0; %Z3 N
-1,   0,  -1,    0,    0,   0,   0,   0,   1,   0,   1,   0,   0,   0,  0,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dx00  %4
1,   0,   0,   -1,    0,   0,   0,   0,  -1,   0,   0,   1,   0,   0,  0,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dx10 %5
0,  -1,   1,    0,    0,   0,   0,   0,   0,   1,  -1,   0,   0,   0,  0,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dx01  %6
0,   1,   0,    1,    0,   0,   0,   0,   0,  -1,   0,  -1,   0,   0,  0,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dx11  %7
0,   0,   0,    0,   -1,   0,  -1,   0,   0,   0,   0,   0,   1,   0,  1,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dy00  %8
0,   0,   0,    0,    1,   0,   0,  -1,   0,   0,   0,   0,  -1,   0,  0,  1,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dy10 %9
0,   0,   0,    0,    0,  -1,   1,   0,   0,   0,   0,   0,   0,   1, -1,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dy01  %10
0,   0,   0,    0,    0,   1,   0,   1,   0,   0,   0,   0,   0,  -1,  0, -1,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dy11  %11
0,   0,   0,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0, -1,  1, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dz0  %12
0,   0,   0,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  1, -1, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;%Dz1 %13
0,   0,   0,    0,    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -2,  2;%NT  %14
-1,  -1,   0,    0,   -1,  -1,   0,   0,   1,   1,   0,   0,   1,   1,  0,  0,  0,  0, 0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1];%NT2  %15
syms k1 k2 k3 k4 k5 k6 k7
k1=0.2; k2=0.3; k3=0.1; k4=1;
e=1; % make large for fast

a=[e  ;e   ;e   ;e   ;e   ;e   ;e   ;e   ;
e*k1;e*k1;e*k2;e*k2;e*k1;e*k1;e*k2;e*k2;
e;e*k2;
0;0;0;1;
0;0;0;1;
0;1;
1;1;1;
1*e;k3*e;
1*e;k4*e];


R=@(x) [
a(1)*x(4)*x(15);
a(2)*x(6)*x(15);
a(3)*x(4)*x(1);
a(4)*x(5)*x(1);
a(5)*x(8)*x(15);
a(6)*x(10)*x(15);
a(7)*x(8)*x(1);
a(8)*x(9)*x(1);%
a(9)*x(5);
a(10)*x(7);
a(11)*x(6);
a(12)*x(7);
a(13)*x(9);
a(14)*x(11);
a(15)*x(10);
a(16)*x(11);%
a(17)*x(12)*x(1);
a(18)*x(13);%
a(19)*x(4);
a(20)*x(5);
a(21)*x(6);
a(22)*x(7);%
a(23)*x(8);
a(24)*x(9);
a(25)*x(10);
a(26)*x(11);%
a(27)*x(12);
a(28)*x(13);%
a(29)*x(1);
a(30)*x(2);
a(31)*x(3);%
a(32)*x(2)*x(3);
a(33)*x(14);%
a(34)*x(14)^2;
a(35)*x(15)];


syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16
x=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12 x13 x14 x15].';





F=G*R(x)+[0.00;0.05;0;zeros(12,1)];
%% subs with my variables
syms Dn1 Dn0 N T O NT NT2 Do00 Do01 Do10 Do11 Dt00 Dt01 Dt10 Dt11 Dn0 Dn1 NT NT2
F = subs(F,'x1',O)
F = subs(F,'x2',T)
F = subs(F,'x3',N)
F = subs(F,'x4',Do00)
F = subs(F,'x5',Do10)
F = subs(F,'x6',Do01)
F = subs(F,'x7',Do11)
F = subs(F,'x8',Dt00)
F = subs(F,'x9',Dt10)
F = subs(F,'x10',Dt01)
F = subs(F,'x11',Dt11)
F = subs(F,'x12',Dn0)
F = subs(F,'x13',Dn1)
F = subs(F,'x14',NT)
F = subs(F,'x15',NT2)


%%
X7=M1-x5-x6-x4;
X11=M2-x9-x10-x8;
X13=M3-x12;
F=simplify(subs(F,'x7',X7))
F=simplify(subs(F,'x11',X11))
F=simplify(subs(F,'x13',X13))

X4=solve(F(4),x4)
F=simplify(subs(F,'x4',X4))
X5=solve(F(5),x5)
F=simplify(subs(F,'x5',X5))
X6=solve(F(6),x6)
F=simplify(subs(F,'x6',X6))

X8=solve(F(8),x8)
F=simplify(subs(F,'x8',X8))
X9=solve(F(9),x9)
F=simplify(subs(F,'x9',X9))
X10=solve(F(10),x10)
F=simplify(subs(F,'x10',X10))

X12=solve(F(12),x12)
F=simplify(subs(F,'x12',X12))
X15=solve(F(15),x15)
F=simplify(subs(F,'x15',X15))
X14=solve(F(14),x14)
F=simplify(subs(F,'x14',X14));

%%% TRY BERTINI

[f1,gg]=numden(F(1));
[f2,gg]=numden(F(2));
[f3,gg]=numden(F(3)); 
f1=expand(f1);
f2=expand(f2);
f3=expand(f3);
f=fopen('function32','w');
fprintf(f,'variable_group x1,x2,x3;\n');
fprintf(f,'function f1, f2, f3;\n');
fprintf(f,['f1=' char(f1) '; \n']);
fprintf(f,['f2=' char(f2) '; \n']);
fprintf(f,['f3=' char(f3) '; \n']);
fprintf(f,'END;');