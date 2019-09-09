%% Steady states solutions
% Variables and Parameters
syms Dn1 Dn0 N T O NT NT2 Do00 Do01 Do10 Do11 Dt00 Dt01 Dt10 Dt11 Dn0 Dn1 NT NT2
params = {0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0., 0., 0.05, 0.};
[K_o, K_nt, Kd, a1, d, a_nt, a_0, alpha_t, alpha_o, alpha_n, delta, beta, m1, m2, m3] =  params{:}



% make vars and params all symbolic %option
syms Dn1 Dn0 N T O NT NT2 Do00 Do01 Do10 Do11 Dt00 Dt01 Dt10 Dt11 Dn0 Dn1 NT NT2 K_o  K_nt  Kd  a1  d  a_nt  a_0  alpha_t  alpha_o alpha_n delta beta m1 m2 m3
% Equations for SSS
fN   = Dn1*alpha_n + Kd*NT*a1 - N*T*a1 - N*delta + m1
fT   = Dt11*alpha_t + Kd*NT*a1 - N*T*a1 - T*delta + m2
fO   = -Dn0*O*a_0 + Dn1*K_o*a_0 - Do00*O*a_0 + Do01*K_o*a_0 - Do10*O*a_0 + Do11*K_o*a_0 + Do11*alpha_o - Dt00*O*a_0 + Dt01*K_o*a_0 - Dt10*O*a_0 + Dt11*K_o*a_0 - O*delta + m3
fDo00 = -Do00*NT2*a_nt - Do00*O*a_0 + Do01*K_o*a_0 + Do10*K_nt*a_nt
fDo10 = Do00*NT2*a_nt - Do10*K_nt*a_nt - Do10*O*a_0 + Do11*K_o*a_0
fDo01 = Do00*O*a_0 - Do01*K_o*a_0 - Do01*NT2*a_nt + Do11*K_nt*a_nt
fDo11 = Do01*NT2*a_nt + Do10*O*a_0 - Do11*K_nt*a_nt - Do11*K_o*a_0
fDt00 = -Dt00*NT2*a_nt - Dt00*O*a_0 + Dt01*K_o*a_0 + Dt10*K_nt*a_nt
fDt10 = Dt00*NT2*a_nt - Dt10*K_nt*a_nt - Dt10*O*a_0 + Dt11*K_o*a_0
fDt01 = Dt00*O*a_0 - Dt01*K_o*a_0 - Dt01*NT2*a_nt + Dt11*K_nt*a_nt
fDt11 = Dt01*NT2*a_nt + Dt10*O*a_0 - Dt11*K_nt*a_nt - Dt11*K_o*a_0
fDn0  = -Dn0*O*a_0 + Dn1*K_o*a_0
fDn1  = Dn0*O*a_0 - Dn1*K_o*a_0
fNT  = -Kd*NT*a1 + N*T*a1 - NT^2*d - NT*beta + 2*NT2*d
fNT2 = -Do00*NT2*a_nt - Do01*NT2*a_nt + Do10*K_nt*a_nt + Do11*K_nt*a_nt - Dt00*NT2*a_nt - Dt01*NT2*a_nt + Dt10*K_nt*a_nt + Dt11*K_nt*a_nt + 0.5*NT^2*d - NT2*beta - NT2*d
ConsDo = 1 - Do10 - Do01 - Do00 - Do11;
ConsDt = 1 - Dt10 - Dt01 - Dt00 - Dt11;
ConsDn = 1 - Dn0 - Dn1;



% Reduce DF to 3 (NTO)
F2 =[fN; fT; fO; fDo00; fDo10; fDo01; fDo11; fDt00; fDt10;  fDt01;  fDt11;  fDn0;  fDn1;  fNT;  fNT2; ConsDn;  ConsDo; ConsDt];
% Conservation law
Do11_s=solve(ConsDo,Do11)
F2=simplify(subs(F2,'Do11',Do11_s))
Dt11_s=solve(ConsDt,Dt11)
F2=simplify(subs(F2,'Dt11',Dt11_s))
Dn1_s=solve(ConsDn,Dn1)
F2=simplify(subs(F2,'Dn1',Dn1_s))

% Do promoter elimination
Do00_s=solve(F2(4),Do00)
F2=simplify(subs(F2,'Do00',Do00_s))
Do10_s=solve(F2(5),Do10)
F2=simplify(subs(F2,'Do10',Do10_s))
Do01_s=solve(F2(6),Do01)
F2=simplify(subs(F2,'Do01',Do01_s))

% Dt promoter elimination
Dt00_s=solve(F2(8),Dt00)
F2=simplify(subs(F2,'Dt00',Dt00_s))
Dt10_s=solve(F2(9),Dt10)
F2=simplify(subs(F2,'Dt10',Dt10_s))
Dt01_s=solve(F2(10),Dt01)
F2=simplify(subs(F2,'Dt01',Dt01_s))


% Dn promoter elimination
Dn0_s=solve(F2(12),Dn0)
F2=simplify(subs(F2,'Dn0',Dn0_s))

% NT, NT2 promoter elimination
NT2_s=solve(F2(15),NT2)
F2=simplify(subs(F2,'NT2',NT2_s))
NT_s=solve(F2(14),NT)

F2=simplify(subs(xp,'NT',NT_s))

% with all para symbolic
F2_1=simplify(subs(F2,'NT',NT_s(1)))
F2_2=simplify(subs(F2,'NT',NT_s(2)))


%% Homotopy Continuation find SSS
[f1,gg]=numden(F2_1(1));
[f2,gg]=numden(F2_1(2));
[f3,gg]=numden(F2_1(3)); 
f1=expand(f1)
f2=expand(f2)
f3=expand(f3)
[f1;f2;f3]  


f=fopen('~/Desktop/bertini_test/tc_crn_ode.txt','w');
fprintf(f,'variable_group N,T,O;\n');
fprintf(f,'function f1, f2, f3;\n');
fprintf(f,['f1=' char(f1) '; \n']);
fprintf(f,['f2=' char(f2) '; \n']);
fprintf(f,['f3=' char(f3) '; \n']);
fprintf(f,'END;');
fclose(f);

system('/usr/local/bin/bertini-serial ~/Desktop/bertini_test/tc_crn_ode.txt')

f22=fopen('~/Desktop/bertini_test/real_finite_solutions','r');
type '~/Desktop/bertini_test/real_finite_solutions'

%% jacobian 
rd = F2(1:3)
jac = jacobian(rd, [N,T,O])

% 3 steady states
SD1 = {0.304,0.181,0.131}
SD2 = {0.,0.05,0.}
SD3 = {0.6954,0.7349,0.6849}

[N,T,O] = SD1 {:}
j1 = subs(jac)
eigen1 = vpa(eig(j1),2)

[N,T,O] = SD2 {:}
j2 = subs(jac)
eigen2 = vpa(eig(j2),2)

[N,T,O] = SD3 {:}
j3 = subs(jac)
eigen3 = vpa(eig(j3),2)

eigen1,eigen2,eigen3

%% ODE solve

syms t x m1    
params = {0.3, 0.2, 0.1, 1.0, 1.0, 1000.0, 1000.0, 1.0, 1.0, 1.0, 1.0, 1e-6, 1e-6, 0.05, 1e-6};
[K_o, K_nt, Kd, a1, d, a_nt, a_0, alpha_t, alpha_o, alpha_n, delta, beta, m1, m2, m3] =  params{:}


rd = @(t,x,m1)[m1 - x(1)*delta - x(1)*x(2)*a1 + (x(3)*alpha_n)/(K_o + x(3)) - (Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d);
    -(K_nt*K_o*x(2)*delta - K_nt*x(3)*m2 - K_nt*K_o*m2 + K_nt*x(3)*x(2)*delta + K_nt*K_o*x(1)*x(2)*a1 + K_nt*x(1)*x(3)*x(2)*a1 - (x(3)*alpha_t*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (x(3)*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (Kd*x(3)*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_o*x(2)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(3)*x(2)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*Kd*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_nt*K_o*Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_nt*Kd*x(3)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_o*x(1)*x(2)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(1)*x(3)*x(2)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + x(3)));
    (K_nt*K_o*m3 + K_nt*x(3)*m3 - K_nt*x(3)^2*delta - K_nt*K_o*x(3)*delta + (x(3)*alpha_o*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(3)*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (x(3)^2*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*x(3)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + x(3)))];
        
syms N T O
rd_NTO = [m1 - N*delta - N*T*a1 + (O*alpha_n)/(K_o + O) - (Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d);
    -(K_nt*K_o*T*delta - K_nt*O*m2 - K_nt*K_o*m2 + K_nt*O*T*delta + K_nt*K_o*N*T*a1 + K_nt*N*O*T*a1 - (O*alpha_t*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (O*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (Kd*O*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_o*T*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (O*T*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*Kd*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_nt*K_o*Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_nt*Kd*O*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_o*N*T*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (N*O*T*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + O));
    (K_nt*K_o*m3 + K_nt*O*m3 - K_nt*O^2*delta - K_nt*K_o*O*delta + (O*alpha_o*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (O*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (O^2*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*O*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + O))]

jac = vpa(jacobian(rd_NTO, [N,T,O]),2)
eigen = vpa(eig(jac),2)

SD1 = {0.304,0.181,0.131};[N,T,O] = SD1 {:};subs(jac)
SD2 = {0.,0.05,0.};[N,T,O] = SD2{:};subs(eigen)
SD3 = {0.6954,0.7349,0.6849};[N,T,O] = SD3{:};subs(eigen)



u0= [0.1,0.5,0.1]
[t,xp] = ode45(@(t,x) rd(t,x), [0 1e3], u0 );


figure()
plot(t,xp,'LineWidth',2)
legend('N','T','O')

%---------------------------------%---------------------------------%---------------------------------%---------------------------------
%% make for loop for params sweep for oct4 overexpress
params = {0.3, 0.2, 0.1, 1.0, 1.0, 1000.0, 1000.0, 1.0, 1.0, 1.0, 1.0, 1e-6, 1e-6, 0.05, 1e-6};
[K_o, K_nt, Kd, a1, d, a_nt, a_0, alpha_t, alpha_o, alpha_n, delta, beta, m1, m2, m3] =  params{:}

syms t x m1 m2 m3 N T O
rd_NTO = [m1 - N*delta - N*T*a1 + (O*alpha_n)/(K_o + O) - (Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d);
    -(K_nt*K_o*T*delta - K_nt*O*m2 - K_nt*K_o*m2 + K_nt*O*T*delta + K_nt*K_o*N*T*a1 + K_nt*N*O*T*a1 - (O*alpha_t*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (O*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (Kd*O*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_o*T*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (O*T*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*Kd*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_nt*K_o*Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_nt*Kd*O*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_o*N*T*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (N*O*T*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + O));
    (K_nt*K_o*m3 + K_nt*O*m3 - K_nt*O^2*delta - K_nt*K_o*O*delta + (O*alpha_o*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (O*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (O^2*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*O*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*N*T*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + O))]
jac = vpa(jacobian(rd_NTO, [N,T,O]),2);
set(0,'DefaultFigureVisible','off')
for K_o = 0.3:0.05:0.5
    for K_nt = 0.2:0.05:0.5
        for Kd = 0.1:0.05:0.5
            
            fprintf('K_o := %0.2f   K_nt :=  %0.2f   Kd :=  %0.2f\n ', K_o, K_nt, Kd)
            
            rd = @(t,x,m1,m2,m3)[m1 - x(1)*delta - x(1)*x(2)*a1 + (x(3)*alpha_n)/(K_o + x(3)) - (Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d);
            -(K_nt*K_o*x(2)*delta - K_nt*x(3)*m2 - K_nt*K_o*m2 + K_nt*x(3)*x(2)*delta + K_nt*K_o*x(1)*x(2)*a1 + K_nt*x(1)*x(3)*x(2)*a1 - (x(3)*alpha_t*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (x(3)*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (Kd*x(3)*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_o*x(2)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(3)*x(2)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*Kd*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_nt*K_o*Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_nt*Kd*x(3)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_o*x(1)*x(2)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(1)*x(3)*x(2)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + x(3)));
            (K_nt*K_o*m3 + K_nt*x(3)*m3 - K_nt*x(3)^2*delta - K_nt*K_o*x(3)*delta + (x(3)*alpha_o*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(3)*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (x(3)^2*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*x(3)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + x(3)))];


            fig = figure('visible','off');
            set(fig,'Position',[100 100 1200 1300]) % for dynamics
            for i = 1:12             
                m1 = 0; m2 = 0.05; m3 = 0;
                u0= rand(3,1);
                [t,xp] = ode45(@(t,x) rd(t,x,m1,m2,m3), [0 50], u0 );

                % calculate jacobian

                N = xp(end,1); T = xp(end,2); O = xp(end,3);
                jac_num = subs(jac);
                eigen = vpa(eig(jac_num),2);


                m1 = 0.06; m3 = 0.36;
                u0_t_on = xp(end,:);
                [t_on,xp_on] = ode45(@(t,x) rd(t,x,m1,m2,m3), [0 5], u0_t_on );

                m1 = 0; m3 = 0;
                u0_t_off = xp_on(end,:);
                [t_off,xp_off] = ode45(@(t,x) rd(t,x,m1,m2,m3), [0 50], u0_t_off );


                tf = [t; t_on + t(end); t_off + t(end) + t_on(end)];
                xpf = [xp; xp_on; xp_off];

                subplot(4,3,i)
                plot(tf, xpf,'LineWidth',2)
                ylim([0,1.5])
                title(sprintf("%0.2f %0.2f %0.2f \n %0.2f %0.2f %0.2f",K_o, K_nt, Kd, eigen(1), eigen(2), eigen(3)))
            end
            legend('N','T','O')
            suptitle('Varying K_o, K_{nt}, Kd. Each set: 12 random initials')
            saveas(gca, sprintf('~/Desktop/epi_plots/params_explore/same_param/K_o %0.2f K_nt %0.2f Kd %0.2f.png',K_o, K_nt, Kd), 'png')
        end
    end
end

%% Same Randomized initial with different parameter set comparison
rand_initial_set = rand(50,3);

for i = 1:size(rand_initial_set,1)
    
    u0 = rand_initial_set(i,:);
    % now iterate different parameters with K_o, K_nt, Kd
    fig = figure('visible','off');
%     set(fig,'Position',[100 100 800 1300]) % for eigen spectral
    set(fig,'Position',[100 100 1200 1300]) % for dynamics
    n = 1;
    for K_o = 0.3:0.05:0.4
        for K_nt = 0.2:0.05:0.3
            for Kd = 0.1:0.05:0.2
                fprintf('K_o := %0.2f   K_nt :=  %0.2f   Kd :=  %0.2f\n ', K_o, K_nt, Kd)
                
                rd = @(t,x,m1,m2,m3)[m1 - x(1)*delta - x(1)*x(2)*a1 + (x(3)*alpha_n)/(K_o + x(3)) - (Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d);
                -(K_nt*K_o*x(2)*delta - K_nt*x(3)*m2 - K_nt*K_o*m2 + K_nt*x(3)*x(2)*delta + K_nt*K_o*x(1)*x(2)*a1 + K_nt*x(1)*x(3)*x(2)*a1 - (x(3)*alpha_t*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (x(3)*m2*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (Kd*x(3)*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_o*x(2)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(3)*x(2)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*Kd*a1*(beta + d)^2*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^3)/(16*beta^3*d^2) + (K_nt*K_o*Kd*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_nt*Kd*x(3)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2)))/(2*beta*d) + (K_o*x(1)*x(2)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(1)*x(3)*x(2)*a1*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + x(3)));
                (K_nt*K_o*m3 + K_nt*x(3)*m3 - K_nt*x(3)^2*delta - K_nt*K_o*x(3)*delta + (x(3)*alpha_o*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (K_o*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) + (x(3)*m3*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (x(3)^2*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d) - (K_o*x(3)*delta*(beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))/((K_nt + ((beta + d)*(beta + Kd*a1 - ((Kd^2*a1^2*beta + d*Kd^2*a1^2 + 2*Kd*a1*beta^2 + 2*d*Kd*a1*beta + 4*x(1)*x(2)*d*a1*beta + beta^3 + d*beta^2)/(beta + d))^(1/2))^2)/(8*beta^2*d))*(K_o + x(3)))];
                

                m1 = 0; m2 = 0.05; m3 = 0;
%                 u0= rand(3,1);
                [t,xp] = ode45(@(t,x) rd(t,x,m1,m2,m3), [0 50], u0 );

                % calculate jacobian

                N = xp(end,1); T = xp(end,2); O = xp(end,3);
                jac_num = subs(jac);
                eigen = vpa(eig(jac_num),2);
                
%                 % plot eigen spectral
%                 subplot(9,3,n)
%                 n = n + 1
%                 bar(eigen)
%                 title(sprintf("%0.2f %0.2f %0.2f \n %0.2f %0.2f %0.2f",K_o, K_nt, Kd, eigen(1), eigen(2), eigen(3)))

                % plot dynamics
                m1 = 0.06; m3 = 0.36;
                u0_t_on = xp(end,:);
                [t_on,xp_on] = ode45(@(t,x) rd(t,x,m1,m2,m3), [0 5], u0_t_on );

                m1 = 0; m3 = 0;
                u0_t_off = xp_on(end,:);
                [t_off,xp_off] = ode45(@(t,x) rd(t,x,m1,m2,m3), [0 50], u0_t_off );


                tf = [t; t_on + t(end); t_off + t(end) + t_on(end)];
                xpf = [xp; xp_on; xp_off];

                subplot(9,3,n)
                n = n + 1;
                plot(tf, xpf,'LineWidth',2)
                
                ylim([0,1.5])
%                 title(sprintf("%0.2f %0.2f %0.2f",eigen(1), eigen(2), eigen(3)))
                title(sprintf("%0.2f %0.2f %0.2f \n %0.2f %0.2f %0.2f",K_o, K_nt, Kd, eigen(1), eigen(2), eigen(3)))
            end
        end
    end
    legend('N','T','O')
%     saveas(gca, sprintf('~/Desktop/epi_plots/params_explore/same_init/%d.png',i), 'png')
    sprintf('current cycle %d',i)
end


    
    
    










