using Reduce


fDn0  = :(Ko*a0*Dn1 - O*a0*Dn0)
fDn1  = :(-Ko*a0*Dn1 + O*a0*Dn0)
con = :(1+ (-Dn0)+ (-Dn1))
X = [fDn0; fDn1; con]

Dn0_s =Algebra.solve(:(0  == $con),:Dn0)
X = Algebra.sub(:Dn0 => Dn0_s[1].args[2],X) |> mat

Dn1_s =Algebra.solve(:(0  == $(X[2])),:Dn1)
@show X = Algebra.sub(:Dn1 => Dn1_s[1].args[2],X) |> mat

X[2]






#  ------- paper reduce model -------
fN    = :(-N*delta + alphaN*Dn1 - N*T*a1 + NT*Kd*a1 + m1)
fT    = :(-T*delta + alphaT*Dt11 - N*T*a1 + NT*Kd*a1 + m2)
fO    = :(-O*delta + alphaO*Do11 + Ko*ao*Do01 + Ko*ao*Dt01 + Ko*ao*Dn1 + Ko*ao*Do11 + Ko*ao*Dt11 - O*ao*Dn0 - O*ao*Do00 - O*ao*Dt00 - O*ao*Do10 - O*ao*Dt10 + m3)
fDo00 = :(Ko*ao*Do01 - NT2*a_nt*Do00 - O*ao*Do00 + a_nt*K_nt*Do10)
fDo10 = :(Ko*ao*Do11 + NT2*a_nt*Do00 - O*ao*Do10 - a_nt*K_nt*Do10)
fDo01 = :(-Ko*ao*Do01 - NT2*a_nt*Do01 + O*ao*Do00 + a_nt*K_nt*Do11)
fDo11 = :(-Ko*ao*Do11 + NT2*a_nt*Do01 + O*ao*Do10 - a_nt*K_nt*Do11)
fDt00 = :(Ko*ao*Dt01 - NT2*a_nt*Dt00 - O*ao*Dt00 + a_nt*K_nt*Dt10)
fDt10 = :(Ko*ao*Dt11 + NT2*a_nt*Dt00 - O*ao*Dt10 - a_nt*K_nt*Dt10)
fDt01 = :(-Ko*ao*Dt01 - NT2*a_nt*Dt01 + O*ao*Dt00 + a_nt*K_nt*Dt11)
fDt11 = :(-Ko*ao*Dt11 + NT2*a_nt*Dt01 + O*ao*Dt10 - a_nt*K_nt*Dt11)
fDn0  = :(Ko*ao*Dn1 - O*ao*Dn0)
fDn1  = :(-Ko*ao*Dn1 + O*ao*Dn0)
fNT   = :(-NT*beta - d*NT^2 + 2*d*NT2 + N*T*a1 - NT*Kd*a1)
fNT2  = :((1/2)*d*NT^2 - d*NT2 - beta*NT2 - NT2*a_nt*Do00 - NT2*a_nt*Dt00 - NT2*a_nt*Do01 - NT2*a_nt*Dt01 + a_nt*K_nt*Do10 + a_nt*K_nt*Dt10 + a_nt*K_nt*Do11 + a_nt*K_nt*Dt11)



Eq = [fN; fT; fO; fDo00; fDo10; fDo01; fDo11; fDt00; fDt10; fDt01; fDt11; fDn0; fDn1; fNT; fNT2]


O11 = :(1 - Do00 - Do10 - Do01)
T11 = :(1 - Dt00 - Dt10 - Dt01)
N1  = :(1 - Dn0)

function subs(old, new, Eq)
    Algebra.sub(old => new, Eq) |> mat
end
function solve(Eq, var)
    Algebra.solve(:(0  == $Eq),var)
end

Eq = subs(:Do11, O11, Eq)
Eq = subs(:Dt11, T11, Eq)
Eq = subs(:Dn1, N1, Eq)

Do00_s = solve(Eq[4], :Do00)
Eq = subs(:Do00, Do00_s[1].args[3], Eq)
Do10_s = solve(Eq[5], :Do10)
Eq = subs(:Do10, Do10_s[1].args[2], Eq)
Do01_s = solve(Eq[6], :Do01)
Eq = subs(:Do01, Do01_s[1].args[3], Eq)

Dt00_s = solve(Eq[8], :Dt00)
Eq = subs(:Dt00, Dt00_s[1].args[3], Eq)
Dt10_s = solve(Eq[9], :Dt10)
Eq = subs(:Dt10, Dt10_s[1].args[3], Eq)
Dt01_s = solve(Eq[10], :Dt01)
Eq = subs(:Dt01, Dt01_s[1].args[3], Eq)

Dn0_s = solve(Eq[12], :Dn0)
Eq = subs(:Dn0, Dn0_s[1].args[3], Eq)
