using ApproxFun, Interact


# -------------- Symbolic ------------------
using Reduce
@force using Reduce.Algebra

function subs(old, new, Eq)
    Algebra.sub(old => new, Eq) |> mat
end
function solve(Eq, var)
    Algebra.solve(:(0  == $Eq),var)
end


fO  = :( -O*d + b*Da + 2*ep*Da - O^2*a*Di)
fDi = :( -Di*gama + ep*Da + th*Dm + (-1/2)*O^2*a*Di)
fDa = :( -ep*Da + (1/2)*O^2*a*Di)
fDm = :( Di*gama - th*Dm)
con = :( 1 -Di -Da -Dm)
@vars O Di Da Dm a ep b gama th d




Dms = solve(fDm, Dm)
fDos = subs(Dm, Dms[1].args[3], fDi)






# final subs for fO
fO
Das = solve(:($fDa ==0), Da)

fO = subs(Da, Das[1].args[3],fO)

# Dms = solve(:($fDm==0), Dm)
# cons1 = subs(Da, Das[1].args[3], con)
# cons2 = subs(Dm, Dms[1].args[3],cons1)
# solve(:($cons2 ==0), Di)

# nl =:((a*(b+ep)/(ep) -2)*(1-Di-gama*Di/th)*2ep/a -d*O)
# solve(nl,Di)

Dis = :(d*a*th*O/(2*(th-gama)*((b+ep)*a-2ep)))

final = subs(Di, Dis ,fO)

O_rslt= solve(:($final ==0), O)







a, ep, b, gama, th, d  =[1,0.1,1,1,1,1]

O1 = (2 * sqrt(((((-a * b * ep * gama + a * b * ep * th) - a * ep ^ 2 * gama) + a * ep ^ 2 * th + 2 * ep ^ 2 * gama) - 2 * ep ^ 2 * th) / (b * th))) / abs(a)

O2 = (-2 * sqrt(((((-a * b * ep * gama + a * b * ep * th) - a * ep ^ 2 * gama) + a * ep ^ 2 * th + 2 * ep ^ 2 * gama) - 2 * ep ^ 2 * th) / (b * th))) / abs(a)
