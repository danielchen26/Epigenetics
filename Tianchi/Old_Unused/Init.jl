# Initial condition

mutable struct SType{T} <: DEDataVector{T}
    x::Array{T,1}
    p::Array{T,1}
end

function crn_init(DO_tot,DT_tot, i)
    # Give initial TET(u00[2]) some values
    u00 = zeros(15); u00[2] = 0.5;
    # Distribute the promotors for O, T, N for the same conservation class
    u00[5]  = DO_tot[i][1]; u00[6] =DO_tot[i][2]; u00[7] =DO_tot[i][3]; u00[8] =DO_tot[i][4];
    u00[10] = DT_tot[i][1]; u00[11] = DT_tot[i][2]; u00[12] = DT_tot[i][3]; u00[13] = DT_tot[i][4];
    u00[14] = 1;
    # kinectic params, keep constant ratio of 3 Kds
    p00 = [0.3, 0.2, 0.1, 1, 1, 1000, 1000, 1.0, 1.0, 1.0, 1, 0., 0., 0.05, 0.]
    u0  = SType(u00, p00)
    p = p00
    return u0, p
end
