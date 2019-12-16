# ================ CRN ======================
using DiffEqBiological, Latexify, ParameterizedFunctions

# === version 2. NT control the demethylation
Demethy_TF_MC = @reaction_network begin
    # ============== 💠 ===============
    # N T complex
    (a1, Kd*a1),              N + T ↔ NT
    # --- Promoter binding
    # Oct4
    (a_nt, K_nt*a_nt),        NT + Do00 ↔ Do10
    (a_nt, K_nt*a_nt),        NT + Do01 ↔ Do11
    (a_O, K_O*a_O),              O + Do00 ↔ Do01
    (a_O, K_O*a_O),              O + Do10 ↔ Do11
    # TET
    (a_nt, K_nt*a_nt),        NT + Dt00 ↔ Dt10
    (a_nt, K_nt*a_nt),        NT + Dt01 ↔ Dt11
    (a_O, K_O*a_O),              O + Dt00 ↔ Dt01
    (a_O, K_O*a_O),              O + Dt10 ↔ Dt11
    # Nanog
    (a_O, K_O*a_O),              O + Dn0 ↔ Dn1

    # --- Protein production
    α_T,                   Dt11 → Dt11 + T
    α_O,                   Do11 → Do11 + O
    α_N,                   Dn1 → Dn1 + N

    # --- Dilution and Degradation
    (δ,δ,δ),           (N, T, O) → ∅

    # ============   🔺 ===============
    # # Oct4 de-Methylation cycle  ---- 💚NT or NT2?
    # a_dn,                Do00 → D5mc                # D + DNMT ↔ C₁ → D5mc + DNMT
    # kh,                  NT + D5mc → D5hmc + NT     # NT oxidize 5mc -> 5hmC
    # beta,                D5hmc →  Do00              # 5hmC -> C by
    # ============ ⟺ ===================
    γ,             Do00 → Dm   # simplification of the 🔺
    θ,             Dm + NT ⇀ Do00 + NT
    # ---- NTO rate control m1 m2 m3-----
    m1,              ∅ → N
    m2,              ∅ → T
    m3,              ∅ → O
end K_O K_nt Kd a1 a_nt a_O α_T α_O α_N δ γ θ m1 m2 m3





# 3 Genes: {N,T,O}, but with slow dynamics Dm. We could defind the BOA in this 4d dimensions
reduced_ODE_4d = @ode_def_bare begin
    dN  = m1 - N*δ + (O*α_N)/(K_O + O)
    dT  = m2 - T*δ + (N*O*T*α_T)/((N*T + K_nt*Kd)*(K_O + O))
    dO  = -(K_nt*Kd*O^2*δ + N*O^2*T*δ - K_O*K_nt*Kd*m3 - N*O*T*α_O - K_nt*Kd*O*m3 - K_O*N*T*m3 - N*O*T*m3 + K_O*N*O*T*δ + K_O*K_nt*Kd*O*δ + Dm*N*O*T*α_O)/((N*T + K_nt*Kd)*(K_O + O))
    dDm = -(Dm*K_O*K_nt*Kd^2*γ - K_O*K_nt*Kd^2*γ + Dm*K_O*N^2*T^2*θ + Dm*N^2*O*T^2*θ + Dm*K_O*K_nt*Kd*N*T*θ + Dm*K_nt*Kd*N*O*T*θ)/(Kd*(N*T + K_nt*Kd)*(K_O + O))
end K_O K_nt Kd a1 a_nt a_O α_T α_O α_N δ γ θ m1 m2 m3











 # -------------------- Latex in the paper --------------------
\begin{align}
\frac{dN(t)}{dt} =&  - a N T + a k_{d} [NT]  + {\alpha}_N D^{N}_{1} - \delta N + m_{1} \\
\frac{dT(t)}{dt} =&  - a N T + a k_{d} [NT] + {\alpha}_T D^{T}_{11} - \delta T + m_{2} \\
\frac{d[NT](t)}{dt} =& a N T - a k_{d} [NT] - a_{nt} [NT] D^{O}_{00} + K_{nt} a_{nt} D^{O}_{10}  - a_{nt} [NT] D^{O}_{01} + K_{nt} a_{nt} D^{O}_{11}  \\
                  & - a_{nt} [NT] D^{T}_{00} + K_{nt} a_{nt} D^{T}_{10}  - a_{nt} [NT] D^{T}_{01} + K_{nt} a_{nt} D^{T}_{11}  \\
\frac{dD^{O}_{00}(t)}{dt} =&  - a_{nt} [NT] D^{O}_{00} + K_{nt} a_{nt} D^{O}_{10}  - a_{O} O D^{O}_{00} + K_{O} a_{O} D^{O}_{01} - \gamma D^{O}_{00} + \theta D_{m} [NT] \\
\frac{dD^{O}_{10}(t)}{dt} =& a_{nt} [NT] D^{O}_{00} - K_{nt} a_{nt} D^{O}_{10}  - a_{O} O D^{O}_{10} + K_{O} a_{O} D^{O}_{11} \\
\frac{dD^{O}_{01}(t)}{dt} =&  - a_{nt} [NT] D^{O}_{01} + K_{nt} a_{nt} D^{O}_{11}  + a_{O} O D^{O}_{00} - K_{O} a_{O} D^{O}_{01} \\
\frac{dD^{O}_{11}(t)}{dt} =& a_{nt} [NT] D^{O}_{01} - K_{nt} a_{nt} D^{O}_{11}  + a_{O} O D^{O}_{10} - K_{O} a_{O} D^{O}_{11} \\
\frac{dO(t)}{dt} =&  - a_{O} O D^{O}_{00} + K_{O} a_{O} D^{O}_{01} - a_{O} O D^{O}_{10} + K_{O} a_{O} D^{O}_{11} - a_{O} O D^{T}_{00} + K_{O} a_{O} D^{T}_{01}\\
                  & - a_{O} O D^{T}_{10} + K_{O} a_{O} D^{T}_{11} - a_{O} O D^{N}_{0} + K_{O} a_{O} D^{N}_{1} + {\alpha}_O D^{O}_{11} - \delta O + m_{3} \\
\frac{dD^{T}_{00}(t)}{dt} =&  - a_{nt} [NT] D^{T}_{00} + K_{nt} a_{nt} D^{T}_{10}  - a_{O} O D^{T}_{00} + K_{O} a_{O} D^{T}_{01} \\
\frac{dD^{T}_{10}(t)}{dt} =& a_{nt} [NT] D^{T}_{00} - K_{nt} a_{nt} D^{T}_{10}  - a_{O} O D^{T}_{10} + K_{O} a_{O} D^{T}_{11} \\
\frac{dD^{T}_{01}(t)}{dt} =&  - a_{nt} [NT] D^{T}_{01} + K_{nt} a_{nt} D^{T}_{11}  + a_{O} O D^{T}_{00} - K_{O} a_{O} D^{T}_{01} \\
\frac{dD^{T}_{11}(t)}{dt} =& a_{nt} [NT] D^{T}_{01} - K_{nt} a_{nt} D^{T}_{11}  + a_{O} O D^{T}_{10} - K_{O} a_{O} D^{T}_{11} \\
\frac{dD^{N}_{0}(t)}{dt} =&  - a_{O} O D^{N}_{0} + K_{O} a_{O} D^{N}_{1} \\
\frac{dD^{N}_{1}(t)}{dt} =& a_{O} O D^{N}_{0} - K_{O} a_{O} D^{N}_{1} \\
\frac{dD_{m}(t)}{dt} =& \gamma D^{O}_{00} - \theta D_{m} [NT]
\end{align}




\begin{tiny}
\begin{align}
\frac{dN}{dt} = & m_{1} - \delta N + \frac{\alpha_{N} O}{K_{O} + O} \\
\frac{dT}{dt} = & m_{2} - \delta T + \frac{\alpha_{T} N T O}{\left( N T + K_{nt} K_{d} \right) \left( K_{O} + O \right)} \\
\frac{dO}{dt} = & \frac{ - \left( K_{nt} K_{d} \delta O^{2} + \delta N T O^{2} - K_{O} K_{nt} K_{d} m_{3} - \alpha_{O} N T O - K_{nt} K_{d} m_{3} O - K_{O} m_{3} N T - m_{3} N T O + K_{O} \delta N T O + K_{O} K_{nt} K_{d} \delta O  +  \alpha_{O} D_{m} N T O \right)}{\left( N T + K_{nt} K_{d} \right) \left( K_{O} + O \right)} \\
\frac{dD_{m}}{dt} = & \frac{ - \left( K_{O} K_{nt} K_{d}^{2} \gamma D_{m} - K_{O} K_{nt} K_{d}^{2} \gamma + \theta  K_{O} N^{2} T^{2} D_{m} +  \theta N^{2} T^{2} O D_{m} +  \theta K_{O} K_{nt} K_{d} N T D_{m} + \theta K_{nt} K_{d} N T O D_{m} \right)}{K_{d} \left( N T + K_{nt} K_{d} \right) \left( K_{O} + O \right)}
\end{align}
\end{tiny}
