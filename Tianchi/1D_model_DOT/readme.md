# Single gene(Oct4) regulation.
## CRN model
```julia
rn1d = @reaction_network begin
    (a, ep),           2O + Di ↔ Da
     b,                Da → Da + O
    (gama,th),         Di ↔ Dm
     d,                O → ∅
 end a ep b gama th d
```

Description:

1. Oct4 binds to its promoter $D_i$ and promoter will be activated ($D_a$)
2. Activaed $D_a$ will produce Oct4 protein O and with some leakage rate
