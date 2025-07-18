#* Michaels-Menten version of base_rn
const MM_rn = @reaction_network MM_rn begin
    #- Parameters to be optimized, both rate constants and DF (dimensionality factor)
    @parameters kfᴸᴬ = $(INPUT_DEFAULTS[:kfᴸᴬ]) [bounds = $(INPUT_BOUNDS[:kfᴸᴬ])] krᴸᴬ = $(INPUT_DEFAULTS[:krᴸᴬ]) [bounds = $(INPUT_BOUNDS[:krᴸᴬ])] Kmᴸᴷ = $(INPUT_DEFAULTS[:kcatᴸᴷ]) [bounds = $(INPUT_BOUNDS[:kcatᴸᴷ])] kcatᴸᴷ = $(INPUT_DEFAULTS[:kcatᴸᴷ]) [bounds = $(INPUT_BOUNDS[:kcatᴸᴷ])] Kmᴸᴾ = $(INPUT_DEFAULTS[:kcatᴸᴾ]) [bounds = $(INPUT_BOUNDS[:kcatᴸᴾ])] kcatᴸᴾ = $(INPUT_DEFAULTS[:kcatᴸᴾ]) [bounds = $(INPUT_BOUNDS[:kcatᴸᴾ])] kfᴬᴷ = $(INPUT_DEFAULTS[:kfᴬᴷ]) [bounds = $(INPUT_BOUNDS[:kfᴬᴷ])] krᴬᴷ = $(INPUT_DEFAULTS[:krᴬᴷ]) [bounds = $(INPUT_BOUNDS[:krᴬᴷ])] kfᴬᴾ = $(INPUT_DEFAULTS[:kfᴬᴾ]) [bounds = $(INPUT_BOUNDS[:kfᴬᴾ])] krᴬᴾ = $(INPUT_DEFAULTS[:krᴬᴾ]) [bounds = $(INPUT_BOUNDS[:krᴬᴾ])] DF = $(INPUT_DEFAULTS[:DF]) [description = "Unitless dimensionality factor"; bounds = $(INPUT_BOUNDS[:DF])] 


    #- Species
    @species L(t) = $(INPUT_DEFAULTS[:L]) [description = "PIP"; bounds = $(INPUT_BOUNDS[:L])] K(t) = $(INPUT_DEFAULTS[:K]) [description = "PIP5K"; bounds = $(INPUT_BOUNDS[:K])] P(t) = $(INPUT_DEFAULTS[:P]) [description = "Synaptojanin"; bounds = $(INPUT_BOUNDS[:P])] A(t) = $(INPUT_DEFAULTS[:A]) [description = "AP2"; bounds = $(INPUT_BOUNDS[:A])] Lp(t) = 0.0 [description = "PIP2"; tunable = false] LpA(t) = 0.0 [description = "PIP2-AP2"; tunable = false] LpAK(t) = 0.0 [description = "PIP2-AP2-Kinase"; tunable = false] LpAP(t) = 0.0 [description = "PIP2-AP2-Phosphatase"; tunable = false] 

    #$ Reaction 1
    (kfᴸᴬ,krᴸᴬ), Lp + A <--> LpA # Lp binding to AP2 adaptor 

    #$ Reaction 2
    (kfᴬᴷ,krᴬᴷ), LpA + K <--> LpAK # Membrane-bound adaptor binding to kinase
    #$ Reaction 3
    (kfᴬᴾ,krᴬᴾ), LpA + P <--> LpAP # Membrane-bound adaptor binding to phosphatase 

    #$ Reaction 4a
    # (kfᴸᴷ,krᴸᴷ), L + K <--> LK # L binding to kinase
    # kcatᴸᴷ, LK --> Lp + K # L phosphorylation by kinase into Lp
    ((L * kcatᴸᴷ * K)/ (L + Kmᴸᴷ)), L + K --> Lp + K

    #$ Reaction 5a
    # (kfᴸᴾ,krᴸᴾ), Lp + P <--> LpP # Lp binding to phosphatase 
    # kcatᴸᴾ, LpP --> L + P # L dephosphorylation by phosphatase
    ((Lp * kcatᴸᴾ * P)/ (Lp + Kmᴸᴾ)), Lp + P --> L + P

    #$ Reaction 4b
    # (kfᴸᴷ*DF,krᴸᴷ), LpAK + L <--> LpAKL # 2D reaction: Membrane-bound kinase binds to L with greater affinity as determined by DF
    # kcatᴸᴷ, LpAKL --> Lp + LpAK # L phosphorylation by kinase into Lp, same as 3D: first order reactions aren't dependent on dimensionality 
    ((L * kcatᴸᴷ * LpAK)/ (L + (Kmᴸᴷ/DF))), L + LpAK --> Lp + LpAK

    #$ Reaction 5b
    # (kfᴸᴾ*DF,krᴸᴾ), Lp + LpAP <--> LpAPLp # 2D reaction: Membrane-bound phosphatase binds to Lp with greater affinity as determined by DF
    # kcatᴸᴾ, LpAPLp --> L + LpAP # L dephosphorylation by phosphatase, same as 3D: first order reactions aren't dependent on dimensionality
    #- MM approximation. Vmax = kcatᴸᴾ * LpAP
    ((Lp * kcatᴸᴾ * LpAP)/ (Lp + (Kmᴸᴾ/DF))), Lp + LpAP --> L + LpAP
end  


