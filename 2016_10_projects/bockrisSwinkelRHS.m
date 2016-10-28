function RHSreproc = bockrisSwinkelRHS(theta,c,f,DeltaG0m)
    T = 298.15;
    import jlh.Constants;
    RHSreproc = c/55.5 .* exp(f.*theta)*exp(-DeltaG0m/(jlh.Constants.R*T));
end
    