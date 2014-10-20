    // Peclet number
    volScalarField Pec  = phase1.d()*mag(U1-U2)/D2;
    
    // Reynold's number
    volScalarField Re   = 
    (
        phase1.d()*mag(U1-U2)*rho2/phase2.turbulence().muEff()
    );
    
    // Mass transfer coefficient
    volScalarField kCH4 = 
    (
        (1 + cbrt(1 + Pec)*(1 + 0.096*cbrt(Re)/(1 + 7*pow(Re+SMALL,-2))))
      * D2/(phase1.d())
    );
    
    // Interfacial area
    volScalarField aInter = 6*phase1/(phase1.d());
    
    // Absorption rate coefficient
    volScalarField kAbs = aInter*kCH4*rho2;
    
    // Mole fraction of methane in phase 1
    volScalarField omegaSatByOmega1 = 1/(omega1 + 16.042/28.96*(1-omega1));
    // Equilibrium mole fraction of methane in phase 2
    omegaSatByOmega1 *= p/kHe*kSalt;
    // Convert mole fraction to mass fraction
    omegaSatByOmega1 /= 
    (
        omegaSatByOmega1*omega1+(1-omegaSatByOmega1*omega1)*18.28/16.042
    );


    fvScalarMatrix omega1Eqn
    (
        fvm::ddt(alpha1, rho1, omega1) + fvm::div(alphaRhoPhi1, omega1)
      - fvm::Sp(contErr1, omega1)
        ==
        fvm::laplacian(alpha1*rho1*D1, omega1)
      - fvm::Sp(kAbs*omegaSatByOmega1,omega1)
      + kAbs*omega2
    );
    
    omega1Eqn.relax();
    omega1Eqn.solve();
    omega1 = min(1.0, max(0.0, omega1));
    
    
    fvScalarMatrix omega2Eqn
    (
        fvm::ddt(alpha2, rho2, omega2) + fvm::div(alphaRhoPhi2, omega2)
      - fvm::Sp(contErr2, omega2)
        ==
        fvm::laplacian(alpha2*rho2*D2, omega2)
      - fvm::Sp(kAbs, omega2)
      + kAbs*omegaSatByOmega1*omega1
    );
    
    omega2Eqn.relax();
    omega2Eqn.solve();
    omega2 = min(1.0, max(0.0, omega2));
    
    absorption = kAbs*(omegaSatByOmega1*omega1 - omega2);
