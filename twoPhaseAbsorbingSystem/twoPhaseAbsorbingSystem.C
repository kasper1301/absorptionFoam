/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "twoPhaseAbsorbingSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "BlendedInterfacialModel.H"
#include "dragModel.H"
#include "virtualMassModel.H"
#include "heatTransferModel.H"
#include "liftModel.H"
#include "wallLubricationModel.H"
#include "turbulentDispersionModel.H"
#include "fvMatrix.H"
#include "surfaceInterpolate.H"
#include "MULES.H"
#include "subCycle.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcCurl.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvmDiv.H"
#include "fixedValueFvsPatchFields.H"

#include "blendingMethod.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseAbsorbingSystem::twoPhaseAbsorbingSystem
(
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(mesh),

    phase1_
    (
        *this,
        *this,
        wordList(lookup("phases"))[0]
    ),

    phase2_
    (
        *this,
        *this,
        wordList(lookup("phases"))[1]
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->calcPhi()
    ),

    dgdt_
    (
        IOobject
        (
            "dgdt",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dgdt", dimless/dimTime, 0)
    ),
    
    absorption_
    (
        IOobject
        (
            "absorption",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("absorption", dimensionSet(1, -3, -1, 0, 0), 0)
    ),
    
    wCH4sat_
    (
        IOobject
        (
            "wCH4sat",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("wCH4sat", dimless, 0)
    ),
    
    ks_
    (
        "ks",
        dimless,
        this->lookup("ks")
    ),
    
    I_
    (
        "I",
        dimless,
        this->lookup("I")
    ),
    
    kHe_
    (
        "kHe",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        this->lookup("kHe")
    ),
    
    MCH4_
    (
        "MCH4",
        dimless,
        this->lookup("MCH4")
    ),
    
    Mair_
    (
        "Mair",
        dimless,
        this->lookup("Mair")
    ),
    
    Msw_
    (
        "Msw",
        dimless,
        this->lookup("Msw")
    )
    
{
    phase2_.volScalarField::operator=(scalar(1) - phase1_);


    // Blending
    forAllConstIter(dictionary, subDict("blending"), iter)
    {
        blendingMethods_.insert
        (
            iter().dict().dictName(),
            blendingMethod::New
            (
                iter().dict(),
                wordList(lookup("phases"))
            )
        );
    }


    // Pairs

    phasePair::scalarTable sigmaTable(lookup("sigma"));
    phasePair::dictTable aspectRatioTable(lookup("aspectRatio"));

    pair_.set
    (
        new phasePair
        (
            phase1_,
            phase2_,
            g,
            sigmaTable
        )
    );

    pair1In2_.set
    (
        new orderedPhasePair
        (
            phase1_,
            phase2_,
            g,
            sigmaTable,
            aspectRatioTable
        )
    );

    pair2In1_.set
    (
        new orderedPhasePair
        (
            phase2_,
            phase1_,
            g,
            sigmaTable,
            aspectRatioTable
        )
    );


    // Models

    drag_.set
    (
        new BlendedInterfacialModel<dragModel>
        (
            lookup("drag"),
            (
                blendingMethods_.found("drag")
              ? blendingMethods_["drag"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    virtualMass_.set
    (
        new BlendedInterfacialModel<virtualMassModel>
        (
            lookup("virtualMass"),
            (
                blendingMethods_.found("virtualMass")
              ? blendingMethods_["virtualMass"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    heatTransfer_.set
    (
        new BlendedInterfacialModel<heatTransferModel>
        (
            lookup("heatTransfer"),
            (
                blendingMethods_.found("heatTransfer")
              ? blendingMethods_["heatTransfer"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    lift_.set
    (
        new BlendedInterfacialModel<liftModel>
        (
            lookup("lift"),
            (
                blendingMethods_.found("lift")
              ? blendingMethods_["lift"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    wallLubrication_.set
    (
        new BlendedInterfacialModel<wallLubricationModel>
        (
            lookup("wallLubrication"),
            (
                blendingMethods_.found("wallLubrication")
              ? blendingMethods_["wallLubrication"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );

    turbulentDispersion_.set
    (
        new BlendedInterfacialModel<turbulentDispersionModel>
        (
            lookup("turbulentDispersion"),
            (
                blendingMethods_.found("turbulentDispersion")
              ? blendingMethods_["turbulentDispersion"]
              : blendingMethods_["default"]
            ),
            pair_,
            pair1In2_,
            pair2In1_
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseAbsorbingSystem::~twoPhaseAbsorbingSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::twoPhaseAbsorbingSystem::rho() const
{
    return phase1_*phase1_.thermo().rho() + phase2_*phase2_.thermo().rho();
}

void Foam::twoPhaseAbsorbingSystem::updateAbsorption()
{
    const volScalarField& alpha1 = phase1_;
    const volScalarField& wCH4_1 = phase1_.wCH4();

    const volScalarField& rho2 = phase2_.thermo().rho();
    const volScalarField& wCH4_2 = phase2_.wCH4();
    const dimensionedScalar& DCH4_2 = phase2_.DCH4();

    const volScalarField& p = phase1_.thermo().p();

    volScalarField magUr = mag(phase1_.U() - phase2_.U());

    volScalarField Pe = phase1_.d()*magUr/DCH4_2;

    volScalarField Re = phase1_.d()*magUr/phase2_.nu();

    volScalarField kCH4
    (
        (
            1.0 + cbrt(1.0+Pe)
          * (
                1.0+0.096*cbrt(Re)/(1.0+7.0/sqr(max(Re,SMALL)))
            )
        )
      * (DCH4_2/phase1_.d())
    );

    volScalarField a(6.0*alpha1/phase1_.d());

    volScalarField xCH4_1
    (
        wCH4_1/MCH4_ /(wCH4_1/MCH4_ + (1.0-wCH4_1)/Mair_)
    );

    volScalarField xCH4_0(xCH4_1*p/kHe_);

    volScalarField xCH4sat(xCH4_0*exp(I_*ks_));

    wCH4sat_ = xCH4sat*MCH4_/(xCH4sat*MCH4_ + (1-xCH4sat*Msw_));

    absorption_ = a*kCH4*rho2*(wCH4sat_ - wCH4_2);
}

Foam::tmp<Foam::volVectorField> Foam::twoPhaseAbsorbingSystem::U() const
{
    return phase1_*phase1_.U() + phase2_*phase2_.U();
}


Foam::tmp<Foam::surfaceScalarField> Foam::twoPhaseAbsorbingSystem::calcPhi() const
{
    return
        fvc::interpolate(phase1_)*phase1_.phi()
      + fvc::interpolate(phase2_)*phase2_.phi();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseAbsorbingSystem::dragCoeff() const
{
    return drag_->K();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseAbsorbingSystem::virtualMassCoeff() const
{
    return virtualMass_->K();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseAbsorbingSystem::heatTransferCoeff() const
{
    return heatTransfer_->K();
}


Foam::tmp<Foam::volVectorField> Foam::twoPhaseAbsorbingSystem::liftForce() const
{
    return lift_->F<vector>();
}


Foam::tmp<Foam::volVectorField>
Foam::twoPhaseAbsorbingSystem::wallLubricationForce() const
{
    return wallLubrication_->F<vector>();
}


Foam::tmp<Foam::volVectorField>
Foam::twoPhaseAbsorbingSystem::turbulentDispersionForce() const
{
    return turbulentDispersion_->F<vector>();
}


void Foam::twoPhaseAbsorbingSystem::solve()
{
    const Time& runTime = mesh_.time();

    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    const surfaceScalarField& phi1 = phase1_.phi();
    const surfaceScalarField& phi2 = phase2_.phi();

    const dictionary& alphaControls = mesh_.solverDict
    (
        alpha1.name()
    );

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));
    Switch implicitPhasePressure
    (
        alphaControls.lookupOrDefault<Switch>("implicitPhasePressure", false)
    );

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    alpha1.correctBoundaryConditions();


    surfaceScalarField phic("phic", phi_);
    surfaceScalarField phir("phir", phi1 - phi2);

    surfaceScalarField alpha1f(fvc::interpolate(max(alpha1, scalar(0))));

    tmp<surfaceScalarField> pPrimeByA;

    if (implicitPhasePressure)
    {
        const volScalarField& rAU1 = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName("rAU", phase1_.name())
        );
        const volScalarField& rAU2 = mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName("rAU", phase2_.name())
        );

        pPrimeByA =
            fvc::interpolate(rAU1*phase1_.turbulence().pPrime())
          + fvc::interpolate(rAU2*phase2_.turbulence().pPrime());

        surfaceScalarField phiP
        (
            pPrimeByA()*fvc::snGrad(alpha1, "bounded")*mesh_.magSf()
        );

        phic += alpha1f*phiP;
        phir += phiP;
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", dgdt_.dimensions(), 0.0)
        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi_)*min(alpha1, scalar(1))
        );

        forAll(dgdt_, celli)
        {
            if (dgdt_[celli] > 0.0)
            {
                Sp[celli] -= dgdt_[celli]/max(1.0 - alpha1[celli], 1e-4);
                Su[celli] += dgdt_[celli]/max(1.0 - alpha1[celli], 1e-4);
            }
            else if (dgdt_[celli] < 0.0)
            {
                Sp[celli] += dgdt_[celli]/max(alpha1[celli], 1e-4);
            }
        }

        surfaceScalarField alphaPhic1
        (
            fvc::flux
            (
                phic,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            )
        );

        // Ensure that the flux at inflow BCs is preserved
        forAll(alphaPhic1.boundaryField(), patchi)
        {
            fvsPatchScalarField& alphaPhic1p =
                alphaPhic1.boundaryField()[patchi];

            if (!alphaPhic1p.coupled())
            {
                const scalarField& phi1p = phi1.boundaryField()[patchi];
                const scalarField& alpha1p = alpha1.boundaryField()[patchi];

                forAll(alphaPhic1p, facei)
                {
                    if (phi1p[facei] < 0)
                    {
                        alphaPhic1p[facei] = alpha1p[facei]*phi1p[facei];
                    }
                }
            }
        }

        if (nAlphaSubCycles > 1)
        {
            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhic10(alphaPhic1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi_,
                    alphaPhic10,
                    (alphaSubCycle.index()*Sp)(),
                    (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                    phase1_.alphaMax(),
                    0
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase1_.alphaPhi() = alphaPhic10;
                }
                else
                {
                    phase1_.alphaPhi() += alphaPhic10;
                }
            }
            
            phase1_.alphaPhi() /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhic1,
                Sp,
                Su,
                phase1_.alphaMax(),
                0
            );

            phase1_.alphaPhi() = alphaPhic1;
        }

        if (implicitPhasePressure)
        {
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian(alpha1f*pPrimeByA(), alpha1, "bounded")
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            phase1_.alphaPhi() += alpha1Eqn.flux();
        }

        phase1_.alphaRhoPhi() =
            fvc::interpolate(phase1_.rho())*phase1_.alphaPhi();

        phase2_.alphaPhi() = phi_ - phase1_.alphaPhi();
        alpha2 = scalar(1) - alpha1;
        phase2_.alphaRhoPhi() =
            fvc::interpolate(phase2_.rho())*phase2_.alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
            
        correct();
            
    }
    
}



void Foam::twoPhaseAbsorbingSystem::correct()
{
    const volScalarField& alpha1 = phase1_;
    const volScalarField& alpha2 = phase2_;
    volScalarField& wCH4_1 = phase1_.wCH4();
    volScalarField& wCH4_2 = phase2_.wCH4();
    const dimensionedScalar& DCH4_1 = phase1_.DCH4();
    const dimensionedScalar& DCH4_2 = phase2_.DCH4();
    
    volScalarField wCH4_1old(wCH4_1);
    
    // Gas phase
    fvScalarMatrix wCH4_1Eqn
    (
        fvm::ddt(alpha1, phase1_.rho(), wCH4_1)
      + fvm::div(phase1_.alphaRhoPhi(), wCH4_1)
      ==
        fvm::laplacian(alpha1*phase1_.rho()*DCH4_1, wCH4_1)
      - absorption_
    );
    
    wCH4_1Eqn.relax();
    wCH4_1Eqn.solve();
    wCH4_1 = max(0.0, min(1.0, wCH4_1));
    
    phase1_.correctRho(wCH4_1old, Mair_);
    
    Info<< wCH4_1.name() 
        << "  Min(wCH4) = " << min(wCH4_1).value()
        << "  Max(wCH4) = " << max(wCH4_1).value()
        << endl;
        
    // Liquid phase
    fvScalarMatrix wCH4_2Eqn
    (
        fvm::ddt(alpha2, phase2_.rho(), wCH4_2)
      + fvm::div(phase2_.alphaRhoPhi(), wCH4_2)
      ==
        fvm::laplacian(alpha2*phase2_.rho()*DCH4_2, wCH4_2)
      + absorption_
    );
    
    wCH4_2Eqn.relax();
    wCH4_2Eqn.solve();
    wCH4_2 = max(0.0, min(wCH4sat_, wCH4_2));
    
    Info<< wCH4_2.name() 
        << "  Min(wCH4) = " << min(wCH4_2).value()
        << "  Max(wCH4) = " << max(wCH4_2).value()
        << endl;
    
    updateAbsorption();
    phase1_.correct();
    phase2_.correct();
}


void Foam::twoPhaseAbsorbingSystem::correctTurbulence()
{
    phase1_.turbulence().correct();
    phase2_.turbulence().correct();
}


bool Foam::twoPhaseAbsorbingSystem::read()
{
    if (regIOobject::read())
    {
        bool readOK = true;

        readOK &= phase1_.read(*this);
        readOK &= phase2_.read(*this);

        // models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


const Foam::dragModel&
Foam::twoPhaseAbsorbingSystem::drag(const phaseModel& phase) const
{
    return drag_->phaseModel(phase);
}


const Foam::virtualMassModel&
Foam::twoPhaseAbsorbingSystem::virtualMass(const phaseModel& phase) const
{
    return virtualMass_->phaseModel(phase);
}


const Foam::dimensionedScalar& Foam::twoPhaseAbsorbingSystem::sigma() const
{
    return pair_->sigma();
}


// ************************************************************************* //
