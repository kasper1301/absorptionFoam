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

#include "averageDispersedDiameter.H"
#include "twoPhaseAbsorbingSystem.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcAverage.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(ADD, 3);

    addToRunTimeSelectionTable
    (
        diameterModel,
        ADD,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::ADD::ADD
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase_.U().time().timeName(),
            phase_.U().mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase_.U().mesh()
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        diameterProperties_.lookup("residualAlpha")
    ),
    sigma_(phase.fluid().sigma()),
    n_
    (
        "n",
        dimless,
        diameterProperties_.lookup("n")
    ),
    m_
    (
        "m",
        dimless,
        diameterProperties_.lookup("m")
    ),
    C1_
    (
        "C1",
        dimless,
        diameterProperties_.lookup("C1")
    ),
    C2_
    (
        "C2",
        dimLength,
        diameterProperties_.lookup("C2")
    ),
    Cb_
    (
        "Cb",
        dimless,
        diameterProperties_.lookup("Cb")
    ),
    Cc_
    (
        "Cc",
        dimless,
        diameterProperties_.lookup("Cc")
    ),
    Cmu_
    (
        "Cmu",
        dimless,
        diameterProperties_.lookup("Cmu")
    ),
    alphaMax_
    (
        "Cmu",
        dimless,
        diameterProperties_.lookup("alphaMax")
    ),
    Deff_
    (
        IOobject
        (
            "Deff",
            phase_.U().time().timeName(),
            phase_.U().mesh()
        ),
        phase_.U().mesh(),
        dimensionedScalar("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0)
    ),
    deq_
    (
        IOobject
        (
            "deq",
            phase_.U().time().timeName(),
            phase_.U().mesh()
        ),
        phase_.U().mesh(),
        dimensionedScalar("zero", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0)
    ),
    tauRel_
    (
        IOobject
        (
            "tauRel",
            phase_.U().time().timeName(),
            phase_.U().mesh()
        ),
        phase_.U().mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    )
{}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::ADD::~ADD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::ADD::Deff()
{
    //- Phase fraction
    const volScalarField& alpha = phase_;
    //- Density of dispersed phase
    const volScalarField& rhod = phase_.rho();
    //- Density of continous phase
    const volScalarField& rhoc = phase_.otherPhase().rho();
    //- Turbulent viscosity
    const volScalarField& mut = phase_.otherPhase().turbulence().mut();
    
    Deff_ = alpha*rhod*mut/rhoc;
}

void Foam::diameterModels::ADD::deq()
{
    //- Phase fraction
    const volScalarField& alpha = phase_;
    //- Density of continous phase
    const volScalarField& rhoc = phase_.otherPhase().rho();
    //- Viscosity of continous phase
    const volScalarField& muc = phase_.otherPhase().mu();
    //- Viscosity of dispersed phase
    const volScalarField& mud = phase_.mu();
    //- Turbulent dispersion
    const volScalarField& epsilon = phase_.otherPhase().turbulence().epsilon();

    deq_ = C1_*pow(alpha,n_)*pow(sigma_/rhoc,0.6)*pow(mud/muc,m_)*epsilon + C2_;
}

void Foam::diameterModels::ADD::tauRel()
{
    const volScalarField& nu = phase_.otherPhase().turbulence().nu();
    const volScalarField& epsilon = phase_.otherPhase().turbulence().epsilon();
    const volScalarField tauK = 6.0*sqrt(nu/epsilon);
    const volScalarField tauBreakUp = tauB();
    const volScalarField tauCoalescence = tauC();
    
    forAll(d_, celli)
    {
        if (d_[celli] > deq_[celli])
        {
            tauRel_[celli] = tauBreakUp[celli];
        }
        else if (d_[celli] < deq_[celli])
        {
            tauRel_[celli] = tauCoalescence[celli];
        }
    }
    
    tauRel_ = max(tauRel_, tauK);
    
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::ADD::tauB() const
{
    //- Turbulent dispersion
    const volScalarField& epsilon = phase_.otherPhase().turbulence().epsilon();
    return Cb_*pow(d_,2.0/3.0)/cbrt(epsilon);
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::ADD::tauC() const
{
    //- Phase fraction
    const volScalarField& alpha = phase_;
    //- Density of dispersed phase
    const volScalarField& rho = phase_.rho();
    //- Viscosity of continous phase
    const volScalarField& mu = phase_.otherPhase().mu();
    //- Turbulent kinetic energy
    const volScalarField& k = phase_.otherPhase().turbulence().k();
    //- Turbulent dispersion
    const volScalarField& epsilon = phase_.otherPhase().turbulence().epsilon();
    
    const volScalarField tauL = 1.5*Cmu_*k/epsilon;
    const volScalarField tauP = rho*sqr(d_)/(18*mu);
    const volScalarField St = tauP/tauL;
    const volScalarField kd = k/(1 + St);
    const volScalarField uRMS = sqrt(2.0/3.0*kd);
    const volScalarField s = d_*cbrt
    (
        constant::mathematical::pi/6.0*(alphaMax_ - alpha)/alpha
    );
    return Cc_*s/uRMS;
    
    return tauL;
}

void Foam::diameterModels::ADD::correct()
{
    Deff();
    deq();
    tauRel();
    
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField rhoPhi = fvc::interpolate(rho)*phase_.phi();
    
    volScalarField& d = d_;

    volScalarField& Deff = Deff_;
    volScalarField& deq = deq_;
    volScalarField& tauR = tauRel_;
    
    fvScalarMatrix dEqn
    (
        fvm::ddt(rho, d)
      + fvm::div(rhoPhi, d) 
        ==
        fvm::laplacian(Deff, d)
      - fvm::Sp(rho/tauR, d)
      + rho/tauR*deq
    );
    
    dEqn.relax();
    dEqn.solve();
    
}


bool Foam::diameterModels::ADD::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    //~ diameterProperties_.lookup("dMax") >> dMax_;
    //~ diameterProperties_.lookup("dMin") >> dMin_;


    return true;
}


// ************************************************************************* //
