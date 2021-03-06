/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "phaseModel.H"
#include "twoPhaseAbsorbingSystem.H"
#include "diameterModel.H"
#include "fvMatrix.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "dragModel.H"
#include "heatTransferModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const twoPhaseAbsorbingSystem& fluid,
    const dictionary& phaseProperties,
    const word& phaseName
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("alpha", dimless, 0)
    ),
    fluid_(fluid),
    name_(phaseName),
    phaseDict_
    (
        phaseProperties.subDict(name_)
    ),
    alphaMax_(phaseDict_.lookupOrDefault("alphaMax", 1.0)),
    thermo_(rhoThermo::New(fluid.mesh(), name_)),
    U_
    (
        IOobject
        (
            IOobject::groupName("U", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    psi_
    (
        thermo_().psi()
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    ),
    alphaRhoPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaRhoPhi", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0)
    ),
    wCH4_
    (
        IOobject
        (
            IOobject::groupName("wCH4", name_),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh()
    ),
    DCH4_
    (
        "DCH4",
        dimensionSet(0, 2, -1, 0, 0),
        phaseDict_.lookup("DCH4")
    ),
    MR_
    (
        "MR",
        dimless,
        phaseDict_.lookup("MR")
    )
{
    thermo_->validate("phaseModel " + name_, "h", "e");

    const word phiName = IOobject::groupName("phi", name_);

    IOobject phiHeader
    (
        phiName,
        fluid_.mesh().time().timeName(),
        fluid_.mesh(),
        IOobject::NO_READ
    );

    if (phiHeader.headerOk())
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluid_.mesh()
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U_) & fluid_.mesh().Sf(),
                phiTypes
            )
        );
    }

    dPtr_ = diameterModel::New
    (
        phaseDict_,
        *this
    );

    turbulence_ =
        PhaseCompressibleTurbulenceModel<phaseModel>::New
        (
            *this,
            thermo_->rho(),
            U_,
            alphaRhoPhi_,
            phi(),
            *this
        );
        
    thermo_().rho() /= (1.0 + (MR_ - 1.0)*wCH4_);
    psi_ /= (1.0 + (MR_ - 1.0)*wCH4_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseModel& Foam::phaseModel::otherPhase() const
{
    return fluid_.otherPhase(*this);
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::d() const
{
    return dPtr_().d();
}

Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::phaseModel::turbulence()
{
    return turbulence_();
}

const Foam::PhaseCompressibleTurbulenceModel<Foam::phaseModel>&
Foam::phaseModel::turbulence() const
{
    return turbulence_();
}

void Foam::phaseModel::correct()
{
    return dPtr_->correct();
}

void Foam::phaseModel::correctRho
(
    const volScalarField& wCH4_old,
    const dimensionedScalar& Mw
)
{
    dimensionedScalar R
    (
        "R",
        dimensionSet(0, 2, -2, -1, 0),
        8.3145
    );
    
    volScalarField dpsidw = 
    (
        Mw*(1.0 - MR_)/(sqr(1.0 + (MR_ - 1.0)*wCH4_old))
      / (R*thermo_().T())
    );
    
    psi_ += dpsidw*(wCH4_ - wCH4_old);
    
    thermo_().rho() += dpsidw*thermo_().p()*(wCH4_ - wCH4_old);
}

bool Foam::phaseModel::read(const dictionary& phaseProperties)
{
    phaseDict_ = phaseProperties.subDict(name_);
    return dPtr_->read(phaseDict_);
}


// ************************************************************************* //
