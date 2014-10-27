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
    Deff_(Deff()),
    deq_(deq())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::ADD::~ADD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::ADD::Deff() const
{
    //- Phase fraction
    const volScalarField& alpha = phase_;
    //- Density of dispersed phase
    const volScalarField& rhod = phase_.rho();
    //- Density of continous phase
    const volScalarField& rhoc = phase_.otherPhase().rho();
    //- Turbulent viscosity
    const volScalarField& mut = phase_.otherPhase().turbulence().mut();
    
    return alpha*rhod*mut/rhoc;
}

Foam::tmp<Foam::volScalarField> Foam::diameterModels::ADD::deq() const
{
    //- Phase fraction
    const volScalarField& alpha = phase_;
    //- Density of continous phase
    const volScalarField& rhoc = phase_.otherPhase().rho();
    //- Viscosity of continous phase
    const volScalarField& muc = phase_.otherPhase().turbulence().mu();
    //- Viscosity of dispersed phase
    const volScalarField& mud = phase_.turbulence().mu();
    //- Turbulent dispersion
    const volScalarField& epsilon = phase_.otherPhase().turbulence().epsilon();
    
    return C1_*pow(alpha,n_)*pow(sigma_/rhoc,0.6)*pow(mud/muc,m_)*epsilon + C2_;
}



void Foam::diameterModels::ADD::correct()
{
    Deff_ = Deff();
    deq_ = deq();
    
}


bool Foam::diameterModels::ADD::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    //~ diameterProperties_.lookup("dMax") >> dMax_;
    //~ diameterProperties_.lookup("dMin") >> dMin_;


    return true;
}


// ************************************************************************* //
