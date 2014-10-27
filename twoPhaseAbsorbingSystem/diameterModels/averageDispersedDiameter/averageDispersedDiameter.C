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
    
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::ADD::~ADD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Placeholder for the nucleation/condensation model
// Foam::tmp<Foam::volScalarField> Foam::diameterModels::IATE::Rph() const
// {
//     const volScalarField& T = phase_.thermo().T();
//     const volScalarField& p = phase_.thermo().p();
//
//     scalar A, B, C, sigma, vm, Rph;
//
//     volScalarField ps(1e5*pow(10, A - B/(T + C)));
//     volScalarField Dbc
//     (
//         4*sigma*vm/(constant::physicoChemical::k*T*log(p/ps))
//     );
//
//     return constant::mathematical::pi*sqr(Dbc)*Rph;
// }

void Foam::diameterModels::ADD::correct()
{
    // Initialise the accumulated source term to the dilatation effect
    //~ volScalarField R
    //~ (
        //~ (
            //~ (1.0/3.0)
           //~ /max
            //~ (
                //~ fvc::average(phase_ + phase_.oldTime()),
                //~ residualAlpha_
            //~ )
        //~ )
       //~ *(fvc::ddt(phase_) + fvc::div(phase_.alphaPhi()))
    //~ );
//~ 
    //~ // Accumulate the run-time selectable sources
    //~ forAll(sources_, j)
    //~ {
        //~ R -= sources_[j].R();
    //~ }
//~ 
    //~ // Construct the interfacial curvature equation
    //~ fvScalarMatrix kappaiEqn
    //~ (
        //~ fvm::ddt(kappai_) + fvm::div(phase_.phi(), kappai_)
      //~ - fvm::Sp(fvc::div(phase_.phi()), kappai_)
     //~ ==
      //~ - fvm::SuSp(R, kappai_)
    //~ //+ Rph() // Omit the nucleation/condensation term
    //~ );
//~ 
    //~ kappaiEqn.relax();
    //~ kappaiEqn.solve();
//~ 
    //~ // Update the Sauter-mean diameter
    //~ d_ = dsm();
}


bool Foam::diameterModels::ADD::read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    //~ diameterProperties_.lookup("dMax") >> dMax_;
    //~ diameterProperties_.lookup("dMin") >> dMin_;


    return true;
}


// ************************************************************************* //
