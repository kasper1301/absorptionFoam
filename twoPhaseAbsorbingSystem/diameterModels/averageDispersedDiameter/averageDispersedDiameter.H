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

Class
    Foam::diameterModels::ADD

Description
    ADD (Averaged Dispersed Diameter) bubble diameter model.

    Solves a transport equation for the bubble diameter based on an equilbrium
    diameter. 

    \verbatim
        "A CFD Analysis  of the Air Entrainment Rate Due to a Plunging Steel Jet
        Combining Mathematical Models for Dispersed and Separated Multiphase 
        Flows"
        H. Laux,
        S. T. Johansen,
        Fluid Flow Phenomena in Metals Processing: Proceedings of a Symposium
    \endverbatim

SourceFiles
    IATE.C

\*---------------------------------------------------------------------------*/

#ifndef ADD_H
#define ADD_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{


/*---------------------------------------------------------------------------*\
                           Class ADD Declaration
\*---------------------------------------------------------------------------*/

class ADD
:
    public diameterModel
{
    // Private data

        //- The average dispersed diameter of the phase (Sauter-mean)
        volScalarField d_;

        //- Maximum diameter
        dimensionedScalar dMax_;
        
        //- Minimum diameter
        dimensionedScalar dMin_;

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;
        
        //- Constant n in equilibrium diameter
        dimensionedScalar n_;

        //- Constant m in equilibrium diameter
        dimensionedScalar m_;
        
        //- Constant C1 in equilibrium diameter
        dimensionedScalar C1_;
        
        //- Constant C2 in equilibrium diameter
        dimensionedScalar C2_;

        //- Constant Cb in break-up relaxation time
        dimensionedScalar Cb_;
        
        //- Constant Cc in coalescence relaxation time
        dimensionedScalar Cc_;
        
        //- Constant Cmu in turbulence model
        dimensionedScalar Cmu_;
        
        //- Maximum dispersed phase fraction
        dimensionedScalar alphaMax_;
        
        //- Effective bubble diffusivity
        volScalarField Deff_;
        
        //- Equilibrium diameter
        volScalarField deq_;
        
        //- Relaxation time
        volScalarField tauRel_;

    // Private member functions
        //- Effective bubble diffusivity
        void Deff();
        
        //- Equilibrium diameter
        void deq();
        
        //- Relaxation time
        void tauRel();

        tmp<volScalarField> tauB() const;
        tmp<volScalarField> tauC() const;
        tmp<volScalarField> tauK() const;
public:

    //- Runtime type information
    TypeName("ADD");


    // Constructors

        //- Construct from components
        ADD
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~ADD();


    // Member Functions

        //- Return the interfacial curvature
        //~ const volScalarField& kappai() const
        //~ {
            //~ return kappai_;
        //~ }

        //- Return the interfacial area
        //~ tmp<volScalarField> a() const
        //~ {
            //~ return phase_*kappai_;
        //~ }

        //- Return the diameter
        virtual tmp<volScalarField> d() const
        {
            return d_;
        }

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
