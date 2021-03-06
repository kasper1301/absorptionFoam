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
    Foam::twoPhaseAbsorbingSystem

Description

SourceFiles
    twoPhaseAbsorbingSystem.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseAbsorbingSystem_H
#define twoPhaseAbsorbingSystem_H

#include "IOdictionary.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dragModel;
class virtualMassModel;
class heatTransferModel;
class liftModel;
class wallLubricationModel;
class turbulentDispersionModel;

class blendingMethod;
template <class modelType> class BlendedInterfacialModel;

/*---------------------------------------------------------------------------*\
                      Class twoPhaseAbsorbingSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseAbsorbingSystem
:
    public IOdictionary
{
private:

    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Phase model 1
        phaseModel phase1_;

        //- Phase model 2
        phaseModel phase2_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //-  Dilatation term
        volScalarField dgdt_;
        
        //- Absorption rate
        volScalarField absorption_;
        
        //- Saturation mass fraction
        volScalarField wCH4sat_;
        
        //- Salting coefficient
        dimensionedScalar ks_;
        
        //- Ionic strength of seawater
        dimensionedScalar I_;
        
        //- Henry's law constant
        dimensionedScalar kHe_;
        
        //- Molecular weight of methane
        dimensionedScalar MCH4_;
        
        //- Molecular weight of air
        dimensionedScalar Mair_;
        
        //- Molecular weight of seawater
        dimensionedScalar Msw_;
        
        //- Unordered phase pair
        autoPtr<phasePair> pair_;

        //- Phase pair for phase 1 dispersed in phase 2
        autoPtr<orderedPhasePair> pair1In2_;

        //- Phase pair for phase 2 dispersed in phase 1
        autoPtr<orderedPhasePair> pair2In1_;

        //- Blending methods
        HashTable<autoPtr<blendingMethod>, word, word::hash> blendingMethods_;

        //- Drag model
        autoPtr<BlendedInterfacialModel<dragModel> > drag_;

        //- Virtual mass model
        autoPtr<BlendedInterfacialModel<virtualMassModel> > virtualMass_;

        //- Heat transfer model
        autoPtr<BlendedInterfacialModel<heatTransferModel> > heatTransfer_;

        //- Lift model
        autoPtr<BlendedInterfacialModel<liftModel> > lift_;

        //- Wall lubrication model
        autoPtr<BlendedInterfacialModel<wallLubricationModel> >
            wallLubrication_;

        //- Wall lubrication model
        autoPtr<BlendedInterfacialModel<turbulentDispersionModel> >
            turbulentDispersion_;


    // Private member functions

        //- Return the mixture flux
        tmp<surfaceScalarField> calcPhi() const;


public:

    // Constructors

        //- Construct from fvMesh
        twoPhaseAbsorbingSystem(const fvMesh&, const dimensionedVector& g);


    //- Destructor
    virtual ~twoPhaseAbsorbingSystem();


    // Member Functions
    
        //- Update the absorption
        void updateAbsorption();

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Return the drag coefficient
        tmp<volScalarField> dragCoeff() const;

        //- Return the virtual mass coefficient
        tmp<volScalarField> virtualMassCoeff() const;

        //- Return the heat transfer coefficient
        tmp<volScalarField> heatTransferCoeff() const;

        //- Return the lift force
        tmp<volVectorField> liftForce() const;

        //- Return the wall lubrication force
        tmp<volVectorField> wallLubricationForce() const;

        //- Return the wall lubrication force
        tmp<volVectorField> turbulentDispersionForce() const;

        //- Solve for the two-phase-fractions
        void solve();

        //- Correct two-phase properties other than turbulence
        void correct();

        //- Correct two-phase turbulence
        void correctTurbulence();

        //- Read base phaseProperties dictionary
        bool read();

        // Access

            //- Return the drag model for the supplied phase
            const dragModel& drag(const phaseModel& phase) const;

            //- Return the virtual mass model for the supplied phase
            const virtualMassModel& virtualMass(const phaseModel& phase) const;

            //- Return the surface tension coefficient
            const dimensionedScalar& sigma() const;

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return phase model 1
            inline const phaseModel& phase1() const;

            //- Return non-const access to phase model 1
            inline phaseModel& phase1();

            //- Return phase model 2
            inline const phaseModel& phase2() const;

            //- Return non-const access to phase model 2
            inline phaseModel& phase2();

            //- Return the phase not given as an argument
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Return the mixture flux
            inline const surfaceScalarField& phi() const;

            //- Return non-const access to the the mixture flux
            inline surfaceScalarField& phi();

            //- Return the dilatation term
            inline const volScalarField& dgdt() const;

            //- Return non-const access to the dilatation parameter
            inline volScalarField& dgdt();
            
            //- Return the absorption
            inline const volScalarField& absorption() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseAbsorbingSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
