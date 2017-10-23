/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::dynamicMultiFieldRefineBalancedFvMesh

Description
    A fvMesh with built-in refinement and run-time load balancing.

SourceFiles
    dynamicMultiFieldRefineBalancedFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMultiFieldRefineBalancedFvMesh_H
#define dynamicMultiFieldRefineBalancedFvMesh_H


#include "dynamicFvMesh.H"
#include "dynamicMultiFieldRefineFvMesh.H"
#include "hexRef8.H"
#include "PackedBoolList.H"
#include "Switch.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "HashTable.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "PtrDictionary.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class dynamicMultiFieldRefineBalancedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMultiFieldRefineBalancedFvMesh
:
    public dynamicMultiFieldRefineFvMesh
{

private:
        //-
        volScalarField* internalRefinementFieldPtr_;

        //-
        HashTable< Pair<scalar> > gradFields_;

        //-
        HashTable< Pair<scalar> > curlFields_;

        //-
        PtrList<entry> refinedRegions_;

        //-
        Switch enableRefinementControl_;

        //-
        void updateRefinementField();

        //-
        void readRefinementDict();

        //-
        Pair<scalar> readRefinementPoints();

        //- Disallow default bitwise copy construct
        dynamicMultiFieldRefineBalancedFvMesh
        (const dynamicMultiFieldRefineBalancedFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicMultiFieldRefineBalancedFvMesh&);

        label topParentID(label p);
public:

    //- Runtime type information
    TypeName("dynamicMultiFieldRefineBalancedFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicMultiFieldRefineBalancedFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicMultiFieldRefineBalancedFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Template to update all volField boundaries
        template<class Type> void correctBoundaries();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicMultiFieldRefineBalancedFvMeshTemplates.C"
#endif

#endif

// ************************************************************************* //
