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


\*---------------------------------------------------------------------------*/

#include "dynamicMultiFieldRefineBalancedFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        dynamicMultiFieldRefineBalancedFvMesh,
        0
    );

    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        dynamicMultiFieldRefineBalancedFvMesh,
        IOobject
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label
Foam::dynamicMultiFieldRefineBalancedFvMesh::topParentID(label p)
{
    label nextP = meshCutter().history().splitCells()[p].parent_;
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParentID(nextP);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMultiFieldRefineBalancedFvMesh::dynamicMultiFieldRefineBalancedFvMesh
(
    const IOobject& io
)
:
    dynamicMultiFieldRefineFvMesh(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMultiFieldRefineBalancedFvMesh::~dynamicMultiFieldRefineBalancedFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicMultiFieldRefineBalancedFvMesh::update()
{
   //Part 1 - Call normal update from dynamicMultiFieldRefineFvMesh
    bool hasChanged =
        dynamicMultiFieldRefineFvMesh::update();

    // Part 2 - Load Balancing
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicMultiFieldRefineFvMeshCoeffs")
    );

    Switch enableBalancing =
        refineDict.lookup("enableBalancing");

    if ( Pstream::parRun() && hasChanged )
    {
        const scalar allowableImbalance =
            readScalar(refineDict.lookup("allowableImbalance"));

        //First determine current level of imbalance - do this for all
        // parallel runs with a changing mesh, even if balancing is disabled
        label nGlobalCells = globalData().nTotalCells();
        scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
        scalar localImbalance = mag(scalar(nCells()) - idealNCells);
        Foam::reduce(localImbalance, maxOp<scalar>());
        scalar maxImbalance = localImbalance/idealNCells;

        Info<<"Maximum imbalance = " << 100*maxImbalance << " %" << endl;

        //If imbalanced, construct weighted coarse graph (level 0) with node
        // weights equal to their number of subcells. This partitioning works
        // as long as the number of level 0 cells is several times greater than
        // the number of processors.
        if( maxImbalance > allowableImbalance && enableBalancing)
        {
            Info<< "Re-balancing dynamically refined mesh" << endl;

            const labelIOList& cellLevel = meshCutter().cellLevel();

            Map<label> coarseIDmap(100);  //map maps label to labels. So entry should be like, 1 label to 1 label
                                          //construct with initial size 100.

            labelList uniqueIndex(nCells(),0); //unique index, has space for 1 entry per cell.

            label nCoarse = 0;

            //Loop over all cells:
            forAll(cells(), cellI)
            {
                //If cell was refined:
                if( cellLevel[cellI] > 0 )
                {
                    //get unique index for this cell, add to number of cells its top parrent cell ID:
                    uniqueIndex[cellI] = nCells() + topParentID
                    (
                        meshCutter().history().parentIndex(cellI)
                    );
                }
                else
                {
                    //If the cell was not refined, give it its index:
                    uniqueIndex[cellI] = cellI;
                }

                //Insert coarse cell ID into the map:

                if( coarseIDmap.insert(uniqueIndex[cellI], nCoarse) )
                {
                    ++nCoarse;
                }
            }

            // Convert to local sequential indexing and calculate coarse
            // points and weights
            labelList localIndex(nCells(),0);
            pointField coarsePoints(nCoarse,vector::zero);
            scalarField coarseWeights(nCoarse,0.0);

            forAll(uniqueIndex, cellI)
            {
                //Map from current cell to its coarse cell:
                localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];

                // If 2D refinement (quadtree) is ever implemented, this '3'
                // should be set in general as the number of refinement
                // dimensions.
                label w = (1 << (3*cellLevel[cellI]));

                //Add one weight per cell in its coarse cell:
                coarseWeights[localIndex[cellI]] += 1.0;
                coarsePoints[localIndex[cellI]] += C()[cellI]/w;
            }

            //Set up decomposer - a separate dictionary is used here so
            // you can use a simple partitioning for decomposePar and
            // ptscotch for the rebalancing (or any chosen algorithms)
            autoPtr<decompositionMethod> decomposer
            (
                decompositionMethod::New
                (
                    IOdictionary
                    (
                        IOobject
                        (
                            "balanceParDict",
                            time().system(),
                            *this,
                            IOobject::MUST_READ_IF_MODIFIED,
                            IOobject::NO_WRITE
                        )
                    )
                )
            );

            //Decompose, get the list to witch processors should each cell go
            labelList finalDecomp = decomposer().decompose
            (
                *this,
                localIndex,
                coarsePoints,
                coarseWeights
            );

            scalar tolDim = globalMeshData::matchTol_ * bounds().mag();

            //Construct distributer:
            fvMeshDistribute distributor(*this, tolDim);

            //Distribute mesh, after this operation mesh on proc. changes
            //according to decompositionMethod indexing:
            autoPtr<mapDistributePolyMesh> map =
                  distributor.distribute(finalDecomp);

            meshCutter_.distribute(map);

            //Correct values on all cyclic patches
            correctBoundaries<scalar>();
            correctBoundaries<vector>();
            correctBoundaries<sphericalTensor>();
            correctBoundaries<symmTensor>();
            correctBoundaries<tensor>();
        }
    }

    return hasChanged;
}


// ************************************************************************* //
