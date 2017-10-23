/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "dynamicMultiFieldRefineFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
//#include "pointFields.H"
#include "SortableList.H"
#include "string.H"
#include "sigFpe.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMultiFieldRefineFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicMultiFieldRefineFvMesh, IOobject);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::dynamicMultiFieldRefineFvMesh::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val)
        {
            n++;
        }

        // debug also serves to get-around Clang compiler trying to optimsie
        // out this forAll loop under O3 optimisation
        if (debug)
        {
            Info<< "n=" << n << endl;
        }
    }

    return n;
}


void Foam::dynamicMultiFieldRefineFvMesh::calculateProtectedCells
(
    PackedBoolList& unrefineableCell
) const
{
    if (protectedCell_.empty())
    {
        unrefineableCell.clear();
        return;
    }

    const labelList& cellLevel = meshCutter_.cellLevel();

    unrefineableCell = protectedCell_;

    // Get neighbouring cell level
    labelList neiLevel(nFaces()-nInternalFaces());

    for (label facei = nInternalFaces(); facei < nFaces(); facei++)
    {
        neiLevel[facei-nInternalFaces()] = cellLevel[faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(*this, neiLevel);


    while (true)
    {
        // Pick up faces on border of protected cells
        boolList seedFace(nFaces(), false);

        forAll(faceNeighbour(), facei)
        {
            label own = faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            label nei = faceNeighbour()[facei];
            bool neiProtected = unrefineableCell.get(nei);

            if (ownProtected && (cellLevel[nei] > cellLevel[own]))
            {
                seedFace[facei] = true;
            }
            else if (neiProtected && (cellLevel[own] > cellLevel[nei]))
            {
                seedFace[facei] = true;
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            label own = faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            if
            (
                ownProtected
             && (neiLevel[facei-nInternalFaces()] > cellLevel[own])
            )
            {
                seedFace[facei] = true;
            }
        }

        syncTools::syncFaceList(*this, seedFace, orEqOp<bool>());


        // Extend unrefineableCell
        bool hasExtended = false;

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = faceOwner()[facei];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }

                label nei = faceNeighbour()[facei];
                if (unrefineableCell.get(nei) == 0)
                {
                    unrefineableCell.set(nei, 1);
                    hasExtended = true;
                }
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = faceOwner()[facei];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }
            }
        }

        if (!returnReduce(hasExtended, orOp<bool>()))
        {
            break;
        }
    }
}


void Foam::dynamicMultiFieldRefineFvMesh::readDict()
{
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
        ).optionalSubDict(typeName + "Coeffs")
    );

    List<Pair<word>> fluxVelocities = List<Pair<word>>
    (
        refineDict.lookup("correctFluxes")
    );
    // Rework into hashtable.
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }

    dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));
}


// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicMultiFieldRefineFvMesh::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            label oldFacei = map().faceMap()[facei];

            if (oldFacei >= nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << facei
                    << " fc:" << faceCentres()[facei]
                    << " originates from boundary oldFace:" << oldFacei
                    << abort(FatalError);
            }
        }
    }

    //    // Remove the stored tet base points
    //    tetBasePtIsPtr_.clear();
    //    // Remove the cell tree
    //    cellTreePtr_.clear();

    // Update fields
    updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        labelHashSet masterFaces(4*cellsToRefine.size());

        forAll(faceMap, facei)
        {
            label oldFacei = faceMap[facei];

            if (oldFacei >= 0)
            {
                label masterFacei = reverseFaceMap[oldFacei];

                if (masterFacei < 0)
                {
                    FatalErrorInFunction
                        << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << facei << abort(FatalError);
                }
                else if (masterFacei != facei)
                {
                    masterFaces.insert(masterFacei);
                }
            }
        }
        if (debug)
        {
            Pout<< "Found " << masterFaces.size() << " split faces " << endl;
        }

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning."
                    << endl;
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (UName == "NaN")
            {
                Pout<< "Setting surfaceScalarField " << iter.key()
                    << " to NaN" << endl;

                surfaceScalarField& phi = *iter();

                sigFpe::fillNan(phi.primitiveFieldRef());

                continue;
            }

            if (debug)
            {
                Pout<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );

            // Recalculate new internal faces.
            for (label facei = 0; facei < nInternalFaces(); facei++)
            {
                label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inflated/appended
                    phi[facei] = phiU[facei];
                }
                else if (reverseFaceMap[oldFacei] != facei)
                {
                    // face-from-masterface
                    phi[facei] = phiU[facei];
                }
            }

            // Recalculate new boundary faces.
            surfaceScalarField::Boundary& phiBf =
                phi.boundaryFieldRef();
            forAll(phiBf, patchi)
            {
                fvsPatchScalarField& patchPhi = phiBf[patchi];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchi];

                label facei = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    label oldFacei = faceMap[facei];

                    if (oldFacei == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFacei] != facei)
                    {
                        // face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    facei++;
                }
            }

            // Update master faces
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label facei = iter.key();

                if (isInternalFace(facei))
                {
                    phi[facei] = phiU[facei];
                }
                else
                {
                    label patchi = boundaryMesh().whichPatch(facei);
                    label i = facei - boundaryMesh()[patchi].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchi];

                    fvsPatchScalarField& patchPhi = phiBf[patchi];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }



    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            newProtectedCell.set(celli, protectedCell_.get(oldCelli));
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicMultiFieldRefineFvMesh::unrefine
(
    const labelList& splitPoints
)
{
    polyTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setUnrefinement(splitPoints, meshMod);


    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    Map<label> faceToSplitPoint(3*splitPoints.size());

    {
        forAll(splitPoints, i)
        {
            label pointi = splitPoints[i];

            const labelList& pEdges = pointEdges()[pointi];

            forAll(pEdges, j)
            {
                label otherPointi = edges()[pEdges[j]].otherVertex(pointi);

                const labelList& pFaces = pointFaces()[otherPointi];

                forAll(pFaces, pFacei)
                {
                    faceToSplitPoint.insert(pFaces[pFacei], otherPointi);
                }
            }
        }
    }


    // Change mesh and generate map.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified faces.
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << iter.key()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << iter.key() << " none) to suppress this warning."
                    << endl;
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (debug)
            {
                Info<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            surfaceScalarField::Boundary& phiBf =
                phi.boundaryFieldRef();

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );


            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                label oldFacei = iter.key();
                label oldPointi = iter();

                if (reversePointMap[oldPointi] < 0)
                {
                    // midpoint was removed. See if face still exists.
                    label facei = reverseFaceMap[oldFacei];

                    if (facei >= 0)
                    {
                        if (isInternalFace(facei))
                        {
                            phi[facei] = phiU[facei];
                        }
                        else
                        {
                            label patchi = boundaryMesh().whichPatch(facei);
                            label i = facei - boundaryMesh()[patchi].start();

                            const fvsPatchScalarField& patchPhiU =
                                phiU.boundaryField()[patchi];
                            fvsPatchScalarField& patchPhi = phiBf[patchi];
                            patchPhi[i] = patchPhiU[i];
                        }
                    }
                }
            }
        }
    }


    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            if (oldCelli >= 0)
            {
                newProtectedCell.set(celli, protectedCell_.get(oldCelli));
            }
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::scalarField
Foam::dynamicMultiFieldRefineFvMesh::maxPointField(const scalarField& pFld) const
{
    scalarField vFld(nCells(), -GREAT);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        forAll(pCells, i)
        {
            vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointi]);
        }
    }
    return vFld;
}


// Get min of connected cell
Foam::scalarField
Foam::dynamicMultiFieldRefineFvMesh::minCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), GREAT);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        forAll(pCells, i)
        {
            pFld[pointi] = min(pFld[pointi], vFld[pCells[i]]);
        }
    }
    return pFld;
}


// Get max of connected cellFoam::scalarField
Foam::scalarField
Foam::dynamicMultiFieldRefineFvMesh::maxCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), -GREAT);

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        forAll(pCells, i)
        {
            pFld[pointi] = max(pFld[pointi], vFld[pCells[i]]);
        }
    }
    return pFld;
}


Foam::scalarField
Foam::dynamicMultiFieldRefineFvMesh::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(nPoints());

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointi] = sum/pCells.size();
    }
    return pFld;
}


Foam::scalarField Foam::dynamicMultiFieldRefineFvMesh::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);

    forAll(fld, i)
    {
        scalar err = min(fld[i]-minLevel, maxLevel-fld[i]);

        if (err >= 0)
        {
            c[i] = err;
        }
    }
    return c;
}


void Foam::dynamicMultiFieldRefineFvMesh::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    scalarField cellError
    (
        maxPointField
        (
            error
            (
                cellToPoint(vFld),
                lowerRefineLevel,
                upperRefineLevel
            )
        )
    );

    // Mark cells that are candidates for refinement.
    forAll(cellError, celli)
    {
        if (cellError[celli] > 0)
        {
            candidateCell.set(celli, 1);
        }
    }
}


Foam::labelList Foam::dynamicMultiFieldRefineFvMesh::selectRefineCells
(
    const label maxCells,
    const labelList& maxRefinementLimit,
    const PackedBoolList& candidateCell
) const
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    PackedBoolList unrefineableCell;
    calculateProtectedCells(unrefineableCell);

    // Count current selection
    label nLocalCandidates = count(candidateCell, 1);
    label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

    // Collect all cells
    DynamicList<label> candidates(nLocalCandidates);

    if (nCandidates < nTotToRefine)
    {
        forAll(candidateCell, celli)
        {
            if
            (
                cellLevel[celli] < maxRefinementLimit[celli]
             && candidateCell.get(celli)
             && (
                    unrefineableCell.empty()
                 || !unrefineableCell.get(celli)
                )
            )
            {
                candidates.append(celli);
            }
        }
    }
    else
    {
        // Sort by error? For now just truncate.
        /*for (label level = 0; level < maxRefinementLimit[celli]; level++)
        {
            forAll(candidateCell, celli)
            {
                if
                (
                    cellLevel[celli] == level
                 && candidateCell.get(celli)
                 && (
                        unrefineableCell.empty()
                     || !unrefineableCell.get(celli)
                    )
                )
                {
                    candidates.append(celli);
                }
            }

            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }*/
    }

    // Guarantee 2:1 refinement after refinement
    labelList consistentSet
    (
        meshCutter_.consistentRefinement
        (
            candidates.shrink(),
            true               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


Foam::labelList Foam::dynamicMultiFieldRefineFvMesh::selectUnrefinePoints
(
    const scalar unrefineLevel,
    const label maxRefinement,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getSplitPoints());

    DynamicList<label> newSplitPoints(splitPoints.size());

    forAll(splitPoints, i)
    {
        label pointi = splitPoints[i];

        const labelList& pCells = pointCells()[pointi];

        bool isLower = false;
        bool hasMarked = false;

        if (pFld[pointi] < unrefineLevel)
        {
            // Check that all cells are not marked
            //const labelList& pCells = pointCells()[pointi];

            //bool hasMarked = false;
            isLower = true;

            // Check that all cells are not marked
            forAll(pCells, pCelli)
            {
                if (markedCell.get(pCells[pCelli]))
                {
                    hasMarked = true;
                    break;
                }
            }

            /*if (!hasMarked)
            {
                newSplitPoints.append(pointi);
            }*/
        }

        bool hasMarked2 = true;

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(pCells, pCelli)
        {
            if (cellLevel[pCells[pCelli]] <= maxRefinement)
            {
                hasMarked2 = false;
                break;
            }
        }

        if ((isLower && !hasMarked) || hasMarked2) newSplitPoints.append(pointi);
    }


    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement
        (
            newSplitPoints,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}


Foam::labelList Foam::dynamicMultiFieldRefineFvMesh::selectUnrefinePoints2
(
    const scalarList& unrefineLevels2,
    const labelList& maxRefinements2,
    List<PackedBoolList> markedCells2,
    List<scalarField> pFld2
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getSplitPoints());

    DynamicList<label> newSplitPoints(splitPoints.size());

    const labelList& cellLevel = meshCutter_.cellLevel();

    forAll(splitPoints, i)
    {
        label pointi = splitPoints[i];

        const labelList& pCells = pointCells()[pointi];

        label pCellsMin = cellLevel[pCells[0]];

        forAll(pCells, pCelli)
        {
            if (pCellsMin == 0) break;

            if (cellLevel[pCells[pCelli]] < pCellsMin)
            {
                pCellsMin = cellLevel[pCells[pCelli]];
            }
        }

        bool unrefine = true;

        forAll(unrefineLevels2, j)
        {
            if (pCellsMin <= maxRefinements2[j])
            {
                unrefine = false;

                if (pFld2[j][pointi] < unrefineLevels2[j])
                {
                    unrefine = true;

                    forAll(pCells, pCelli)
                    {
                        if (markedCells2[j].get(pCells[pCelli]))
                        {
                            unrefine = false;
                            break;
                        }
                    }
                }
            }

            if (!unrefine) break;
        }

        if (unrefine) newSplitPoints.append(pointi);
    }

    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement
        (
            newSplitPoints,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}


void Foam::dynamicMultiFieldRefineFvMesh::extendMarkedCells
(
    PackedBoolList& markedCell
) const
{
    // Mark faces using any marked cell
    boolList markedFace(nFaces(), false);

    forAll(markedCell, celli)
    {
        if (markedCell.get(celli))
        {
            const cell& cFaces = cells()[celli];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(*this, markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label facei = 0; facei < nInternalFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(faceOwner()[facei], 1);
            markedCell.set(faceNeighbour()[facei], 1);
        }
    }
    for (label facei = nInternalFaces(); facei < nFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(faceOwner()[facei], 1);
        }
    }
}


void Foam::dynamicMultiFieldRefineFvMesh::checkEightAnchorPoints
(
    PackedBoolList& protectedCell,
    label& nProtected
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    labelList nAnchorPoints(nCells(), 0);

    forAll(pointLevel, pointi)
    {
        const labelList& pCells = pointCells(pointi);

        forAll(pCells, pCelli)
        {
            label celli = pCells[pCelli];

            if (pointLevel[pointi] <= cellLevel[celli])
            {
                // Check if cell has already 8 anchor points -> protect cell
                if (nAnchorPoints[celli] == 8)
                {
                    if (protectedCell.set(celli, true))
                    {
                        nProtected++;
                    }
                }

                if (!protectedCell[celli])
                {
                    nAnchorPoints[celli]++;
                }
            }
        }
    }


    forAll(protectedCell, celli)
    {
        if (!protectedCell[celli] && nAnchorPoints[celli] != 8)
        {
            protectedCell.set(celli, true);
            nProtected++;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMultiFieldRefineFvMesh::dynamicMultiFieldRefineFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    meshCutter_(*this),
    dumpLevel_(false),
    nRefinementIterations_(0),
    protectedCell_(nCells(), 0)
{
    // Read static part of dictionary
    readDict();


    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(nCells(), 0);

    label nProtected = 0;

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        forAll(pCells, i)
        {
            label celli = pCells[i];

            if (!protectedCell_.get(celli))
            {
                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    nAnchors[celli]++;

                    if (nAnchors[celli] > 8)
                    {
                        protectedCell_.set(celli, 1);
                        nProtected++;
                    }
                }
            }
        }
    }


    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner since
    // the coupled one is refined.

    {
        labelList neiLevel(nFaces());

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[faceNeighbour()[facei]];
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[faceOwner()[facei]];
        }
        syncTools::swapFaceList(*this, neiLevel);


        boolList protectedFace(nFaces(), false);

        forAll(faceOwner(), facei)
        {
            label faceLevel = max
            (
                cellLevel[faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = faces()[facei];

            label nAnchors = 0;

            forAll(f, fp)
            {
                if (pointLevel[f[fp]] <= faceLevel)
                {
                    nAnchors++;

                    if (nAnchors > 4)
                    {
                        protectedFace[facei] = true;
                        break;
                    }
                }
            }
        }

        syncTools::syncFaceList(*this, protectedFace, orEqOp<bool>());

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCell_.set(faceOwner()[facei], 1);
                nProtected++;
                protectedCell_.set(faceNeighbour()[facei], 1);
                nProtected++;
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCell_.set(faceOwner()[facei], 1);
                nProtected++;
            }
        }

        // Also protect any cells that are less than hex
        forAll(cells(), celli)
        {
            const cell& cFaces = cells()[celli];

            if (cFaces.size() < 6)
            {
                if (protectedCell_.set(celli, 1))
                {
                    nProtected++;
                }
            }
            else
            {
                forAll(cFaces, cFacei)
                {
                    if (faces()[cFaces[cFacei]].size() < 4)
                    {
                        if (protectedCell_.set(celli, 1))
                        {
                            nProtected++;
                        }
                        break;
                    }
                }
            }
        }

        // Check cells for 8 corner points
        checkEightAnchorPoints(protectedCell_, nProtected);
    }

    if (returnReduce(nProtected, sumOp<label>()) == 0)
    {
        protectedCell_.clear();
    }
    else
    {

        cellSet protectedCells(*this, "protectedCells", nProtected);
        forAll(protectedCell_, celli)
        {
            if (protectedCell_[celli])
            {
                protectedCells.insert(celli);
            }
        }

        Info<< "Detected " << returnReduce(nProtected, sumOp<label>())
            << " cells that are protected from refinement."
            << " Writing these to cellSet "
            << protectedCells.name()
            << "." << endl;

        protectedCells.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMultiFieldRefineFvMesh::~dynamicMultiFieldRefineFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicMultiFieldRefineFvMesh::update()
{
    // Re-read dictionary. Choosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
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
        ).optionalSubDict(typeName + "Coeffs")
    );

    List<dictionary> fieldDicts = List<dictionary>
    (
        refineDict.lookup("fields")
    );

    label numFields = fieldDicts.size();

    //label refineInterval = readLabel(refineDict.lookup("refineInterval"));
    List<label> refineIntervals(numFields);

    forAll(refineIntervals, i)
    {
        refineIntervals[i] = readLabel(fieldDicts[i].lookup("refineInterval"));
    }

    bool hasChanged = false;

    /*if (refineInterval == 0)
    {
        topoChanging(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorInFunction
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }*/
    forAll(refineIntervals, i)
    {
        if (refineIntervals[i] < 0)
        {
            FatalErrorInFunction
                << "Illegal refineInterval " << refineIntervals[i] + " in field " + i << nl
                << "The refineInterval setting in the dynamicMeshDict should"
                << " be >= 1." << nl
                << exit(FatalError);
        }
    }
    forAll(refineIntervals, i)
    {
        if (refineIntervals[i] != 0) break;

        if (i == numFields-1)
        {
            topoChanging(hasChanged);

            return false;
        }
    }



    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    bool timeModRefineIntervalsEqZero = false;

    forAll(refineIntervals, i)
    {
        if (time().timeIndex() % refineIntervals[i] == 0)
        {
            timeModRefineIntervalsEqZero = true;
            break;
        }
    }

    if (time().timeIndex() > 0 && timeModRefineIntervalsEqZero)
    {
        label maxCells = readLabel(refineDict.lookup("maxCells"));

        if (maxCells <= 0)
        {
            FatalErrorInFunction
                << "Illegal maximum number of cells " << maxCells << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        //label maxRefinement = readLabel(refineDict.lookup("maxRefinement"));

        List<label> maxRefinements(numFields);

        label sumMaxRefinements(0);

        forAll(maxRefinements, i)
        {
            maxRefinements[i] = readLabel(fieldDicts[i].lookup("maxRefinement"));
            sumMaxRefinements+=maxRefinements[i];
        }

        /*if (maxRefinement <= 0)
        {
            FatalErrorInFunction
                << "Illegal maximum refinement level " << maxRefinement << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }*/

        forAll(maxRefinements, i)
        {
            if (maxRefinements[i] < 0)
            {
                FatalErrorInFunction
                    << "Illegal maximum refinement level "
                    << maxRefinements[i] + " in field " + i << nl
                    << "The maxCells setting in the dynamicMeshDict should"
                    << " be >= 0." << nl
                    << exit(FatalError);
            }
        }

        //const word fieldName(refineDict.lookup("field"));
        List<word> fieldNames(numFields);

        /*const volScalarField& vFld = lookupObject<volScalarField>(fieldName);

        const scalar lowerRefineLevel =
            readScalar(refineDict.lookup("lowerRefineLevel"));
        const scalar upperRefineLevel =
            readScalar(refineDict.lookup("upperRefineLevel"));
        const scalar unrefineLevel = refineDict.lookupOrDefault<scalar>
        (
            "unrefineLevel",
            GREAT
        );
        const label nBufferLayers =
            readLabel(refineDict.lookup("nBufferLayers"));*/

        scalarList lowerRefineLevels(numFields);
        scalarList upperRefineLevels(numFields);
        scalarList unrefineLevels(numFields);
        List<label> nBufferLayersList(numFields);

        forAll(fieldNames, i)
        {
            fieldNames[i] = word(fieldDicts[i].lookup("field"));


            lowerRefineLevels[i] =
                readScalar(fieldDicts[i].lookup("lowerRefineLevel"));
            upperRefineLevels[i] =
                readScalar(fieldDicts[i].lookup("upperRefineLevel"));
            unrefineLevels[i] =
                readScalar(fieldDicts[i].lookup("unrefineLevel"));
            nBufferLayersList[i] =
                readLabel(fieldDicts[i].lookup("nBufferLayers"));
        }

        // Cells marked for refinement or otherwise protected from unrefinement.
        //PackedBoolList refineCell(nCells());
        List<PackedBoolList> refineCells(numFields);

        // Determine candidates for refinement (looking at field only)
        /*selectRefineCandidates
        (
            lowerRefineLevel,
            upperRefineLevel,
            vFld,
            refineCell
        );*/
        forAll(refineCells, i)
        {
            const volScalarField& vFld =
                lookupObject<volScalarField>(fieldNames[i]);

            selectRefineCandidates
            (
                lowerRefineLevels[i],
                upperRefineLevels[i],
                vFld,
                refineCells[i]
            );
        }

        if (globalData().nTotalCells() < maxCells && (sumMaxRefinements > 0.5))
        {
            PackedBoolList refCells0 = refineCells[0];

            labelList maxRefinementLimit0(nCells(), maxRefinements[0]);

            for(int i(1); i < numFields; i++) // ... orOperator für alle Sublisten!
            {
                PackedBoolList refCells1 = refineCells[i-1];
                PackedBoolList refCells2 = refineCells[i];

                if (i > 1)
                {
                    forAll(refCells0, celli)
                    {
                        refCells0[celli] = refCells0[celli] || refCells1[celli];
                    }
                }

                forAll(refCells2, celli)
                {
                    refCells2[celli] = !refCells0[celli] && refCells2[celli];

                    if ( refCells2.get(celli) )
                    {
                        maxRefinementLimit0[celli] =
                            min(maxRefinementLimit0[celli], maxRefinements[i]);
                    }
                }

                refineCells[i] = refCells2;
            }

            forAll(refineCells, i)
            {
                const labelList& maxRefinementLimit(maxRefinementLimit0);
                const PackedBoolList& refineCell = refineCells[i];

                // Select subset of candidates. Take into account max allowable
                // cells, refinement level, protected cells.
                labelList cellsToRefine
                (
                    selectRefineCells
                    (
                        maxCells,
                        maxRefinementLimit,//maxRefinement,
                        refineCell
                    )
                );

                label nCellsToRefine = returnReduce
                (
                    cellsToRefine.size(), sumOp<label>()
                );

                maxCells -= nCellsToRefine*7;

                if (nCellsToRefine > 0)
                {
                    // Refine/update mesh and map fields
                    autoPtr<mapPolyMesh> map = refine(cellsToRefine);

                    // Update refineCell. Note that some of the marked ones have
                    // not been refined due to constraints.
                    {
                        const labelList& cellMap = map().cellMap();
                        const labelList& reverseCellMap = map().reverseCellMap();

                        PackedBoolList newRefineCell(cellMap.size());

                        forAll(cellMap, celli)
                        {
                            label oldCelli = cellMap[celli];

                            if (oldCelli < 0)
                            {
                                newRefineCell.set(celli, 1);
                            }
                            else if (reverseCellMap[oldCelli] != celli)
                            {
                                newRefineCell.set(celli, 1);
                            }
                            else
                            {
                                //newRefineCell.set(celli, refineCell.get(oldCelli));
                                newRefineCell.set(celli, refineCells[i].get(oldCelli));
                            }
                        }
                        //refineCell.transfer(newRefineCell);
                        refineCells[i].transfer(newRefineCell);
                    }

                    // Extend with a buffer layer to prevent neighbouring points
                    // being unrefined.
                    for (label j = 0; j < nBufferLayersList[i]; j++)
                    {
                        //extendMarkedCells(refineCell);
                        extendMarkedCells(refineCells[i]);
                    }

                    hasChanged = true;
                }
            }
        }


        {
            const scalarList& unrefineLevels2 = unrefineLevels;
            const labelList& maxRefinements2 = maxRefinements;
            List<PackedBoolList> markedCells2 = refineCells;

            List<scalarField> minCellFields(numFields);

            forAll(minCellFields, i)
            {
                minCellFields[i] =
                    minCellField(lookupObject<volScalarField>(fieldNames[i]));
            }

            List<scalarField> pFld2 = minCellFields;

            // Select unrefineable points that are not marked in refineCell
            /*labelList pointsToUnrefine
            (
                selectUnrefinePoints
                (
                    unrefineLevel,
                    refineCell,
                    maxCellField(vFld)
                )
            );*/
            labelList pointsToUnrefine2
            (
                selectUnrefinePoints2
                (
                    unrefineLevels2,
                    maxRefinements2,
                    markedCells2,
                    pFld2
                )
            );

            label nSplitPoints = returnReduce
            (
                pointsToUnrefine2.size(),//pointsToUnrefine.size(),
                sumOp<label>()
            );

            if (nSplitPoints > 0)
            {
                // Refine/update mesh
                unrefine
                (
                    pointsToUnrefine2//pointsToUnrefine
                );

                hasChanged = true;
            }
        }


        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occassionally (how often?).
            // Unrefinement causes holes in the refinementHistory.
            const_cast<refinementHistory&>(meshCutter().history()).compact();
        }
        nRefinementIterations_++;
    }

    topoChanging(hasChanged);
    if (hasChanged)
    {
        // Reset moving flag (if any). If not using inflation we'll not move,
        // if are using inflation any follow on movePoints will set it.
        moving(false);
    }

    Info << "numberOfCells: " << globalData().nTotalCells() << endl;

    return hasChanged;
}


bool Foam::dynamicMultiFieldRefineFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    // Force refinement data to go to the current time directory.
    const_cast<hexRef8&>(meshCutter_).setInstance(time().timeName());

    bool writeOk =
    (
        dynamicFvMesh::writeObject(fmt, ver, cmp, valid)
     && meshCutter_.write(valid)
    );

    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("level", dimless, 0)
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(cellLevel, celli)
        {
            scalarCellLevel[celli] = cellLevel[celli];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }

    return writeOk;
}


// ************************************************************************* //
