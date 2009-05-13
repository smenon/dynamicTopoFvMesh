// Rotate all points belonging to the 'patchNames' by 'theta'
// about an axis p1-p2, translate by a distance defined by 't'
void rotatePoints
(
    dynamicTopoFvMesh& mesh,
    wordList& patchNames,
    doubleScalar theta,
    vector p1,
    vector p2,
    vector t
)
{
    label nP = 0;
    vector p, q;
    doubleScalar costheta = Foam::cos(theta);
    doubleScalar sintheta = Foam::sin(theta);

    labelList patchID(patchNames.size(), -1);
    labelHashSet pointSet;

    // Define the rotation axis and normalize it
    vector r = (p2-p1)/mag(p2-p1);

    // Fetch the mesh-points
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const pointField& oldMeshPoints(mesh.points());

    // Copy existing point locations
    pointField meshPoints(oldMeshPoints);

    // Match patch names and add to a HashSet to avoid moving
    // patch points twice
    forAll(patchNames, wordI)
    {
        forAll(mesh.boundaryMesh(), patchI)
        {
            if (bMesh[patchI].name() == patchNames[wordI])
            {
                patchID[nP++] = patchI;

                // Add all points of this patch
                const labelList& patchPoints = bMesh[patchI].meshPoints();

                forAll(patchPoints, index)
                {
                    if (!pointSet.found(patchPoints[index]))
                    {
                        pointSet.insert(patchPoints[index]);
                    }
                }

                break;
            }
        }
    }

    // Move all patch points cumulatively
    labelList allPatchPoints = pointSet.toc();

    forAll(allPatchPoints, index)
    {
        q = vector::zero;

        // Fetch the old point and translate it
        p = oldMeshPoints[allPatchPoints[index]];

        // Translate to the origin
        p -= p1;

        // Apply the rotation matrix
        q.x() += (costheta + (1 - costheta) * r.x() * r.x()) * p.x();
        q.x() += ((1 - costheta) * r.x() * r.y() - r.z() * sintheta) * p.y();
        q.x() += ((1 - costheta) * r.x() * r.z() + r.y() * sintheta) * p.z();

        q.y() += ((1 - costheta) * r.x() * r.y() + r.z() * sintheta) * p.x();
        q.y() += (costheta + (1 - costheta) * r.y() * r.y()) * p.y();
        q.y() += ((1 - costheta) * r.y() * r.z() - r.x() * sintheta) * p.z();

        q.z() += ((1 - costheta) * r.x() * r.z() - r.y() * sintheta) * p.x();
        q.z() += ((1 - costheta) * r.y() * r.z() + r.x() * sintheta) * p.y();
        q.z() += (costheta + (1 - costheta) * r.z() * r.z()) * p.z();

        // Translate back to original location
        q += p1;

        // Assign to the mesh
        meshPoints[allPatchPoints[index]] = q + t;
    }

    // Update the displacement BCs for mesh motion
    setMotionBC(mesh, meshPoints);
}
