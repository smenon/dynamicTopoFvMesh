{
    const labelList& curPointEdges = pointEdges[curPoint];

    label patchID = -1;
    label edgeID = -1;

    forAll (curPointEdges, edgeI)
    {
        label curEdge = curPointEdges[edgeI];

        if (edgeFaces[curEdge].size() == 1)
        {
            forAll(aMesh().boundary(), patchI)
            {
                const labelList& curEdges = aMesh().boundary()[patchI];

                forAll(curEdges, eI)
                {
                    if (curEdges[eI] == curEdge)
                    {
                        if 
                        (
                            aMesh().boundary()[patchI].type()
                         == emptyFaPatch::typeName
                        )
                        {
                            patchID = patchI;
                            edgeID = curEdge;
                            break;
                        }
                    }
                }
            }
        }
    }

    if (patchID != -1)
    {
        vector nE = 
            pointNormals[edges[edgeID].start()]
          + pointNormals[edges[edgeID].end()];

        nE /= mag(nE);

        vector eP =
            controlPoints()[edgeFaces[edgeID][0]]
          - edges[edgeID].centre(points);

        vector mirrorCtrlPoint =
            edges[edgeID].centre(points)
          + ((2.0*nE*nE - I)&eP);

        // Resize weight field
        w.setSize(w.size()+1);

        w[w.size()-1] = 1.0/mag
            (
                mirrorCtrlPoint
              - points[curPoint]
            );

        // Add mirror control point to LS matrix and rhs vector
        ownM[curPatchPoint] += 
            sqr(w[w.size()-1])*sqr(mirrorCtrlPoint);

        ownS[curPatchPoint] += 
            sqr(w[w.size()-1])*mirrorCtrlPoint;
    }
}
