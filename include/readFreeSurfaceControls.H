const dictionary& fs = mesh.solutionDict().subDict("FreeSurface");

Switch smooth
(
    fs.lookup("smooth")
);

// Switch ddtPhiCorr
// (
//     fs.lookup("ddtPhiCorr")
// );

int nFreeSurfCorr = 0;
if (fs.found("nFreeSurfCorr"))
{
    nFreeSurfCorr = readInt(fs.lookup("nFreeSurfCorr"));
}
