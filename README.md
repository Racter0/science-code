science-code
============

Code written at grad school

This is the initial iteration of the code that eventually became FindDQE.

This algorithm uses the "edge method" to find the DQE of an electron detector from an image of the beamstop. Requires the image to be in MRC format, preferably the MRC format as output by IMAGIC's em2em routine.

I wrote FindDQE.m; other m-files included in the repository, e.g. the spline-fitting algorithm, are used by FindDQE.m but are from the MATLAB exchange and were not written by me.

FindDQE, as found on my lab's website grigoriefflab.janelia.org, works more reliably due to using an innovative algorithm rather than the edge method. This code is just deposited here as a sample of my MATLAB coding.
