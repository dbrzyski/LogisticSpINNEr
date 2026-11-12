# Fig. 5: colorbar

The glass-brain panels of Fig. 5 were rendered with BrainNet Viewer, which colors each edge by its median
bootstrap coefficient on the jet colormap, from the smallest (−0.1210) to the largest (0.1433) coefficient of the
41 cortical edges. Its own colorbar is unreadable at print size, so the manuscript crops it and places a redrawn one
next to the brains.

`make_fig5_colorbar.m` reads `../fig5_fig6_real_data/results/EstimatesBootstrap.mat` (create it with the scripts
in that folder; the data are not included), selects the edges of Figs. 5–6 (nonzero in at least 90% of the 300
bootstrap fits, both regions cortical), checks that there are 41 of them and that their median coefficients span
−0.1210 to 0.1433, and writes the colorbar as the vector file `fig5-colorbar.pdf` in this folder. It needs only
base MATLAB.
