# MetapopulationSynchrony

Codes and algorithm used for generating results of manuscript "Dispersal, synchrony, and climate-induced tipping: Mechanisms of persistence and extinction in spatially structured predator-prey metapopulations" by S. Marick et al.

% File Basin of attraction

"Basinofattraction.m" plot of basin boundary with respect to reference domain. This code is used for generating basin boundaries or saddle-separatrix used in Figures where saddle-separatrix is used. 

% File Single patch realization

"single_patch.m" generates and plots realizations of single patch simulation with climatic fluctuations. If a particular realization tips to extinction it provides the tipping time.

% File Two patch realizations

"Two_patch.m" generates realizations of two patch simulation with non-identical climatic fluctuations. The random seed for generating Figure 4 are rng(365) (for fig 4abcd), rng(251) (for fig 4efgh) and rng(2) (for fig 4ijkl).

"Synchrony_plotter.m" generates correlation data of Figure 3b.

% File Relative basin size

"Two_patch_RelbasinSize.m" generates data of relative basin size in Figure 3a. Can be upgraded for relative basin sizes for multipatch system of 100 patches.
 
% File Multipatch realization

"simulation.m" generates realizations of Figure 6 in main text and Video files of multipatch system of supplementary data. Basin boundary data are regenerated before starting this simulation. For generating video of realizations in Fig.6, user must generate the basin boundary data as structured in the "Basin_data.mat", "Basin_data_zoomed" and "Climate_timeseries" data files. For arbitrary realization choose the random seed rng(0). For simulation on arbitrary network structure user must input the respective Laplacian matrix of the network instead of L in line 46.

% File Networks

"scale_free.m" and "erdos_reyni" generates respective random networks and gives Laplacian and adjacency matrix of the generated network as the output.

% File multipatch frequency
To generate data for frequency of events in multipatch scenario. Four cases are defined as accordingly.
