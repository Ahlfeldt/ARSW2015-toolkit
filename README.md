# The Economics of Density: Evidence from the Berlin Wall: Toolkit

© Gabriel M. Ahlfeldt, S. J. Redding, D. Sturm, N. Wolf

Version 0.91, 2024

## General Instructions

This toolkit complements the original replication directory for Ahlfeldt, Redding, Sturm, Wolf (2015). The toolkit does not cover all stages of the analysis presented in the article. Instead, it covers a subset of codes that are crucial for the quantification and simulation of the model using 2006 data. Detailed comments have been added to all scripts to provide intuition and closer links to the model described in the article and supplement. Some codes have been processed to make them more accessible, reducing the memory requirements and facilitating faster execution on desktop computers. For the same reason, we use trimmed versions of the datasets, which contain variables that are essential for the purpose of this teaching directory. Any data set or program that has been substantively altered (other than adding comments) is marked by a "_TD" suffix in the filename.

A program has been added to rescale productivities and amenities recovered using the sequential procedure (using files in the section6/calibration folder) so that they are amenable to the equilibrium solvers used in the counterfactual folders of Sections 6 and 7. Some other smaller programs have been added to make the structure of the code more accessible. A META script has been added from which all other scripts are called. Codes have been updated so that the root directory of this toolkit has to be defined only once in the META script. A MAPIT program has been added that allows for immediate inspection of the outputs as they are being generated by the code (alongside shapefiles indexed in the same way as the MATLAB data). Finally, throughout the various steps of the calibration and simulation, illustrative explorations of recovered productivities and amenities, solved wages, and some simple counterfactuals have been added for didactic purposes.

Before you can start, you need to install the following MATLAB toolboxes:

- Global optimization toolbox
- Mapping toolbox
- Statistics and machine learning toolbox

To install a MATLAB toolbox, open MATLAB and go to the 'Home' tab. Click on 'Add-Ons' in the menu. From there, browse or search for the toolbox you want to install. Click on the toolbox to view its details. Then, proceed by clicking 'Install.' Follow the on-screen instructions, which will include accepting the license agreement and selecting the installation path. After the installation is complete, the toolbox will be ready to use. Make sure you have the appropriate licenses for the toolbox you've selected.

You must also download three large files from external URLs indicated in the file directory below.
