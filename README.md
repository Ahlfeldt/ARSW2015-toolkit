# The Economics of Density: Evidence from the Berlin Wall - toolkit

**© Gabriel M. Ahlfeldt, S. J. Redding, D. Sturm, N. Wolf**

Version 0.91, 2024

## General Instructions

This toolkit complements the original replication directory for Ahlfeldt, Redding, Sturm, Wolf (2015). The toolkit does not cover all stages of the analysis presented in the article. Instead, it covers a subset of codes that are crucial for the quantification and simulation of the model using 2006 data. Detailed comments have been added to all scripts to provide intuition and closer links to the model described in the article and supplement. Some codes have been processed to make them more accessible, reducing the memory requirements and facilitating faster execution on desktop computers. For the same reason, we use trimmed versions of the datasets, which contain variables that are essential for the purpose of this teaching directory. Any data set or program that has been substantively altered (other than adding comments) is marked by a "_TD" suffix in the filename.

A program has been added to rescale productivities and amenities recovered using the sequential procedure (using files in the section6/calibration folder) so that they are amenable to the equilibrium solvers used in the counterfactual folders of Sections 6 and 7. Some other smaller programs have been added to make the structure of the code more accessible. A META script has been added from which all other scripts are called. Codes have been updated so that the root directory of this toolkit has to be defined only once in the META script. A MAPIT program has been added that allows for immediate inspection of the outputs as they are being generated by the code (alongside shapefiles indexed in the same way as the MATLAB data). Finally, throughout the various steps of the calibration and simulation, illustrative explorations of recovered productivities and amenities, solved wages, and some simple counterfactuals have been added for didactic purposes.

Before you can start, you **need to install the following MATLAB toolboxes:**

- Global optimization toolbox
- Mapping toolbox
- Statistics and machine learning toolbox

To install a MATLAB toolbox, open MATLAB and go to the 'Home' tab. Click on 'Add-Ons' in the menu. From there, browse or search for the toolbox you want to install. Click on the toolbox to view its details. Then, proceed by clicking 'Install.' Follow the on-screen instructions, which will include accepting the license agreement and selecting the installation path. After the installation is complete, the toolbox will be ready to use. Make sure you have the appropriate licenses for the toolbox you've selected.

**You must also download three large files from external URLs** indicated in the file directory below.

Note that due to large matrices, the execution of this teaching directory is memory-intensive. Even with large internal memory (RAM), you may receive out-of-memory error messages. To avoid this, ensure you **assign plenty of virtual memory** (e.g., 32GB) on your system, preferably on a solid-state drive. Using virtual memory will increase computational time, but it will at least ensure that the codes execute successfully.

## How to Use the Toolkit

The toolkit is designed as a didactic journey through the codes essential for the quantification and the simulation of the model. The aim is to convey how to quantify the model and conduct simple counterfactuals.

Once you are set up, execute the codes in the order in which they are being called by `META.m`. When you use the toolkit for the first time, it is important that you execute the scripts exactly in the order outlined by `META.m` as intermediate inputs will be generated that will be used at later stages. There are also didactic reasons for proceeding in the given order since the structure broadly follows the use of the model in the paper and codes become progressively more complex. For the best learning experience, we recommend that you open the scripts and go through the code line by line (instead of calling scripts from `META.m`). Plenty of comments have been added to refer to the relevant equations in the paper. It will likely be beneficial to triangulate between the code, the pseudo codes summarized in the **codebook**, and the paper and the supplement to understand how the code implements the model.

Throughout the quantification of the model, we use the `MAPIT` program to illustrate selected variables recovered using the structure of the model. You can use `MAPIT` to conveniently plot any model input or output at any stage of the analysis to develop your intuition for the various endogenous and exogenous objects. Due to the granularity of the data, `MAPIT` is slow. If you want to experiment with the code and reduce the computational time, you may find it convenient to outcomment the use of `MAPIT`.

After recovering the unobserved exogenous location characteristics, we use various solution algorithms to solve for the spatial equilibrium. For illustrative purposes, we perform two simple counterfactuals. First, we implement a ban of private cars by setting bilateral travel times to those that arise when using public transport (this counterfactual is also performed in the paper). Second, we increase the fundamental productivity within the boundaries of former East Berlin. For each counterfactual, we compare the effects on workplace employment, residence employment, wages and floor space prices between the cases with exogenous productivity and amenity and endogenous agglomeration forces. We provide a similar comparison for the closed-city and open city-case. You can use `MAPIT` to inspect the effects on any other endogenous outcome and execute similar counterfactuals by changing other primitives of the model. The toolkit will hopefully be sufficiently transparent for you to conveniently run your own counterfactuals.

When using the toolkit in your work, please cite Ahlfeldt, Redding, Sturm, Wolf (2015)

## MATLAB Files and Data Folders

| Directory | File | Description | Additional Information |
| --- | --- | --- | --- |
| `matlab/data/input` | | Folder containing required data inputs to execute this teaching directory | -|
| `matlab/data/input` | `prepdata_big_TD86.mat` | MATLAB file containing 1986 variable generated using `matlab/prepdata_TD86.m` code. | This is a large data set, you need to download it to your machine from an external URL. [Download](https://box.hu-berlin.de/f/3f364937b4a240e883f8/?dl=1) |
| `matlab/data/input` | `prepdata_big_TD.mat` | MATLAB file containing 2006 variable generated using `matlab/prepdata_TD.m` code. | This is a large data set, you need to download it to your machine from an external URL. [Download](https://box.hu-berlin.de/f/54d2f718ec8644e5888f/?dl=1) |
| `matlab/data/input` | `ttpublic_2006_ren.mat` | MATLAB file generated by `prepdata_TD06ftpub.m` containing travel times without cars used for counterfactuals. | This is a large data set, you need to download it to your machine from an external URL. [Download](https://box.hu-berlin.de/f/9eb052b712204866beda/?dl=1) |
| `matlab/data/input` | `roptimis_all.mat` | MATLAB file containing parameter values estimated using the GMM procedure in Section 7. | -|
| `matlab/data/output` | | Folder containing data files generated by the scripts in this toolkit. | Will be populated while you execute the MATLAB programmes. |
| `matlab/figs` | | Folder containing various subfolders that will be populated as you execute the MATLAB programmes. |-|
| `matlab/section6` | | Folder containing various subfolders and code for the quantification and simulation of the model with exogenous fundamentals. |-|
| `matlab/section7` | | Folder containing various subfolders and code for the quantification and simulation of the model with endogenous agglomeration forces. |-|

## MATLAB Programs
> We highlight main programs in **bold** that call other functions to execute the analyses.

| Script | Description | Special Instructions |
| --- | --- | --- |
| **`matlab/META.m`** | Meta file that calls other code files to execute the analysis. Your journey through the teaching directory starts here! | -|
| `matlab/prepdata_TD.m` | Loads various data files provided in the replication directory. Many files are not read in this teaching directory to save space. Output is saved as Matlab/data/prepdata_big_TD.mat. | You do not need to execute this programme unless you want to recreate the `prepdata_big_TD.mat` file from data provided in the replication directory. |
| `matlab/prepdata_TD86.m` | Loads various data files provided in the replication directory. Many files are not read in this teaching directory to save space. Output is saved as Matlab/data/prepdata_big_TD86.mat. | You do not need to execute this programme unless you want to to recreate the `prepdata_big_TD86.mat` file from data provided in the replication directory. |
| `matlab/prepdata_TD06ftpub.m` | Reads 2006 travel times without cars and saves them as ttpublic_2006_ren.mat for use in counterfactuals. | You do not need to execute this programme unless you want to recreate the `ttpublic_2006_ren.mat` file from data provided in the replication directory. |
| **`matlab/MAPIT`** | Function that can be called to create simple maps that illustrate outcomes by block. Useful to develop an intuition for the variables that are being generated and how they relate to each other economically. Especially if you are familiar with the geography of Berlin. | -|

**One-step epsilon estimation**

| Script | Description | Special Instructions |
| --- | --- | --- |
 `matlab/section6/optimepsilon/optimepsilon_TD86.m` | Script containing the steps involved in the one-step estimation of epsilon. | -|
|  `matlab/section6/optimepsilon/comegaopt0.m` | Function containing the solver that solves for transformed wages omega. | -|
|  `matlab/section6/optimepsilon/cdensityoptren.m` | Objective function that is being minimized when estimating epsilon. | -|
| `matlab/section6/optimepsilon/cdensityoptren.m` | Small program used to keep selected objects |  -|

**Section 6 quantification with exogenous amenities (sequential procedure)** 

| Script | Description | Special Instructions |
| --- | --- | --- |
| `matlab/section6/calibration/calcal_TD.m` | Script containing the steps involved in the quantification. | -|
|   `matlab/section6/calibration/calcal_adj_TD.m` | Program that rescales amenities and productivities so that they generate a population in the data that matches the model. | This program has no equivalent in the original replication directory. It only needs to be run if you are interested in the level of productivities and amenities. It must be run if you want to use the amenities and productivities recovered in this section in the smodemx.m solver in Section 6. | 
|   `/matlab/section6/calibration/comegaoptC.m` | Function that recovers adjusted productivities and adjusted wages. | -|
|  `/matlab/section6/calibration/camen.m` | Function that recovers commuting adjusted amenities. | -|
|  `/matlab/section6/calibration/expincome.m` | Function that recovers total expected worker income. | -|
|  `/matlab/section6/calibration/cdensity.m` | Function that recovers density of development, total floor space, residential floor space share. | -|
|  `/matlab/section6/calibration/modbezirk.m` | Function that generates modern Bezirke identifier. | -|

**Section 6 Counterfactuals with Exogenous Amenities**

| Script | Description | Special Instructions |
| --- | --- | --- |
| **`matlab/section6/exogcfutal/cfprep_TD.m`** | Prepares data for the counterfactuals with exogenous fundamentals (generates exogcfutal_prep_big_TD.mat). Illustrates various recovered variables and performs a comparison to the fundamentals recovered using the sequential procedure in the section6/calibration folder. |  -|
| `matlab/section6/exogcfutal/cmodexog.m` | Function that inverts amenities and productivities and solves for equilibrium wages given observed data simultaneously. Called by cfprep_TD.m. | -|
| `matlab/section6/calibration/cdensityE.m` | Function that recovers density of development, total floor space, commercial floor space, residential floor space, and commercial floor space share. Called by cfprep_TD.m. |- |
| **`matlab/section6/exogcfutal/cftualexog_TD`** | Runs illustrative counterfactuals and illustrates the effects on endogenous outcomes. | -|
| `matlab/section6/exogcfutal/smodex` | Program that solves for the equilibrium for given primitives under exogenous fundamentals in the closed-city case (exogenous employment).  | -|

**Section 7 Counterfactuals with Endogenous Agglomeration Forces**

| Script | Description | Special Instructions |
| --- | --- | --- |
| `cftualprep_end_TD.m` | Prepares data for the counterfactuals with endogenous fundamentals (generates matlab/data/endogcfual_prep_big_TD.mat). Illustrates various recovered variables and performs a comparison to the fundamentals recovered using the sequential procedure in the section6/calibration folder. |  -|
| `Cftalendog_HHbar_TD.m` | Conducts illustrative counterfactuals under endogenous agglomeration forces in a closed city and compares the results to the case with exogenous fundamentals. | - |
| `Cftalendog_Ubar_TD.m` | Conducts illustrative counterfactuals under endogenous agglomeration forces in an open city and compares the results to the closed-city case. |  -|
| `cprod.m` | Function that decomposes productivities recovered by cmodexog.m into endogenous and exogenous components. |  -|
| `cres.m` | Function that decomposes amenities recovered by cmodexog.m into endogenous and exogenous components. |  -|
| `ubar.m` | Function that recovers the exogenous reservation utility level. |  -|
| `smodendog.m` | Function that solves for the equilibrium under endogenous agglomeration forces in the closed-city case (exogenous total employment). |  -|
| `ussmodendog.m` | Function that solves for the equilibrium under endogenous agglomeration forces in the open-city case (exogenous utility). |  -|

## Shapefiles

| Name | Description |
| --- | --- |
| `BerlinAllBlocks` | Full block shapefile from Senatsverwaltung (public administration) |
| `Berlin4matlab` | Trimmed version of the BerlinAllBlocks shapefile, containing only blocks with workplace or residence employment, indexed in the same way as MATLAB data set |
| `WestBerlin4matlab` | Same as Berlin4matlab, but restricted to blocks in former West Berlin |
| `BerlinGreen` | Shapefile of green spaces |
| `BerlinWater` | Shapefile of water bodies |

## Further resources:

Ahlfeldt, Redding, Sturm, Wolf (2015): The Economics of Density: Evidence from the Berlin Wall, 83(6), p. 21272189. [https://doi.org/10.3982/ECTA10876](https://doi.org/10.3982/ECTA10876). The journal website provides the link to the full replication directory.

## Acknowledgements: 

We thank Alexander Hansen and Andrea Herrera for comments and suggestions.

