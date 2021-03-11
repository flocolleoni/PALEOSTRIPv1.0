Flex3D: Mathematically 2D, Physically 3D, deflection of a plate of variable thickness

Eq. 3.83 of Ventsel and Krauthammer, 2001: Thin plates and shells. Solution by centered finite differences

How does the program work?

In order to run the program you need two input text files: One file with the height (in meters) and density (in kg/m^3) of the loads (loads.txt), and another file with the distribution of elastic thickness in meters (te.txt). Loads and elastic thicknesses should be given at the nodes or points of a regular grid with equal x, y spacing (delta) between the nodes.

Besides these files, you need to edit some parameters in flex3d.m before running your problem. These parameters are:

Number of nodes or points in x
Number of nodes or points in y
Distance between nodes or points (delta)
Young Modulus in Pascals
Poisson's ratio
Gravity in m/s^2
Density of the mantle in kg/m^3
Density of the material filling any resultant depression in kg/m^3

That is all you need. After that, just type flex3d in the Matlab command window and the program will output the deflections (in meters) at the nodes (file deflection.txt)

IMPORTANT: BOUNDARY CONDITIONS: Displacement at plate boundaries is zero. Therefore you should locate your loads in the middle of the plate and have enough space between the loads and the plate boundaries to avoid side, boundary effects.

EXAMPLES:

Several examples are provided to show the functionality of the program and to prove that the program works. This examples are test1.m to test7.m. The way you run these examples is as follows:

In the Matlab command window Type:

clear all
test1.m
clear all
flex3d
clear all
plotTest

After this, you should get a map showing the deflection in meters. Do the same for all tests. A brief description of the load and elastic thickness distribution is included in each test1-7.m file. Make sure to read this. The test will give you a clear idea of how to write the input files for flex3d, and in general how to set up the program parameters.

Test 7 includes an irregular distribution of loads and elastic thicknesses, and therefore should be executed within its own folder (test7). Just cd to that folder. Type

test7

This will produce two plots with the interpolated loads and elastic thicknesses from their contour maps (load.jpg and te.jpg), and the files necessary for flex3d. The interpolation is done using script gridfit by John D'Errico (2006). Notice that the loads are put in the center of the plate, and some space is left between the loads and the plate boundaries to avoid side effects.

Then type

clear all
flex3d
clear all
plotTest

You will be presented with a map of the deflection.

The scripts are based on a numerical approximation, and therefore there might be problems with convergence, grid size, etc. You should be careful when using the scripts. Don't use them blindly as a black box.

The scripts (except for gridfit) are copyright of Nestor Cardozo 2009. They can only be used for academic/research purposes and they cannot be distributed by third parties. Nestor Cardozo assumes no liability for damages, direct or consequential, which may result from the use of the scripts.