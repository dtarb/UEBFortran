
# UEB Grid Snowmelt Model

## Version 1.0

David Tarboton, Avirup Sen Gupta

March 23 2014

####Folders

* Model:  This contains the source code and Compiled Windows Binary executable in Release folder

* Data:  This contains example input and output data sets

* Documentation:  This contains documentation



The Utah Energy Balance (UEB) snow model is an energy balance snowmelt model developed by David Tarboton's research group, first in 1994, and updated over the years. The model is written in FORTRAN and uses a lumped representation of the snowpack and keeps track of water and energy balance. The model is driven by inputs of air temperature, precipitation, wind speed, humidity and radiation at time steps sufficient to resolve the diurnal cycle (six hours or less). The model uses physically-based calculations of radiative, sensible, latent and advective heat exchanges. A force-restore approach is used to represent surface temperature, accounting for differences between snow surface temperature and average snowpack temperature without having to introduce additional state variables. Melt outflow is a function of the liquid fraction, using Darcy's law. This allows the model to account for continued outflow even when the energy balance is negative. Because of its parsimony (few state variables - but increasing with later versions) this model is suitable for application in a distributed fashion on a grid over a watershed. There are a number of versions available. 

This repository holds the UEBGrid version, developed for gridded application of UEB. Other versions of UEB may be obtained from [http://hydrology.usu.edu/dtarb/snow/snow.html](http://hydrology.usu.edu/dtarb/snow/snow.html)

UEBGrid added the capability to represent the melting of glaciers and adopted a very structured file based input/output format using ASCII and netCDF files to facilitate its use in a NASA project and its incorporation into the EPA BASINS software.  

UEBGrid is a command line program developed using Intel Fortran XE 2011 on a Windows platform.  We have done our best to use standard platform independent Fortran so the code should work with other compilers.  We welcome comments on making this better.

To run the model

1. Open a command prompt.

2. Create a path to the executable UEB\Model\UEBGrid\Release\UEBGrid.exe (or copy the executable to the folder where input data is saved)

3. Change directory to the folder where the input data is saved

4. At the prompt type: uebgrid twdefControl.dat

Project configuration

All source code is in the folder UEB\Model\

The file UEB\Model\VisualStudioProjConfiguration.txt describes the configuration of UEB for compiling with Intel Fortran XE 2011 and Visual Studio.  This also includes information on obtaining and linking to the Unidata NetCDF library that we use.

The folder UEB\Model\UEBGrid contains the Visual Studio Solution and project files.  
UEB\Model\UEBGrid\Release contains the latest compiled executable for Windows (32 bit).

Copyright (C) 2014  David Tarboton, Utah State University, dtarb@usu.edu.  http://hydrology.usu.edu/dtarb/ 

UEB is open source software: you can redistribute it and/or modify it under the terms of the MIT Open Source License as published by the Open Source Initiative https://opensource.org/licenses/MIT.

UEB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

If you wish to use or incorporate this program (or parts of it) into other software that does not meet the MIT License conditions contact the author to request permission.

David G. Tarboton  
Utah State University  
8200 Old Main Hill  
Logan, UT 84322-8200  
USA  
http://hydrology.usu.edu/dtarb/ 

email:  dtarb@usu.edu 

Acknowledgements

I am grateful to the following funding agencies that have supported the develoopment of UEB

* NASA Grant NNX11AK03G has supported the development of UEBGrid

* USDA-CREES award 2008-34552-19042 Utah Drought Management Project supported the development of the vegetation components.
