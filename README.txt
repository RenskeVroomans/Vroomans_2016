This code was used to produce the results in "Vroomans RMA, Hogeweg P, Ten Tusscher KHWJ. 2016. In silico evo-devo: reconstructing stages in the evolution of animal segmentation. EvoDevo 7:14

The code depends on the following libraries: libpng-dev and zlib1g-dev.
It was developed in 2015, so the makefile may be a bit old-fashioned.

After installing dependencies, simply type "make"; the binary will be placed in bin/
The second makefile, specifying flags and libraries, can be found in obj/
It contains a few flags that can be specified to compile different versions: the main simulation code and different programs to analyse simulation outcomes.

