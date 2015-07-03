adapted code from evodevocode: inclusive of all analysis, features and bugs!

Header:
- GrowGeneNr
- NrFinalCells ->150 (if grow even further than this, you cannot divide, but should also be stopped)


Agent.cc:

- removed intercalation stuff: no longer needed

FormZygote:
- can still specify how many cells to start with. 
- posterior-most cell receives morphogen 100.0

DivideGrowthZone -> DivideCells()
- Now have to check whether a cell divides: if the concentration of the growth gene is above a threshold, you may divide. Can either do this with some probability, or a clock-like system. Can start with the probability and halving the growgene, but perhaps a timer is more desirable: then you do not have to evolve the regulation for the growth gene as thoroughly.
Indeed, now at the beginning of the run, the animal does not regulate the gene as thoroughly, and the tissue explodes unreasonably fast. Instead: add division counter to cell. After x steps, a cell has a high chance of dividing, after which the division counter is reset. 
UPDATE: a threshold works fine, it just has to be higher!

WriteDivisionProfile: makes a picture of cell age.

CellCellSignalling:
could remove the morphogen gradient part: at least shut down the Gradient def.
Instead, added a few other possibilities: Posterior cell high, with or without diffusion (DifCoef, Nrdiffsteps), or free morphogen: only initial bias.
Gradient can still be tested, but will only yield an effect if transition to posterior growth. Can we see that bias?


DetermineFitness:
- added fitness for size
- new parameters: targetsize, sizebonus (benefit for cell gained), sizepen (penalty for overshooting target)
- removed skipzone: no longer a predefined growth zone
- optional penalty for divisions after you start checking for stability: stablesizepen. 

MaintenanceIntraCellularDynamics:
- For now, a daughter cell created within the time window of checking will simply inherit her mother's maint and varmaintproteinstates. (in DivideCells). 
- A penalty for dividing within this window could be added, but I am not sure about the desirability of this. Still, I store the size of the agent at that point: maintsize


Network.cc:
UpdateNetworkState: the maternal gene can now be modified by other genes in the network, if in header FREEMORPH is defined.


Runge Kutta integration added; type 4 used in UpdateNetworkState.

To add:
* speed selection
* decay evolution?

Code issues:


