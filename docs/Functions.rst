Functions
=========

#######
core.py
#######

getComponentData
****************
Inputs:
    *   componentList

In order to simplify the process of importing component data from source files, this function will
automatically generate the correct 2-dimensional array containing the fluid-specific data for
each component specified in the 1-dimensional componentList list input.

getInteractionParameters
************************
Inputs:
    *   componentList
  
This function serves the same purpose as getComponentData, but it does so for Peng-Robinson binary
interaction paramters, also as a 2-dimensional array.

guessPressure
*************
Inputs:
    *   componentData
    *   moleFractions
    *   T
  
Using inputs of a 2-dimensional componentData list (of the same shape provided by the
getComponentData function), the moleFractions (of the same shape provided to the getComponentData
function), and the desired temperature (in K), this function will return a guess pressure in Pa.

guessTemperature
****************
Inputs:
    *   compounds
    *   moleFractions
    *   P
  
This function serves the same purpose as guessTemperature, but it does so for a desired pressure
(in Pa), returning a guess temperature in K.

PRSV
****
Inputs:
    *   compoundData
    *   moleFractions
    *   T
    *   P
    *   interactionParameters
  
Peng-Robinson-Stryjek-Vera equation of state - inputs of 2-dimensional compoundData list,
1-dimensional moleFractions list, T (in K), P (in Pa), and 2-dimensional interactionParameters
list derived from getInteractionParameters function. If κ parameters cannot be found, the function
reverts to a standard Peng-Robinon 1976 equation of state.

hydrationNumber
***************
Inputs:
    *   structure
    *   occupancies
  
The hydration number of a hydrate system is the number of water molecules per guest molecule. This
function takes the structure (as either "I" or "II") and the hydrate occupancies, which have a 
generalized structure of [[θsmall_1, θsmall_2, ...], [θlarge_1, θlarge_2, ...]].

hydrateDensity
**************
Inputs:
    *   structure
    *   occupancies
    *   compoundData
    *   moleFractions
    *   T
    *   P
  
To obtain an estimation of the density of a hydrate system (in kg/m^3), this function is employed
using the same conventions as previously mentioned.

HuLeeSum
********
Inputs:
    *   T
    *   saltConcs
    *   inhibitorConcs
    *   betaGas
    *   freezingPoint
  
This equation utilizes the Hu-Lee-Sum inhibition model for hydrates to obtain the depressed
equilibrium temperature for a given system and salt and/or inhibitor concentration. T is the 
undepressed equilibrium temperature (in K), saltConcs and inhibitorConcs are 1-dimensional lists
using the weight concentrations of the salts (based on water amount) and organic inhibitors (based 
on salt aqueous solutions) given in Component and Inhibitor IDs. betaGas is a value which is
dependent on the guest molecules, and can be determined either from the original Hu-Lee-Sum
publication or obtained with the betaGas function. freezingPoint is an optional variable used as a 
failsafe to ensure that the function does not provide an inhibited temperature result for a 
hydrate which does not form in water under uninhibited conditions. To override this failsafe, 
simply provide an unreasonably low freezing point (such as 200K).

betaGas
*******
Inputs:
    *   temperatures
    *   pressures
  
Using equilibrium temperatures (in K) and pressures (in Pa), this function returns the value for
betaGas for a given equilbrium curve. It suffers from tangible inaccuracy above about 20%
inhibitor weight concentration, so it is advised to pull a value from the original Hu-Lee-Sum
publication if possible, which are fitted to experimental data.

########
model.py
########

getHydrateCellProperties
************************
Inputs:
    *   structure
  
Obtains data about hydrate structures (e.g. number of water molecules per void). Called inside
functions and so generally remains unused; most data is left over from an earlier model iteration.

henrysLawConst
**************
Inputs:
    *   compoundData
    *   T
  
Calculates the Henry's Law Constant (in Pa) based on guest gases dissolved in water using a
4-constant exponential model. compoundData is used as a stand-in for a 1-dimensional slice of 
a standard componentData 2-dimensional list for a single compound.

Z
***
Inputs:
    *   compoundData
    *   T
    *   P
  
Returns the infinite compressibility gas constant; only used in liquid phase composition
calculations.

liqPhaseComposition
*******************
Inputs:
    *   compounds
    *   T
    *   fug_vap
    *   compoundData
    *   P
  
Calculates the composition of the aqueous phase, returning a 1-dimensional list where the index 0
is the fraction of water in the aqueous phase, and subsequent incides contain compositions in the
order in which they were provided in compounds. fug_vap is the vapor fugacity obtained from the
PRSV function, in Pa.

activityCoeff
*************
Inputs:
    *   T
    *   phaseComposition
    *   chemGroups
  
Specially tailored function to return the activity coefficient of the aqueous phase considering
dissolved guest gases utilizing the UNIFAC model from the thermo library. phaseComposition uses 
the output from liqPhaseComposition (including water fraction), and chemGroups are taken from
the componentData.

freezingPointDepression
***********************
Inputs:
    *   compounds
    *   T
    *   fug_vap
    *   compoundData
    *   P
    *   chemGroups
  
Returns the freezing point depression of the aqueous phase from the guest gases *as a negative
number* to be added to 273.15 to obtain the freezing point.

frac(T, vaporFugacities, compoundData, structure, compounds, Ac, Bc, Dc)
************************************************************************

deltaHydratePotential(T, structure, vaporFugacities, compoundData, compounds, Ac, Bc, Dc)
*****************************************************************************************

waterFugacity(T, P, phase, fug_vap, compounds, compoundData)
************************************************************

hydrateFugacity(T, P, PvapConsts, structure, fug_vap, compounds, compoundData, Ac, Bc, Dc)
******************************************************************************************

class KlaudaSandler2003
************************