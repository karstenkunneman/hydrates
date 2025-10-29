import matplotlib.pyplot as plt
from src import core
from src import model

#Define simulation conditions
components = [1,2]
moleFractions = [0.5, 0.5]
Tlist = [265, 270, 275, 280, 285]
saltConcs = [5, 0, 0, 0, 0, 0, 0, 0, 0] #5 wt% NaCl
inhibitorConcs = [0, 5, 0, 0, 0, 0] #5 wt% Ethanol

#Pull data from file
componentData = core.getComponentData(components)
interactionParameters = core.getInteractionParameters(components)

#Usually not called directly, but used here do demonstrate PRSV function
guessPlist = [core.guessPressure(componentData, moleFractions, T) for T in Tlist]
    
#Main simulation loop
simResults = [None for i in range(len(Tlist))]
for i in range(len(simResults)):
    simResults[i] = model.KlaudaSandler2003(components, moleFractions, "T", Tlist[i])

#Pull data from simulation objects
eqPressures = [result.pressure for result in simResults]
eqStructures = [result.eqStructure for result in simResults]
eqFracs = [result.eqFrac for result in simResults]
eqPhases = [result.eqPhase for result in simResults]
hydrationNumbers = [result.hydrationNumber for result in simResults]
hydrateDensities = [result.hydrateDensity for result in simResults]
freezingPoints = [result.freezingPoint for result in simResults]

#Fugacity (in Pa) of gases at equilibrium
gas_fugs = [0. for i in range(len(Tlist))]
for i in range(len(gas_fugs)):
    gas_fugs[i] = core.PRSV(componentData, moleFractions, Tlist[i], eqPressures[i], interactionParameters)

#Determine inhibition of hydrate formation by salts and inhibitors
betaGas = core.betaGas(Tlist, eqPressures)
inhibitedTemps = [core.HuLeeSum(T, saltConcs, inhibitorConcs, betaGas, freezingPoint) for T, freezingPoint in zip(Tlist, freezingPoints)]

#Plot results
plt.plot(Tlist, eqPressures)
plt.plot(inhibitedTemps, eqPressures)
plt.title("Equilibrium Predictions")
plt.xlabel("Temperature (K)")
plt.ylabel("Pressure (Pa)")
plt.show

#Data Output
print("Guessed Pressures (Pa): " + str([round(num) for num in guessPlist]))
print("Equil. Pressures (Pa): " + str([round(num) for num in eqPressures]))
print("Equil. Gas Fugacities (Pa): " + str([round(sum(num[2])) for num in gas_fugs]))
print("Equil. Hydrate Structures: " + str(eqStructures))
print("Equil. Hydrate Fracs: " + str([[[round(value, 3) for value in row] for row in array] for array in eqFracs])) 
print("Hydration Numbers: " + str([round(num, 2) for num in hydrationNumbers]))
print("Hydrate Densities (kg/m^3): " + str([round(num[0], 2) for num in hydrateDensities]))
print("Freezing Points (K): " + str([round(num, 2) for num in freezingPoints]))
