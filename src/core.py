import math
import pandas as pd
import numpy as np
import scipy

R = 8.31446261815324 #m^3 Pa/mol K
boltzmannConstant = 1.380649E-23 #m^2 kg/s^2 K
N_A = 6.02214076E23 #mol^-1

fluidData = pd.read_csv('componentData.csv', encoding='utf-8-sig').to_numpy()
interactionParameters = pd.read_csv('interactionParameters.csv', encoding='utf-8-sig').to_numpy()
saltData = pd.read_csv('Salt Data.csv', encoding='utf-8-sig').to_numpy()
inhibitorData = pd.read_csv('Organic Inhibitor Data.csv', encoding='utf-8-sig').to_numpy()

#Obtain component data for given components present in the simulated system
def getComponentData(componentList):
    mask = np.isin(fluidData[:, 0], componentList)
    componentData = fluidData[mask]
    
    return componentData

def getInteractionParameters(componentList):
    compounds = np.nonzero(componentList)
    
    if len(compounds[0]) > 1:
        mixConstants = np.zeros((len(compounds[0]),len(compounds[0])))
        for i in range(len(compounds[0])):
            for j in range(len(compounds[0])):
                mixConstants[i][j] = interactionParameters[componentList[i]][componentList[j]+1]
    else:
        mixConstants = [[0]]
    
    return mixConstants

#Generate a pressure guess (in Pa) from input temperature and component mole fractions
def guessPressure(componentData, moleFractions, T):
    if T < 273.15:
        guessPressure = 1
        for i in range(len(moleFractions)):
            guessPressure *= moleFractions[i]*componentData[i][14]*math.exp(componentData[i][15]*T)
        
        return guessPressure**(1/len(moleFractions))
    else:
        guessPressure = 1
        for i in range(len(moleFractions)):
            guessPressure *= moleFractions[i]*componentData[i][16]*math.exp(componentData[i][17]*T)

        return guessPressure**(1/len(moleFractions))

#Generate a temperature guess (in K) from input pressure and component mole fractions
def guessTemperature(compounds, moleFractions, P):
    def pressureMatch(guessT, P, compounds, moleFractions):
        return P - guessPressure(compounds, moleFractions, guessT)
    
    guessTemp = scipy.optimize.fsolve(pressureMatch, [273.15], args = (P, compounds, moleFractions))[0]+5

    return guessTemp

#Peng-Robinson-Stryjek-Vera Equation of State
#D. Y. Peng, D. B. Robinson, A new two-constant equation of state, Industrial & Engineering Chemistry Fundamentals (1976). doi:https://doi.org/10.1021/i160057a011.
#R. Stryjek, J. H. Vera, Prsv: An improved peng—robinson equation of state for pure compounds and mixtures, The Canadian Journal of Chemical Engineering (1986). doi:https://doi.org/10.1002/cjce.5450640224.
def PRSV(compoundData, moleFractions, T, P, interactionParameters):
    noComponents = len(moleFractions)
    #Calculate "b" value for each pure substance and add to weighted and unweighted lists
    bmix = 0
    b = [0]*noComponents
    for i in range(noComponents):    
        Tc = compoundData[i][2]
        Pc = compoundData[i][3]*1E6
        b[i] = 0.07780*(R*Tc/Pc)
        bmix += b[i]*moleFractions[i]
        
    #Calculate the "a" values for all pure compounds
    a = [0]*noComponents
    for i in range(noComponents):
        Tc = compoundData[i][2]
        Pc = compoundData[i][3]*1E6
        omega = compoundData[i][4]
        kappa1 = compoundData[i][5]
        kappa = 0.378893 + 1.4897153*omega - 0.17131848*omega**2 + 0.0196554*omega**3 + kappa1*(1+math.sqrt(T/Tc))*(0.7-(T/Tc))
        alpha = (1 + kappa*(1-math.sqrt(T/Tc)))**2
        a[i] = 0.45724*((R**2)*(Tc**2)/Pc)*alpha
        
    #Finally, determine the partial "a" values based on mole ratios
    amix = 0
    xia = np.zeros(noComponents)
    for i in range(noComponents):
        for j in range(noComponents):
            amix += (math.sqrt(a[i]*a[j]))*(1-interactionParameters[i][j])*moleFractions[i]*moleFractions[j]
            xia[i] += (math.sqrt(a[i]*a[j]))*(1-interactionParameters[i][j])*moleFractions[j]
            
    A = float(amix*P/((R*T)**2))
    B = float((bmix*P)/(R*T))
    
    Z = np.roots([1, -1+B, A-3*(B**2)-2*B, -1*A*B+B**2+B**3])
    
    #Analyze roots
    noReal = 0
    realZ = []
    for i in range(len(Z)):
        if np.iscomplex(Z[i]) == False:
            noReal += 1
            realZ.append(Z[i])
    
    if noReal == 3:
        ZVmix = max(realZ).real
        ZLmix = min(realZ).real
    elif noReal == 1:
        ZVmix = realZ[0].real
        ZLmix = realZ[0].real
    
    #Calculate molar volumes for each compound
    VmVap = (ZVmix*R*T/P)
    VmLiq = (ZLmix*R*T/P)
    
    fugVap = np.zeros(noComponents)
    
    #Calculate fugacities for each of the compounds
    for i in range(noComponents):
        fugVap[i] = P*moleFractions[i]*math.exp((b[i]/bmix)*(ZVmix-1)-math.log(ZVmix-B)-(A/(2*math.sqrt(2)*B))*(2*xia[i]/amix-b[i]/bmix)*math.log((ZVmix+(1+math.sqrt(2))*B)/(ZVmix+(1-math.sqrt(2))*B)))
       
    return VmVap, VmLiq, fugVap, ZVmix

#Calculates the hydration number of a hydrate at equilibrium
def hydrationNumber(structure, occupancies):
    if sum(occupancies[0] + occupancies[1]) != 0:
        if structure == "I":
            hydrationNumber = 46/(sum(occupancies[0])*2+sum(occupancies[1])*6)
        else:
            hydrationNumber = 136/(sum(occupancies[0])*16+sum(occupancies[1])*8)
        
        return hydrationNumber
    else:
        return None
    
#Calculates the density of the equilibrium hydrate structure (in kg/m^3)
def hydrateDensity(structure, occupancies, compoundData, moleFractions, T, P):
    noCompounds = len(moleFractions)
    guestMass = 0
    
    molarmass = np.zeros(noCompounds)
    for i in range(noCompounds):
        molarmass[i] = compoundData[i][10]/1000
    
    if structure == "I":
        Vm_water = (11.835+2.217E-5*T+2.242E-6*T**2)**3*(1E-30*N_A/46)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2 #m3/mol
        waterMass = 18.02/1000/Vm_water #kg/m3
        for i in range(noCompounds):
            guestMass += (molarmass[i]*occupancies[0][i]*2 + molarmass[i]*occupancies[1][i]*6)/Vm_water/46
    elif structure == "II":
        Vm_water = (17.13+2.249E-4*T+2.013E-6*T**2-1.009E-9*T**3)**3*(1E-30*N_A/136)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2 #m3/mol
        waterMass = 18.02/1000/Vm_water #kg/m3
        for i in range(noCompounds):
            guestMass += (molarmass[i]*occupancies[0][i]*16 + molarmass[i]*occupancies[1][i]*8)/Vm_water/136
            
    return waterMass + guestMass, 1

#Hu-Lee-Sum Inhibition Model
#Yue Hu et al. (2018)
def HuLeeSum(T, saltConcs, inhibitorConcs, betaGas, freezingPoint):
    noSalts = len(saltConcs)
    noInhibitors = len(inhibitorConcs)
    if T >= freezingPoint:
        if sum(saltConcs) != 0:
            catMols = np.zeros(noSalts)
            anMols = np.zeros(noSalts)
            catMolFracs = np.zeros(noSalts)
            anMolFracs = np.zeros(noSalts)

            #Calculate mole fractions for the anions and cations in solution
            for i in range(noSalts):
                catMols[i] = saltConcs[i]/saltData[i][1]*saltData[i][4]
                anMols[i] = saltConcs[i]/saltData[i][1]*saltData[i][5]
            totalSaltMols = sum(catMols)+sum(anMols)
            waterMols = (100-sum(saltConcs))/18.015
            for i in range(noSalts):
                catMolFracs[i] = saltData[i][2]*catMols[i]/(waterMols+totalSaltMols)
                anMolFracs[i] = saltData[i][3]*anMols[i]/(waterMols+totalSaltMols)

            totalEffSaltMols = sum(catMolFracs) + sum(anMolFracs)
        else:
            waterMols = 100/18.015
            totalSaltMols = 0
            totalEffSaltMols = 0

        #Calculate mole fractions for the organic inhibitors in solution
        inhibitorMols = np.zeros(noInhibitors)
        inhibitorMolFracs = np.zeros(noInhibitors)
        for i in range(noInhibitors):
            inhibitorMols[i] = inhibitorConcs[i]/(1-0.01*inhibitorConcs[i])/inhibitorData[i][1]
        totalMols = totalSaltMols + sum(inhibitorMols)
        for i in range(noInhibitors):
            inhibitorMolFracs[i] = inhibitorMols[i]/(totalMols+waterMols)

        #Calculate the activity of the salt ions and organic inhibitors in solution
        lnawSalts = -1.06152*totalEffSaltMols+3.25726*totalEffSaltMols*totalEffSaltMols-37.2263*totalEffSaltMols*totalEffSaltMols*totalEffSaltMols
        lnawInhibitors = inhibitorMols = np.zeros(noInhibitors)
        for i in range(noInhibitors):
            lnawInhibitors[i] = inhibitorData[i][2]*inhibitorMolFracs[i]+inhibitorData[i][3]*inhibitorMolFracs[i]**2

        Tinhibited = T/(1-betaGas*(lnawSalts+sum(lnawInhibitors))*T)
    else:
        Tinhibited = None
    return Tinhibited

def checkMaxConc(inhibitorConcs):
    exceededInhibitors = ""
    for i in range(len(inhibitorConcs)):
        if inhibitorConcs[i] > inhibitorData[i][5]:
            exceededInhibitors += str(inhibitorData[i][0]) + " "
    return exceededInhibitors

def betaGas(temperatures, pressures):
    for i in range(len(pressures)):
        pressures[i]

    inverseTemp = []
    lnPressure = []
    for i in range(len(temperatures)):
        if temperatures[i] >= 273.15:
            inverseTemp.append(1/temperatures[i])
            lnPressure.append(math.log(pressures[i]))
    
    try:
        slope = np.polyfit(inverseTemp, lnPressure, 1)[0]
    
        betaGas = -10/slope

    except:
        betaGas = 0
        
    return betaGas