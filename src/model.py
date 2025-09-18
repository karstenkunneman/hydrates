import core
import math
import scipy
import numpy
from thermo.unifac import UNIFAC, PSRKSG, PSRKIP
import ast
import pandas as pd

R = 8.31446261815324 #m^3 Pa/mol K
boltzmannConstant = 1.380649E-23 #m^2 kg/s^2 K
N_A = 6.02214076E23 #mol^-1
errorMargin = 1E-5 #Converges to within 1 Pa precision

langParameters = pd.read_csv('KlaudaSandlerLangParams.csv', encoding='utf-8-sig').to_numpy()
PvapConsts = pd.read_csv('KlaudaSandlerP_vap.csv', encoding='utf-8-sig').to_numpy()

hydrateCellProperties = pd.read_csv('hydrateCellProperties.csv', encoding='utf-8-sig').to_numpy()
def getHydrateCellProperties(structure):
    mask = hydrateCellProperties[:, 0] == structure
    cellProperties = hydrateCellProperties[mask]
    return cellProperties

def henrysLawConst(componentData, T):
    constants = componentData[10:14]
    H_i = 101325*math.exp(-1*(constants[0]/1.987 + constants[1]/T/1.987 + constants[2]*math.log(T)/1.987 + constants[3]*T/1.987))
    return H_i 

def Z(compoundData, T, P):
    noComponents = len(compoundData)
    Z = numpy.zeros(noComponents)
    waterData = numpy.array([0, "H2O", 647.3, 22.09, 0.3438, -0.06635, {16: 1}, None, None, 18.02, None, None, None, None, 0, 0, 0,	0])
    for i in range(noComponents):
        localCompoundData = compoundData[i]
        localCompoundData = numpy.column_stack((localCompoundData, waterData))
        localCompoundData = localCompoundData.transpose()
        Z[i] = core.PRSV(localCompoundData, [0.00001, 0.99999], 273.15, P, [[0,0],[0,0]])[3]
    return Z

def liqPhaseComposition(compounds, T, fug_vap, compoundData, P):
    noComponents = len(compounds)
    x = numpy.zeros(noComponents+1)
    for i in range(noComponents):
        Z_inf = Z(compoundData, T, P)
        x[i+1] = fug_vap[i]/(henrysLawConst(compoundData[i], T)*math.exp(Z_inf[i]))
    x[0] = 1-sum(x) #Water composition
    return x

def activityCoeff(T, phaseComposition, chemGroups):
    chemGroupList = [{16: 1}]
    for i in range(len(chemGroups)):
        chemGroupList.append(ast.literal_eval(chemGroups[i]))

    GE = UNIFAC.from_subgroups(T, phaseComposition, chemGroupList, interaction_data=PSRKIP, subgroups=PSRKSG).gammas()[0]
    return GE

def freezingPointDepression(compounds, T, fug_vap, compoundData, P, chemGroups):
    phaseComposition = liqPhaseComposition(compounds, T, fug_vap, compoundData, P)
    deltadT = R*(273.15)**2/6011*math.log(liqPhaseComposition(compounds, T, fug_vap, compoundData, P)[0]*activityCoeff(T, phaseComposition, chemGroups))
    return deltadT

def Lang_Const(T, Ac, Bc, Dc):    
    Cml = math.exp(Ac + Bc/T + Dc/T/T)
    
    return Cml

def Lang_GG_Const(T, compoundData, fracs, structure):
    try:
        I = compoundData[:,7]#eV
        alpha = compoundData[:,8]#angstrom^3
    except:
        I = [compoundData[0][7]]
        alpha = [compoundData[0][8]]
        
    noComponents = len(I)    
        
    Cgg = [[1,1] for i in range(noComponents)]
    
    C6 = numpy.zeros((noComponents,noComponents))
    C8 = numpy.zeros((noComponents,noComponents))
    C10 = numpy.zeros((noComponents,noComponents))
            
    #Dispersion Coefficients for each guest-guest interaction
    for i in range(noComponents):
        for j in range(noComponents):
            C6[i][j] = (3/2)*alpha[i]*alpha[j]*I[i]*I[j]/(I[i]+I[j])*23.05/.001987
            C8[i][j] = 496778.3824*alpha[i]*alpha[j]*(I[i]/(2*I[i]+I[j])+I[j]/(2*I[j]+I[i]))
            C10[i][j] = 13260598.42*alpha[i]*alpha[j]/(I[i]+I[j])
    
    wrgg = [[[[0,0],[0,0]] for i in range(noComponents)] for i in range(noComponents)]
    
    #Unsure exactly what these are relative to r, just dont touch it if you want to stay sane.
    if structure == "I":
        a_lc = 12.03
        for i in range(noComponents):
            for j in range(noComponents):
                r0 = a_lc*0.86602
                wrgg[i][j][0][0] = -1*C6[i][j]*12.25367/(r0**6)-C8[i][j]*10.3552/(r0**8)-C10[i][j]*9.5644/(r0**10)
                r1 = a_lc*0.55901
                wrgg[i][j][0][1] = -1*C6[i][j]*13.41525/(r1**6)-C8[i][j]*12.38994/(r1**8)-C10[i][j]*12.12665/(r1**10)
                
                r0 = a_lc*0.55901
                wrgg[i][j][1][0] = -1*C6[i][j]*4.47175/(r0**6)-C8[i][j]*4.12998/(r0**8)-C10[i][j]*4.02482/(r0**10)
                r1 = a_lc*0.5
                wrgg[i][j][1][1] = -1*C6[i][j]*5.14048/(r1**6)-C8[i][j]*3.74916/(r1**8)-C10[i][j]*3.09581/(r1**10)
    elif structure == "II":
        a_lc = 17.31
        for i in range(noComponents):
            for j in range(noComponents):
                r0 = a_lc*0.35355
                wrgg[i][j][0][0] = -1*C6[i][j]*6.92768/(r0**6)-C8[i][j]*6.23392/(r0**8)-C10[i][j]*6.06724/(r0**10)
                r1 = a_lc*0.41457
                wrgg[i][j][0][1] = -1*C6[i][j]*6.91143/(r1**6)-C8[i][j]*6.28271/(r1**8)-C10[i][j]*6.10214/(r1**10)
                
                r0 = a_lc*0.41457
                wrgg[i][j][1][0] = -1*C6[i][j]*13.82287/(r0**6)-C8[i][j]*12.56542/(r0**8)-C10[i][j]*12.20428/(r0**10)
                r1 = a_lc*0.43301
                wrgg[i][j][1][1] = -1*C6[i][j]*5.11677/(r1**6)-C8[i][j]*4.33181/(r1**8)-C10[i][j]*4.11102/(r1**10)
    
    #Obtain the interaction energy of guestz`
    for i in range(noComponents):
        for j in range(noComponents):
            for k in range(2):
                Cgg[i][k] *= math.exp(-1*(wrgg[i][j][k][0]*fracs[0][j])/T)*math.exp(-1*(wrgg[i][j][k][1]*fracs[1][j])/T)
    
    return Cgg

#Calculates the fractional occupancy of small and large shells by each component
def frac(T, vaporFugacities, compoundData, structure, compounds, Ac, Bc, Dc):
    noComponents = len(vaporFugacities)

    guessFractions = numpy.zeros((2, noComponents))
    oldGuessFractions = [[1 for i in range(noComponents)],[1 for i in range(noComponents)]]
    Cgg = [[1.5, 1.5] for i in range(noComponents)]
    langConsts = [[0, 0] for i in range(noComponents)]
    fracDiff = [[1 for i in range(noComponents)],[1 for i in range(noComponents)]]

    if structure == "I":
        structureIndex = 0
    else:
        structureIndex = 1
    
    #Iterate over shell occupancy until the average fractional occupancy of all shells does not change by a factor greater than the error margin
    while abs(sum(fracDiff[0])/noComponents+sum(fracDiff[1])/noComponents)/2 >= errorMargin: #average fractional occupancy difference of all shells
        fracDiff = numpy.zeros((2, noComponents))
        for i in range(2):
            denominator = 0
            
            for j in range(noComponents):
                langConsts[j][i] = Lang_Const(T, Ac[structureIndex][i][j], Bc[structureIndex][i][j], Dc[structureIndex][i][j])*math.sqrt(1-fracDiff[i][j])
                denominator += Cgg[j][i]*langConsts[j][i]*vaporFugacities[j]        
            
            #Determine the difference in fractional occupancy between the previous guess fraction and the new guess fraction for each component
            for j in range(noComponents):
                guessFractions[i][j] = Cgg[j][i]*langConsts[j][i]*vaporFugacities[j]/(1 + denominator)
                fracDiff[i][j] = abs(guessFractions[i][j]/oldGuessFractions[i][j]-1)/noComponents
                oldGuessFractions[i][j]=guessFractions[i][j]
            
        #Guess a new guest-weighted Langmuir constant with new guess Fractions
        Cgg = Lang_GG_Const(T, compoundData, guessFractions, structure)

    for i in range(2):
        for j in range(noComponents):
            if guessFractions[i][j] < 1E-10:
                guessFractions[i][j] = 0
        
    return guessFractions

def deltaHydratePotential(T, structure, vaporFugacities, compoundData, compounds, Ac, Bc, Dc):
    cellProperties = getHydrateCellProperties(structure) 
    fractions = 0
    Deltamu_H_w = 0

    fractions = frac(T, vaporFugacities, compoundData, structure, compounds, Ac, Bc, Dc)
    
    Deltamu_H_w += cellProperties[0][10]*math.log(1-sum(fractions[0]))
    Deltamu_H_w += cellProperties[1][10]*math.log(1-sum(fractions[1]))
    return Deltamu_H_w, fractions

def waterFugacity(T, P, phase, fug_vap, compounds, compoundData):
    if phase == "ice":
        Vm_water = 1.912E-5 + T*8.387E-10 + (T**2)*4.016E-12
        Psat_water = math.exp(4.6056*math.log(T)-5501.1243/T+2.9446-T*8.1431E-3)
        f_w = Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))
    elif phase == "liquid":
        chemGroups = compoundData[:,6]
        Vm_water = math.exp(-10.921 + 2.5E-4*(T-273.15) - 3.532E-4*(P/1E6-0.101325) + 1.559E-7*((P/1E6-.101325)**2))
        Psat_water = math.exp(4.1539*math.log(T)-5500.9332/T+7.6537-16.1277E-3*T)
        phaseComposition = liqPhaseComposition(compounds, T, fug_vap, compoundData, P)
        f_w = phaseComposition[0]*activityCoeff(T, phaseComposition, chemGroups)*Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))
    return f_w

#Equation A2 in Sandler (2000)
def hydrateFugacity(T, P, PvapConsts, structure, fug_vap, compounds, compoundData, Ac, Bc, Dc):
    noComponents = len(compounds)
    cellProperties = getHydrateCellProperties(structure)

    #Molar volume of water (m^3/mol)
    if structure == "I":
        Vm_water = (11.835+2.217E-5*T+2.242E-6*T**2)**3*(1E-30*N_A/46)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2
    elif structure == "II":
        Vm_water = (17.13+2.249E-4*T+2.013E-6*T**2-1.009E-9*T**3)**3*(1E-30*N_A/136)-8.006E-9*P/1E6+5.448E-12*(P/1E6)**2
    
    dH = deltaHydratePotential(T, structure, fug_vap, compoundData, compounds, Ac, Bc, Dc)
    frac = dH[1]
    
    A = 0
    B = 0
    D = 0
    
    denominator = 0
    for i in range(noComponents):
        denominator += frac[0][i]*cellProperties[0][10]
        denominator += frac[1][i]*cellProperties[1][10]
            
    Z = numpy.zeros(noComponents)
    for i in range(noComponents):
        for j in range(2):
            if denominator != 0:
                Z[i] += (frac[j][i]*cellProperties[j][10])/denominator
            else:
                Z[i] = 0
            
    for i in range(noComponents):
        A += PvapConsts[i,3]*Z[i]
        B += PvapConsts[i,4]*Z[i]
        D += PvapConsts[i,6]*Z[i]
    
    Psat_water = math.exp(A*math.log(T)+B/T+2.7789+D*T)
    
    f_h = Psat_water*math.exp(Vm_water*(P-Psat_water)/(R*T))*math.exp(dH[0])
    return f_h,frac

#J. B. Klauda, S. I. Sandler, Phase behavior of clathrate hydrates:
#a model for single and multiple gas component hydrates, Chem-
#ical Engineering Science (2003). doi:https://doi.org/10.1016/
#S0009-2509(02)00435-9.
class KlaudaSandler2003:
    def __init__(self, componentList, moleFractions, definedVariable = "T", temperature = None, pressure = None):
        self.componentData = core.getComponentData(componentList)
        self.interactionParameters = core.getInteractionParameters(componentList)
        
        self.componentList = componentList
        self.moleFractions = moleFractions
        self.definedVariable = definedVariable
        
        if temperature != None and definedVariable == "T":
            self.temperature = temperature
            if pressure != None:
                self.pressure = pressure
            else:
                self.pressure = core.guessPressure(self.componentData, self.moleFractions, self.temperature)
        
        elif pressure != None and definedVariable == "P":
            self.pressure = pressure
            if temperature != None:
                self.temperature = temperature
            else:
                self.temperature = core.guessTemperature(self.componentData, self.moleFractions, self.pressure)
            
        else:
            raise Exception("Chosen defined variable is not defined")
    
        self.simResults = self.simulation()
    
    def simulation(self):
        noCompounds = len(self.componentList)
        
        freezingPoint = 273.15+freezingPointDepression(self.componentList, self.temperature, core.PRSV(self.componentData, self.moleFractions, 273.15, self.pressure, self.interactionParameters)[2], self.componentData, self.pressure, self.componentData[:,6])
        if self.temperature > freezingPoint:
            waterPhase = "liquid"
        else:
            waterPhase = "ice"
        
        def f_defT(P, T):
            vaporFugacities = core.PRSV(self.componentData, self.moleFractions, abs(T), abs(P), self.interactionParameters)[2]
            f_w = waterFugacity(abs(T), abs(P), waterPhase, vaporFugacities, self.componentList, self.componentData)
            f_h = hydrateFugacity(abs(T), abs(P), localPvapConsts, structure, vaporFugacities, self.componentList, self.componentData, Ac, Bc, Dc)[0]
            return math.log(f_h/f_w)
        
        def f_defP(T, P):
            vaporFugacities = core.PRSV(self.componentData, self.moleFractions, abs(T), abs(P), self.interactionParameters)[2]
            f_w = waterFugacity(abs(T), abs(P), waterPhase, vaporFugacities, self.componentList, self.componentData)
            f_h = hydrateFugacity(abs(T), abs(P), localPvapConsts, structure, vaporFugacities, self.componentList, self.componentData, Ac, Bc, Dc)[0]
            return math.log(f_h/f_w)
          
        structure = "I"

        #Get Langmuir parameters for each compound (preloaded for performance)
        Ac = numpy.zeros((2,2,noCompounds))
        Bc = numpy.zeros((2,2,noCompounds))
        Dc = numpy.zeros((2,2,noCompounds))
        for i in range(2):
            for j in range(2):
                for k in range(noCompounds):
                    compound = self.componentList[k]
                    condition = (langParameters[:,0] == structure) & \
                    (langParameters[:,1] == j) & \
                    (langParameters[:,2] == compound)

                    try:
                        row_index = numpy.where(condition)[0][0]

                        Ac[i][j][k] = langParameters[row_index, 3]
                        Bc[i][j][k] = langParameters[row_index, 4]
                        Dc[i][j][k] = langParameters[row_index, 5]
                    except:
                        Ac[i][j][k] = -100
                        Bc[i][j][k] = 0
                        Dc[i][j][k] = 0
        
        mask = [0 for i in range(len(PvapConsts[:,0]))]
        for i in range(len(PvapConsts[:,0])):
            if sum(numpy.isin(element = PvapConsts[:,1],test_elements=self.componentList)) > 1:
                structureMask = numpy.isin(element = PvapConsts[:,0],test_elements="I")
                componentMask = numpy.isin(element = PvapConsts[:,1],test_elements=self.componentList)
                mask[i] = numpy.logical_and(structureMask, componentMask)[i]
            else:
                mask[i] = True
            
        localPvapConsts = PvapConsts[mask]
        
        if self.definedVariable == "T":
            pGuess = self.pressure
            try:
                if "I" in localPvapConsts[:,0]:
                    SIEqPressure = abs(scipy.optimize.fsolve(f_defT,pGuess,xtol=errorMargin,args=self.temperature)[0])
                    SIEqFrac = hydrateFugacity(self.temperature, SIEqPressure, localPvapConsts, structure, core.PRSV(self.componentData, self.moleFractions, self.temperature, pGuess, self.interactionParameters)[2], self.componentList, self.componentData, Ac, Bc, Dc)[1]
                else:
                    raise
            except:
                SIEqPressure = math.inf
                SIEqFrac = numpy.zeros((2,len(self.moleFractions)))
        elif self.definedVariable == "P":
            tGuess = self.temperature
            try:
                if "I" in localPvapConsts[:,0]:
                    SIEqTemperature = abs(scipy.optimize.fsolve(f_defP,tGuess,xtol=errorMargin,args=self.pressure)[0])
                    SIEqFrac = hydrateFugacity(SIEqTemperature, self.pressure, localPvapConsts, structure, core.PRSV(self.componentData, self.moleFractions, tGuess, self.pressure, self.interactionParameters)[2], self.componentList, self.componentData, Ac, Bc, Dc)[1]
                else:
                    raise
            except:
                SIEqTemperature = 0
                SIEqFrac = numpy.zeros((2,noCompounds))
        
        structure = "II"

        #Get Langmuir parameters for each compound (preloaded for performance)
        Ac = numpy.zeros((2,2,noCompounds))
        Bc = numpy.zeros((2,2,noCompounds))
        Dc = numpy.zeros((2,2,noCompounds))
        for i in range(2):
            for j in range(2):
                for k in range(noCompounds):
                    compound = self.componentList[k]
                    condition = (langParameters[:,0] == structure) & \
                    (langParameters[:,1] == j) & \
                    (langParameters[:,2] == compound)

                    row_index = numpy.where(condition)[0][0]

                    Ac[i][j][k] = langParameters[row_index, 3]
                    Bc[i][j][k] = langParameters[row_index, 4]
                    Dc[i][j][k] = langParameters[row_index, 5]
        
        mask = [0 for i in range(len(PvapConsts[:,0]))]
        for i in range(len(PvapConsts[:,0])):
            if sum(numpy.isin(element = PvapConsts[:,1],test_elements=self.componentList)) > 1:
                structureMask = numpy.isin(element = PvapConsts[:,0],test_elements="II")
                componentMask = numpy.isin(element = PvapConsts[:,1],test_elements=self.componentList)
                mask[i] = numpy.logical_and(structureMask, componentMask)[i]
            else:
                mask[i] = True
            
        localPvapConsts = PvapConsts[mask]

        if self.definedVariable == "T":
            pGuess = self.pressure
            try:
                if "II" in localPvapConsts[:,0]:
                    SIIEqPressure = abs(scipy.optimize.fsolve(f_defT,pGuess,xtol=errorMargin,args=self.temperature)[0])
                    SIIEqFrac = hydrateFugacity(self.temperature, SIEqPressure, localPvapConsts, structure, core.PRSV(self.componentData, self.moleFractions, self.temperature, pGuess, self.interactionParameters)[2], self.componentList, self.componentData, Ac, Bc, Dc)[1]
                else:
                    raise
            except:
                SIIEqPressure = math.inf
                SIIEqFrac = numpy.zeros((2,len(self.moleFractions)))
                
            if SIIEqPressure >= SIEqPressure:
                self.eqStructure = "I"
                self.EqFrac = SIEqFrac
            else:
                self.eqStructure = "II"
                self.EqFrac = SIIEqFrac
            self.pressure = min(SIEqPressure, SIIEqPressure)
            
        elif self.definedVariable == "P":
            tGuess = self.temperature
            try:
                if "II" in localPvapConsts[:,0]:
                    SIIEqTemperature = abs(scipy.optimize.fsolve(f_defP,tGuess,xtol=errorMargin,args=self.pressure)[0])
                    SIIEqFrac = hydrateFugacity(SIEqTemperature, self.pressure, localPvapConsts, structure, core.PRSV(self.componentData, self.moleFractions, tGuess, self.pressure, self.interactionParameters)[2], self.componentList, self.componentData, Ac, Bc, Dc)[1]
                else:
                    raise
            except:
                SIIEqTemperature = 0
                SIIEqFrac = numpy.zeros((2,noCompounds))
        
            if SIIEqTemperature <= SIEqTemperature and SIEqTemperature != math.inf:
                self.eqStructure = "I"
                self.EqFrac = SIEqFrac
            else:
                self.eqStructure = "II"
                self.EqFrac = SIIEqFrac
        
            if SIEqTemperature != math.inf and SIIEqTemperature != math.inf:
                self.temperature = max(SIEqTemperature, SIIEqTemperature)
            elif SIEqTemperature == math.inf and SIIEqTemperature != math.inf:
                self.temperature = SIIEqTemperature
            else:
                self.temperature = SIEqTemperature

        if waterPhase == "ice":
            self.equilPhase = "I-H-V"
        else:
            self.equilPhase = "L-H-V"
            
        self.hydrationNumber = core.hydrationNumber(self.eqStructure, self.EqFrac)
        self.hydrateDensity = core.hydrateDensity(self.eqStructure, self.EqFrac, self.componentData, self.moleFractions, self.temperature, self.pressure)
        self.freezingPoint = freezingPoint

        return 0