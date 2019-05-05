import random
import pylab
import sys

class NoChildException(Exception):
    """
    NoChildException is raised by the reproduce() method in the SimpleVirus
    and ResistantVirus classes to indicate that a virus particle does not
    reproduce.
    """
    
class SimpleVirus(object):

    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    def __init__(self, maxBirthProb, clearProb):
        """
        Initialize a SimpleVirus instance, saves all parameters as attributes
        of the instance.        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        clearProb: Maximum clearance probability (a float between 0-1).
        """

        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb

    def getMaxBirthProb(self):
        """
        Returns the max birth probability.
        """
        return self.maxBirthProb

    def getClearProb(self):
        """
        Returns the clear probability.
        """
        return self.clearProb

    def doesClear(self):
        """ Stochastically determines whether this virus particle is cleared from the
        patient's body at a time step. 
        returns: True with probability self.getClearProb and otherwise returns
        False.
        """

        if self.clearProb >= random.random():
            return True
        else:
            return False

    
    def reproduce(self, popDensity):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient and
        TreatedPatient classes. The virus particle reproduces with probability
        self.maxBirthProb * (1 - popDensity).
        
        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring SimpleVirus (which has the same
        maxBirthProb and clearProb values as its parent).         

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population.         
        
        returns: a new instance of the SimpleVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.               
        """

        if self.maxBirthProb * (1 - popDensity) >= random.random():
            return SimpleVirus(self.maxBirthProb,self.clearProb)
        else:
            raise NoChildException()



class Patient(object):
    """
    Representation of a simplified patient. The patient does not take any drugs
    and his/her virus populations have no drug resistance.
    """    

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes.

        viruses: the list representing the virus population (a list of
        SimpleVirus instances)

        maxPop: the maximum virus population for this patient (an integer)
        """
        self.viruses = viruses
        self.maxPop = maxPop
        

    def getViruses(self):
        """
        Returns the viruses in this Patient.
        """
        return self.viruses


    def getMaxPop(self):
        """
        Returns the max population.
        """
        return self.maxPop


    def getTotalPop(self):
        """
        Gets the size of the current total virus population. 
        returns: The total virus population (an integer)
        """

        return len(self.viruses)       


    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute the following steps in this order:
        
        - Determine whether each virus particle survives and updates the list
        of virus particles accordingly.   
        
        - The current population density is calculated. This population density
          value is used until the next call to update() 
        
        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.                    

        returns: The total virus population at the end of the update (an
        integer)
        """
        self.survivorViruses = []
        self.newViruses = []
        
        for i in self.viruses:
               if not i.doesClear():
                  self.survivorViruses.append(i)
        
        self.popDensity = len(self.survivorViruses) / self.maxPop
        for i in self.survivorViruses:
                try:
                    pot = i.reproduce(self.popDensity)
                    if pot:
                        self.newViruses.append(i)
                except NoChildException:
                    pass
        
        self.viruses = self.survivorViruses + self.newViruses
        return len(self.viruses)
        


def simulationWithoutDrug(numViruses, maxPop, maxBirthProb, clearProb,
                          numTrials):
    """
    Run the simulation and plot the graph (no drugs are used,
    viruses do not have any drug resistance).    
    For each of numTrials trial, instantiates a patient, runs a simulation
    for 300 timesteps, and plots the average virus population size as a
    function of time.

    numViruses: number of SimpleVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: Maximum clearance probability (a float between 0-1)
    numTrials: number of simulation runs to execute (an integer)
    """
    yValues = []
    listViruses = []
    for i in range(numViruses):
       listViruses.append(SimpleVirus(maxBirthProb,clearProb))
       
    for i in range(numTrials):
        initial = Patient(listViruses, maxPop)
        for i in range(300):
            initial.update()
            if len(yValues) == i:
                 yValues.append(initial.getTotalPop())
            else:
                yValues[i] = yValues[i] + initial.getTotalPop()
    for y in range(len(yValues)):
        yValues[y] = yValues[y] / numTrials

    pylab.plot(yValues, label = "SimpleVirus")
    pylab.title("SimpleVirus simulation")
    pylab.xlabel("Time Steps")
    pylab.ylabel("Average Virus Population")
    pylab.legend(loc = "best")
    pylab.show()



class ResistantVirus(SimpleVirus):
    """
    Representation of a virus which can have drug resistance.
    """   

    def __init__(self, maxBirthProb, clearProb, resistances, mutProb):
        """
        Initialize a ResistantVirus instance, saves all parameters as attributes
        of the instance.

        maxBirthProb: Maximum reproduction probability (a float between 0-1)       

        clearProb: Maximum clearance probability (a float between 0-1).

        resistances: A dictionary of drug names (strings) mapping to the state
        of this virus particle's resistance (either True or False) to each drug.
        e.g. {'guttagonol':False, 'srinol':False}, means that this virus
        particle is resistant to neither guttagonol nor srinol.

        mutProb: Mutation probability for this virus particle (a float). This is
        the probability of the offspring acquiring or losing resistance to a drug.
        """

        SimpleVirus.__init__(self, maxBirthProb,clearProb)
        self.resistances = resistances
        self.mutProb = mutProb


    def getResistances(self):
        """
        Returns the resistances for this virus.
        """
        return self.resistances

    def getMutProb(self):
        """
        Returns the mutation probability for this virus.
        """
        return self.mutProb

    def isResistantTo(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This method
        is called by getResistPop() in TreatedPatient to determine how many virus
        particles have resistance to a drug.       

        drug: The drug (a string)

        returns: True if this virus instance is resistant to the drug, False
        otherwise.
        """
        if drug in self.resistances:
            if self.resistances[drug]:
                return True
            else:
                return False


    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the TreatedPatient class.

        A virus particle will only reproduce if it is resistant to ALL the drugs
        in the activeDrugs list. For example, if there are 2 drugs in the
        activeDrugs list, and the virus particle is resistant to 1 or no drugs,
        then it will NOT reproduce.

        Hence, if the virus is resistant to all drugs
        in activeDrugs, then the virus reproduces with probability:      

        self.maxBirthProb * (1 - popDensity).                       

        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring ResistantVirus (which has the same
        maxBirthProb and clearProb values as its parent). The offspring virus
        will have the same maxBirthProb, clearProb, and mutProb as the parent.

        For each drug resistance trait of the virus (i.e. each key of
        self.resistances), the offspring has probability 1-mutProb of
        inheriting that resistance trait from the parent, and probability
        mutProb of switching that resistance trait in the offspring.       

        For example, if a virus particle is resistant to guttagonol but not
        srinol, and self.mutProb is 0.1, then there is a 10% chance that
        that the offspring will lose resistance to guttagonol and a 90%
        chance that the offspring will be resistant to guttagonol.
        There is also a 10% chance that the offspring will gain resistance to
        srinol and a 90% chance that the offspring will not be resistant to
        srinol.

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population       

        activeDrugs: a list of the drug names acting on this virus particle
        (a list of strings).

        returns: a new instance of the ResistantVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.
        """
        resist = True
        newResistances = {}
        for drug in activeDrugs:
          if not len(self.getResistances()) == 0:
            if not self.resistances[drug] or not drug in self.resistances:
              resist = False

        if resist:            
            if not len(self.resistances) == 0:
              for drug in self.resistances:
                  if 1-self.mutProb >= random.random():
                     newResistances[drug] = self.resistances[drug]
                     
                  else:
                     newResistances[drug] = not self.resistances[drug]
                              
            if self.maxBirthProb * (1 - popDensity) >= random.random():
                return ResistantVirus(self.maxBirthProb,self.clearProb,newResistances,self.mutProb)
            else:
                raise NoChildException()
        else:
            raise NoChildException()
        
        
                

class TreatedPatient(Patient):
    """
    Representation of a patient. The patient is able to take drugs and his/her
    virus population can acquire resistance to the drugs he/she takes.
    """

    def __init__(self, viruses, maxPop, prescription):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered
        (which should initially include no drugs).              

        viruses: The list representing the virus population (a list of
        virus instances)

        maxPop: The  maximum virus population for this patient (an integer)
        
        prescription: The list of drugs for the treatment, but not implemented until simulation
        """

        Patient.__init__(self,viruses,maxPop)
        self.prescription = prescription
        self.Prescriptions = []


    def addPrescription(self, newDrugs):
        """
        Administer a drug to this patient. After a prescription is added, the
        drug acts on the virus population for all subsequent time steps. If the
        newDrug is already prescribed to this patient, the method has no effect.

        newDrug: The name of the drug to administer to the patient (a string).

        postcondition: The list of drugs being administered to a patient is updated
        """
        for d in newDrugs:
           if not newDrugs in self.Prescriptions:
              self.Prescriptions.append(newDrugs)


    def getPrescriptions(self):
        """
        Returns the drugs that are being administered to this patient.

        returns: The list of drug names (strings) being administered to this
        patient.
        """

        return self.Prescriptions


    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed in
        drugResist.       

        drugResist: Which drug resistances to include in the population (a list
        of strings - e.g. ['guttagonol'] or ['guttagonol', 'srinol'])

        returns: The population of viruses (an integer) with resistances to all
        drugs in the drugResist list.
        """
        self.resistantViruses = 0
        for i in self.viruses:
            resists = True
            for d in drugResist:
                 if not i.isResistantTo(d):
                     resists = False
            if resists:
                self.resistantViruses += 1
        return self.resistantViruses
        


    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:

        - Determine whether each virus particle survives and update the list of
          virus particles accordingly

        - The current population density is calculated. This population density
          value is used until the next call to update().

        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.
          The list of drugs being administered should be accounted for in the
          determination of whether each virus particle reproduces.

        returns: The total virus population at the end of the update (an
        integer)
        """

        self.survivorViruses = []
        self.newViruses = []

        for i in self.viruses:
               if not i.doesClear():
                  self.survivorViruses.append(i)
        

        self.popDensity = len(self.survivorViruses) / self.maxPop
        
        for i in self.survivorViruses:
            reproductiv = True
            for d in self.Prescriptions:
                if not len(i.getResistances()) == 0:
                   if not i.resistances[d]:
                      reproductiv = False
            try:
              if reproductiv or len(i.getResistances()) == 0:
                  try:
                      pot = i.reproduce(self.popDensity, self.Prescriptions)
                      if pot:
                          self.newViruses.append(pot)
                  except NoChildException:
                      pass
            except NoChildException:
                      pass
        
        self.viruses = self.survivorViruses + self.newViruses
        return len(self.viruses)
        

def simulationWithDrug(virus, patient, numTrials):
    """
    Runs simulations and plots graphs. Viruses have drug resistance and the 
    patient has a prescription.

    For each of numTrials trials, instantiates the patient, runs a simulation for
    150 timesteps, adds the prescription, and runs the simulation for an additional
    150 timesteps.  At the end plots the average virus population size
    (for both the total virus population and the drug-resistant virus
    population) as a function of time.

    numTrials: number of simulation runs to execute (an integer)
    
    """
    yValues = []
    yDrugValues = []
    listViruses = []
    for v in patient.viruses:

        listViruses.append(virus)     
            
    for t in range(numTrials):
      p = TreatedPatient(patient.viruses, patient.maxPop, [])

      for i in range(150):
        p.update()
        for adder in range(1):
            if len(yValues) == i:
               yValues.append(p.getTotalPop())
            else:
               yValues[i] = yValues[i] + p.getTotalPop()
        for adder in range(1):
            if len(yDrugValues) == i:
               yDrugValues.append(p.getResistPop(patient.getPrescriptions()))
            else:
               yDrugValues[i] = yDrugValues[i] + p.getResistPop(patient.getPrescriptions())
      
      for d in patient.prescription:
         p.addPrescription(d)
         
      for i in range(150):
        p.update()
        for adder in range(1):
            if len(yValues) == i+150:
               yValues.append(p.getTotalPop())
            else:
               yValues[i+150] = yValues[i+150] + p.getTotalPop()
        for adder in range(1):
            if len(yDrugValues) == i+150:
               yDrugValues.append(p.getResistPop(p.getPrescriptions()))
            else:
               yDrugValues[i+150] = yDrugValues[i+150] + p.getResistPop(p.getPrescriptions())
   
    for y in range(len(yValues)):
       yValues[y] = yValues[y] / numTrials
    for e in range(len(yDrugValues)):
       yDrugValues[e] = yDrugValues[e] / numTrials   

        
    pylab.plot(yValues, label = "ResistantVirus population")
    pylab.plot(yDrugValues, label = "prescription-resistant virus population")
    pylab.title("ResistantVirus simulation")
    pylab.xlabel("time step")
    pylab.ylabel("# viruses")
    pylab.legend(loc = "best")
    
    pylab.show()
    


#User Interface#

print("------------Welcome to the Virus Growth Simulation App!-------")
print("")
print("This program permits the user to create personalized Virus objects, as well as Patient objects, in order to test 100 simulations of virus growth with multiple variables and show the results in a graph. The graphs take some time to load.")
print("")
print("-Simple Virus: Without drug resistance, this virus reproduces by itself in a time step. Its properties are: Maximun Reproduction Probability, and Elimination Probability.")
print("")
print("-Patient: Representation of a simplified patient. The patient does not take any drugs and his/her virus populations have no drug resistance. It has a starting virus population and a max population variable.")
print("")
print("-Resistant Virus: A virus with a list of drug resistance(s) and mutation probability (which may or not alter its drug resistances after reproduction).")
print("")
print("-Treated Patient: A Patient object treated with one or more drugs. The prescription eliminates unresistant viruses, but Resistant viruses may acquire or lose drug resistances in the process.")
print("")

navigator = 4
while navigator != 1 or 2 or 3: 
    navigator = int(input("Enter 1 to create a simulation with Simple Viruses and a Patient, or enter 2 to experiment with Resistant Viruses and a Treated Patient. To exit, enter 3: "))
    print("")
    if navigator == 1:
       print("1. Creating our Simple Virus")
       print("")
       maxBirthProb = -1
       while maxBirthProb < 0 or maxBirthProb > 1:
          try:
             maxBirthProb = float(input("Enter the Maximun Reproduction Probability (float between 0 and 1): "))
             if maxBirthProb < 0 or maxBirthProb > 1:
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")

       clearProb = -1
       while clearProb < 0 or clearProb > 1:
          try:
             clearProb = float(input("Enter the Elimination Probability (float number between 0 and 1): "))
             if clearProb < 0 or clearProb > 1 :
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")        

       ourSimpleVirus = SimpleVirus(maxBirthProb, clearProb)
       print("--------------------------------------------------------------")
       print("Simple Virus creation was successful!")
       print("Maximun Reproduction Probability: ",str(maxBirthProb)," / ","Elimination Probability: ",str(clearProb))
       print("")
    
       print("2. Creating our Patient")
       print("")
       numVirus = 0
       listVirus = []
       while numVirus < 1:
          try:
             numVirus = int(input("Enter the starting Simple Virus population (int number greater than 0. Recommended to be lower than 400): "))      
             if numVirus < 1:
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")    
       for i in range(numVirus):
          listVirus.append(ourSimpleVirus)
       print("")
       maxPop = 0
       while maxPop < 1:
          try:
             maxPop = int(input("Enter the maximun Simple Virus population (int number greater than 0. Recommended to be lower than 400): "))      
             if maxPop < 1:
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")    

       ourPatient = Patient(listVirus,maxPop)
       print("--------------------------------------------------------------")
       print("Patient creation was successful")
       print("Initial Simple Virus population: ",str(numVirus)," / ","Maximun Simple Virus population: ",str(maxPop))
       print("")
       print("The program performs 100 simulations and takes the Average Virus Population of all of them at every Time Step. Increasing averages indicate virus growth in the patient over time, while decreasing averages indicate virus decrease over time.")
       simulationWithoutDrug(numVirus, maxPop, maxBirthProb, clearProb,100)
       print("")
       navigator = 4
    
    if navigator == 2:
       print("1. Creating our Resistant Virus")
       print("")
       maxBirthProb = -1
       while maxBirthProb < 0 or maxBirthProb > 1:
          try:
             maxBirthProb = float(input("Enter the Maximun Reproduction Probability (float between 0 and 1): "))
             if maxBirthProb < 0 or maxBirthProb > 1:
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")

       clearProb = -1
       while clearProb < 0 or clearProb > 1:
          try:
             clearProb = float(input("Enter the Elimination Probability (float number between 0 and 1): "))
             if clearProb < 0 or clearProb > 1 :
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")  
       
       resistances = {}
       resistanceAdded = "x"
       while resistanceAdded != "0":
           resistanceAdded = str(input("Write the name of the drug you want the virus to be initially resistant to (it can be any name. Some examples: amoxicillin, doxycycline, guttagonol, etc.). Enter 0 to stop adding drugs: "))
           if resistanceAdded != "0":
               resistances[resistanceAdded] = True
             
       mutProb = -1
       while mutProb < 0 or mutProb > 1:
          try:
             mutProb = float(input("Enter the Mutation Probability (float number between 0 and 1): "))
             if mutProb < 0 or mutProb > 1 :
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")                     
       ourResistantVirus = ResistantVirus(maxBirthProb, clearProb, resistances, mutProb)   
       print("--------------------------------------------------------------")
       print("Resistant Virus creation was successful!")
       print("Maximun Reproduction Probability: ",str(maxBirthProb)," / ","Elimination Probability: ",str(clearProb)," / ","Drug Resistances: ",str(resistances)," / ","Mutation Probability: ",str(mutProb))
       print("")       

       print("2. Creating our Treated Patient")
       print("")
       numVirus = 0
       listVirus = []
       while numVirus < 1:
          try:
             numVirus = int(input("Enter the starting Resistant Virus population (int number greater than 0. Recommended to be lower than 400): "))      
             if numVirus < 1:
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")    
             
       for i in range(numVirus):
          listVirus.append(ourResistantVirus)
          
       print("")
       maxPop = 0
       while maxPop < 1:
          try:
             maxPop = int(input("Enter the maximun Resistant Virus population (int number greater than 0. Recommended to be lower than 400): "))      
             if maxPop < 1:
                print("Please enter a valid number")
                print("")
          except ValueError:
             print("Please enter a valid number")
             print("")    

       drugAdded = "x"
       prescription = []
       while drugAdded != "0":
          drugAdded = str(input("Write the name of the drug you want to add to the Patient's prescription (it can be any name, and it kills viruses that don't have a resistance with the EXACT name. Some examples: amoxicillin, doxycycline, guttagonol, etc.). Enter 0 to stop adding drugs: "))
          if drugAdded != "0":
             prescription.append(drugAdded)
    
       ourTreatedPatient = TreatedPatient(listVirus,maxPop,prescription)
       
       for d in prescription:
           if d not in ourResistantVirus.getResistances():
              ourResistantVirus.resistances[d] = False
       
    
       print("--------------------------------------------------------------")
       print("Treated Patient creation was successful")
       print("Initial Resistant Virus population: ",str(numVirus)," / ","Maximun Simple Virus population: ",str(maxPop)," / ","Prescription: ",str(prescription))   
       print("")
       print("This is a graph of the average results of 100 simulations with two parts. First, from time steps 0 to 150, the patient starts without prescription. Then, since time step 150, the prescription is added to fight the viruses. The graph illustrates its efficiency. The blue line represents the entire virus population, and the orange line represents the virus population resistant to the prescription.")
       simulationWithDrug(ourResistantVirus, ourTreatedPatient, 100)
       print("")
    
       navigator = 4
       
    if navigator == 3:
       sys.exit("Thanks for using the Virus Growth Simulation App!")








#simulationWithDrug(100, 1000, 0.1, 0.05, {'guttagonol': False},0.005, 10)
#simulationWithDrug(1, 10, 1.0, 0.0, {}, 1.0, 5)
#simulationWithDrug(1, 20, 1.0, 0.0, {"guttagonol": True}, 1.0, 5)                #done
#simulationWithDrug(75, 100, .8, 0.1, {"guttagonol": True}, 0.8, 1)