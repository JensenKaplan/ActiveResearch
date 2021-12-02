# WELCOME
This GitHub is a showcase of code that I've written and used during the ongoing course of my PhD program at the Georgia Institute of Technology.
As a condensed matter experimentalist I explore quantum magnetism in frusatrated lattice systems, with Lanthanides as the central ion. (I study crystals).

## What's in this Git?
Here you'll find Jupyter notebooks I've created that I regularly use to carry out my research. They will be complete and ready to run as-is, carrying out statistical analysis on the associated data files.

## What types of measurements do I use?
There are some common measurements we use to characterize our crystalline systems:

#### QuantumDesign PPMS Dynacool with Vibrating Sample Magnetometer
1. Magnetization (M vs H) measured at different temperatures.
2. Susceptibility (M vs T) Field Cooling / Zero Field Cooling (with an applied .01T magnetic field)
#### OakRidge National Lab (ORNL) Neutron Beamlines
1. Inelastic Neutron Scattering

#### What do I do with the measured data?
1. Thermodynamic (M vs T) data as measured by the PPMS is analyzed to find the Curie and Weiss constants. The Weiss Constant give us insight into the ferromagnetic / antiferromagnetic nature exhibited by our systems. The Curie Constant is used to calculate the effective magnetic moment of our system.
2. Thermodynamic (M vs H) data as measured by the PPMS is analyzed to find saturation points of my compounds.
3. Neutron data is analyzed to find low lying (in energy) magnetic modes in our systems. Analysis is carried out to determine the Crystal Electric Field Hamiltonian.

# Important Details of My Code (kwargs)
I've been implementing changes that allows all of my scripts to be more general. As well as adjusting the scripts in this repo I have conventionalized how I handle data in my DropBox. There are a few important kwargs that should be defined in every script so that the proper JensenTools functions are called.

1. "LS_on" a boolean. If True then all calculations done in LS basis; if False then calculations performed in J basis.
1. "comp" a string. This allows the script to find the proper data in my Dropbox.
2. "ion" a string. This is needed if working the J basis. Not necessary for LS basis. This is due to how PCF handles LS vs J calculations.
3. "KMeans" a boolean. If True then Kmeans sorting will be used to handle the eigenvalues; if False then all eigenvalues as produced by PCF will be stored.
4. "LSValue" a float. The LS coupling strength of the ion in meV.
5. "numlevels" an int. The number of excited states expected for the system. Works in conjunction with "KMeans". If KMeans is True then the number of levels to be calculated/stored is given by this variable.
6. "who" a string. Etiher 'Arun' or 'PPMS' this determines how M vs H / M vs T data is read and stored as well as how filenames are read. This is due to a difference in convention between Arun's filenaming and measurement files / PPMS filenaming and measurement files.
7. "molweight" is a float. It is the molarweight of the given compound. This is retrieved from a hardcoded dictionary imported from JensenTools.py
8. "per" is a string. Either "mol" or "spin" for normalizing measured thermodynamic data. We tend to look at data normalized per mol; PCF calculations are done per spin, therefore for comparing data you would use "spin".

Below is a codeblock example of what should be near the top of each script.
```
comp = 'Sr2PrO4'
ion = 'Ce3+'
LS_on = True
Kmeans = True
LSValue = 100
numlevels = 4
who = 'Arun'
molweight = molweight[comp]
per = 'spin'
```