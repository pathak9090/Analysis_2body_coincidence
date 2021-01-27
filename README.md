# Double_sided_VMI_coincidence_analysis
## Coincidence Analysis For a Double-Sided Velocity Map Imaging (VMI) Setup

### Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)

### General info
This project analyses and calculates the momentum and energies based on positions and time-of-flight for individual particles hitting the position sensitive detector. The input data file is a large binary file with positions and time for several hundred million particles. This analysis works for two-particle and three-particle coincidence, it can also be easily extended to multi-body coincidences.
	
### Technologies
Project is created with:
* fast-histogram version: 0.9 (Offers 25 times more speed compared to numpy.histogram2d)
* h5py version: 3.1.0
* SciPy version: 1.6.0
* NumPy version: 1.19.0
* tkinter version: 3.9.1

### Setup
To run this project, install it locally. The input data files are generated using COBOLD PC 2008 software.

