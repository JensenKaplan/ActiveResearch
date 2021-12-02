## The Cubic Crystal Field Hamiltonian
Cubic symmetry allows for a simple Crystal Field Hamiltonian which consists of only 2 independent Stevens Coefficients:

H<sub>CUB</sub> = B<sub>4</sub>(O<sub>4</sub><sup>0</sup> + 5B<sub>4</sub><sup>4</sup>) + B<sub>6</sub>(O<sub>6</sub><sup>0</sup> - 21B<sub>6</sub><sup>4</sup>)

I parameterize this hamiltonian by x = B<sub>6</sub><sup>0</sup>/B<sub>4</sub><sup>0</sup> and bpf, which is an overall pre factor applied to all the coefficients resulting in the following form:

H<sub>CUB</sub> = bpf(O<sub>4</sub><sup>0</sup> + 5B<sub>4</sub><sup>4</sup>) + x\*bpf(O<sub>6</sub><sup>0</sup> - 21B<sub>6</sub><sup>4</sup>)

This was inspired by [Boothroyd's paper on PrO<sub>2</sub>](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.2082).

## Cubic Grid Creation
Create_Grid.py is used to create .mat files that contain the eigenvalues calculated by PCF, enforcing cubic constraints, for all (x,bpf).
The user defines the range for both x and bpf as well as the number of steps for each. It's recommended to keep the dimensions equal for more subtle reasons. 

I use x = \[-1,1\] and bpf = \[-1,1\]. 

I parallelized the process for eigenvalue calculation; however, I use K-Means sorting to keep track of the eigenvalues (crystal field levels), which is a bit intensive so this could take a while to run depending on size of x/bpf arrays. The idea is to only need to run this once to create the data, then the same files can be read to find good starting points for coefficients for a given compound.

Unique, non-zero eigenvalues for all (x,bpf) combinations are stored. The number of unique levels is given by the user, and this should be determined before hand. For example in the J basis for Pr<sup>4+</sup> 6 eigenvalues are produced. With cubic constraints, either a groundstate quartet with a doublet excited state or a doublet groundstate with a quartet excited state are calculated. I define the number of levels to be 1. PCF is used to calculate the eigenvalues and K-means is implemented to track the energy levels.

## Cubic Grid Search
Grid_Search.py is used to find the starting point Stevens Coefficients to be used for fitting.

1. The user declares the crystal field energy levels as measured by inelastic neutron scattering.
2. The user declares an allowed tolerance between measured levels and calculated levels. This tolerance may need to be bigger if the grid of (x,bpf) is not so fine. 
3. Eigenvalues and measured levels are matched within that tolerance level. If "compatible coordinates" are found that produce matching eigenvalues they are printed and these can be used as starting points for the given compound.