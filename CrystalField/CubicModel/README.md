## The Cubic Crystal Field Hamiltonian
Cubic symmetry allows for a simple Crystal Field Hamiltonian which consists of only 2 independent Stevens Coefficients:

H<sub>CUB</sub> = B<sub>4</sub>(O<sub>4</sub><sup>0</sup> + 5B<sub>4</sub><sup>4</sup>) + B<sub>6</sub>(O<sub>6</sub><sup>0</sup> - 21B<sub>6</sub><sup>4</sup>)

I parameterize this hamiltonian by x = B<sub>6</sub><sup>0</sup>/B<sub>4</sub><sup>0</sup> and bpf, which is an overall pre factor applied to all the coefficients resulting in the following form:

H<sub>CUB</sub> = bpf(O<sub>4</sub><sup>0</sup> + 5B<sub>4</sub><sup>4</sup>) + x\*bpf(O<sub>6</sub><sup>0</sup> - 21B<sub>6</sub><sup>4</sup>)

The Hamiltonian is diagonalied for each (x,bpf) coordinate and the eigenvalues are stored.

This approach of beginning with a cubic model and then lifting restraints was inspired by [Boothroyd's paper on PrO<sub>2</sub>](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.2082).

## Cubic Grid Creation
Create_Grid.py is used to create .mat files that contain the eigenvalues calculated by PCF, enforcing cubic constraints, for all (x,bpf).

1) Define the bounds of x and bpf as well as the number of steps between.\
	a) x = \[-1,1\] and bpf = \[-1,1\] with 200 steps each
2) Choose LS or J basis. If using LS Basis, provide LS value.
3) Create all possible (x,bpf) coordinates.
4) Diagonalize Hamiltonian for all (x,bpf) coordinates.
5) Use K-means sorting to store unique, non-zero eigenvalues. (These are the crystal field levels).
6) Save the data as a .mat file.

The user defines the range for both x and bpf as well as the number of steps for each. It's recommended to keep the dimensions equal for more subtle reasons. 

I parallelized the process for eigenvalue calculation; however, I use K-Means sorting to keep track of the eigenvalues. This is a bit intensive so this could take a while to run depending on the size of your x/bpf arrays. The idea is to only need to run this once to create and save the data, then the data can be scanned through quickly to find good starting points for the Stevens' coefficients for a given compound. Creating a 200x200 grid takes about 10 minutes to create.

Unique, non-zero eigenvalues for all (x,bpf) combinations are stored. These are the crystal field levels, so the user needs insight as what to expect. For example in the J basis for the Pr<sup>4+</sup> 6 eigenvalues are produced. With cubic constraints, either a groundstate quartet with a doublet excited state or a doublet groundstate with a quartet excited state is calculated. I figured this out by choosing a few random (x,bpf) coordinate, diagonalizing the Hamiltonian, then displaying all eigenvalues. For this J basis example I define the number of levels to be 1. PCF is used to calculate the eigenvalues and K-means is implemented to track the energy level.

In the LS basis for Pr<sup>4+</sup> there are 14 eigenvalues. I see by experimenting with a (x,bpf) coordinates that there are 4 excited levels that can occur. So the number of levels to pass to the K-means sorting function is 4.



## Cubic Grid Search
1. Declare the crystal field energy levels as measured by inelastic neutron scattering.\
	a) It's crucical to make this list in ascending order.
2. Declare an allowed tolerance between measured levels and calculated levels.\
	a) For my 200x200 grid I was able to find good starting points with 2.5% tolerance.\
	b) This tolerance may need to be bigger if the grid of (x,bpf) is not as fine. \
3. Eigenvalues and measured levels are matched within that tolerance level. If "compatible coordinates" are found that produce matching eigenvalues they are printed and these can be used as starting points for the given compound.

Grid_Search.py is used to find the starting point Stevens' Coefficients to be used for fitting. It uses the previously generated .mat file(s) to scan through each level and make sure it's within a certain allowed tolerance of that energy level as measured by inelastic neutron scattering. If compatible coordinates are found the output will be one pair of compatible (x,bpf) coordinates, as well as the diagonalized Hamiltonian as produced by the coordinates. If no compatible coordinates are found, then try increasing the tolerance level or create a finer mesh grid (more x and bpf points). It's okay if the match isn't great,  we just need decent starting points that we will later let vary.

