## The Cubic Crystal Field Hamiltonian
Cubic symmetry allows for a simple Crystal Field Hamiltonian which consists of only 2 independent Stevens Coefficients:

### H<sub>CUB</sub> = B<sub>4</sub>(O<sub>4</sub><sup>0</sup> + 5B<sub>4</sub><sup>4</sup>) + B<sub>6</sub>(O<sub>6</sub><sup>0</sup> - 21B<sub>6</sub><sup>4</sup>)

I parameterize this hamiltonian by x = B<sub>6</sub><sup>0</sup>/B<sub>4</sub><sup>0</sup> and bpf, which is an overall pre factor applied to all the coefficients resulting in the following form:

### H<sub>CUB</sub> = bpf(O<sub>4</sub><sup>0</sup> + 5B<sub>4</sub><sup>4</sup>) + x\*bpf(O<sub>6</sub><sup>0</sup> - 21B<sub>6</sub><sup>4</sup>)



## Cubic Grid Creation
Create_Grid.py is used to create .mat files that contain the eigenvalues calculated by PCF, enforcing cubic constraints, for all (x,bpf).
