import numpy as np
import matplotlib.pyplot as plt 
import os

# This script works if the dielectric tensor is symmetric about the diagonal  
# Run this script in the same directory as the other script 

# Inputs
level='IP' # IP, RPA, TDDFT, TDDFT_LRC, BSE_SEX_TDA, BSE_SEX_bTDA
mode = '4' # Phonon mode 
polarizations = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]],dtype=float)
pol_strings = ['XX','YY','ZZ','XY','XZ','YZ']
n_points = 5000 # Number of points on energy scale (controlled by Yambo)
x_array = np.array([0.0,0.001,0.002]) # Displacements along phonon eigenvector
# Inputs end

fd_image = np.arange(0,len(x_array),1,dtype=int)
colours = ['red','blue','green','turqouise','brown','brown','darkorange','black']
spectra_data = np.zeros((len(polarizations),len(fd_image),n_points,3)) # stores yambo comuted spectra 
#reading data from Yambo output
for j, pol in enumerate(polarizations):
    for i in fd_image:
        file_string = f'mode_{mode}/displacement_{i}/Si.save/o-{level}_{pol_strings[j]}.eps'
        spectra_data[j,i,:,0] = np.loadtxt(file_string,usecols=(0))
        spectra_data[j,i,:,1] = np.loadtxt(file_string,usecols=(1))
        spectra_data[j,i,:,2] = np.loadtxt(file_string,usecols=(2))
    
epsilon = np.zeros((len(polarizations),len(fd_image),n_points,3))
# The diagonal elements of epsilon 
epsilon[0,:,:,:] = spectra_data[0,:,:,:] # The diagonal elements 
epsilon[1,:,:,:] = spectra_data[1,:,:,:]
epsilon[2,:,:,:] = spectra_data[2,:,:,:]
# The off-diagonal elements of epsilon
epsilon[3,:,:,:] = 0.5 * (4.0 * spectra_data[3,:,:,:] - spectra_data[0,:,:,:] - spectra_data[1,:,:,:] )
epsilon[4,:,:,:] = 0.5 * (4.0 * spectra_data[4,:,:,:] - spectra_data[0,:,:,:] - spectra_data[2,:,:,:] ) 
epsilon[5,:,:,:] = 0.5 * (4.0 * spectra_data[5,:,:,:] - spectra_data[1,:,:,:] - spectra_data[2,:,:,:] )

eps_slope_real = np.zeros((len(polarizations),n_points)) # arrays to store the derivatives
eps_slope_imag = np.zeros((len(polarizations),n_points)) # of real and imag parts

#calculating slope using central difference formula (higher accuracy)
for j in range(len(polarizations)):
    for i in range(n_points):
        eps_slope_imag[j,i] = ( epsilon[j,1,i,1] - epsilon[j,0,i,1] ) * (1.0/(x_array[1]-x_array[0])) 
        eps_slope_real[j,i] = ( epsilon[j,1,i,2] - epsilon[j,0,i,2] ) * (1.0/(x_array[1]-x_array[0])) 
    
alpha = np.zeros((len(polarizations),n_points)) # Square of Raman Tensor 
alpha  = eps_slope_imag**2 + eps_slope_real**2  # (proportional to scattering rate)

fig, ax = plt.subplots()

for i in range(len(polarizations)):
    np.savetxt(f'alpha_{level}_{pol_strings[i]}.dat',np.c_[spectra_data[0,0,:,0],alpha[i,:]])
    ax.plot(spectra_data[0,0,:,0],alpha[i,:],label=f'{pol_strings[i]}')

plt.xlim(left=2.5,right=4.5)
#plt.yscale('log')
#plt.ylim(bottom=10000)
plt.legend()
plt.savefig(f'raman_{level}.pdf')
os.system(f'evince raman_{level}.pdf')    
