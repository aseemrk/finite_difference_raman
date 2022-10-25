# This script automates calculation of Resonant/off-resonant Raman spectra
# within finite-difference approach

import numpy as np 
import os 

# INPUTS:
UseMPI=True
UseOMP=True
MPI_processes=4
OMP_threads=2
prefix='Si' # Prefix as already in the pw_input etc
pw_input = 'scf.inp' #This file needs to be already available 
pw_output = 'scf.out'
pw_nscf_input = 'nscf.inp' #This file needs to be already available
pw_nscf_output = 'nscf.out'
ph_input = 'ph.inp'
ph_output = 'ph.out'
dynmat_input = 'dynmat.inp'
run_qe_ph=False # Whether to calculate the phonon eigenvectors 
run_qe_yambo=False # Whether to perform DFT for subsequent calculation of epsilon 
yambo_initialize=True # Whether to intialize Yambo databases 
run_yambo=True # Whether to run Yambo and calculate epsilon
inputs_path = '/home/aseem/work/ramfun_tests/finite_difference/Si/test_fd_script' # DFT and Yambo inputs
level='IP' # IP, RPA, TDDFT, TDDFT_LRC, BSE_SEX_TDA, BSE_SEX_bTDA
polarizations = np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1]],dtype=float) # Set the direction of photon electric field
pol_strings = ['XX','YY','ZZ','XY','XZ','YZ'] # Sttrings describing the polarizations (do not correspond to the epsilon matrix elements)
delta_r = np.array([0.0,0.001,0.002]) # Different displacements along phonon eigenvector in Angstrom
num_fd_images = len(delta_r)
bohr_ang = 0.529177249

###########################
######## functions ########
###########################

def create_scf(index_image,input_scf_name,new_atomic_pos): #Creates scf DFT input for each displacement
    # input_scf_name = template DFT scf input file name 
    # new_atomic_pos = atomic positions after diplacement
    input_file = open(input_scf_name,'r')  # Reads template file 
    lines = input_file.readlines()
    new_input_scf_name = f'scf_displacement_{index_image}.inp' # Creates new file 
    new_input_file = open(new_input_scf_name,'w')
    k=[] #Dummy variable list (useful to index lines with atomic positions)
    for i, line in enumerate(lines):
        if "ATOMIC_POSITIONS" in line: 
            new_input_file.write("ATOMIC_POSITIONS angstrom\n") # Force switches the units to angstrom
            for j in range(len(new_atomic_pos[:,0])):
                k.append(i+j+1) # Stores indices of all lines with atomic positions
                # For this dummy variable to work properly, there should be no empty line after ATOMIC_POSITIONS
        elif i in k:
            atom_name = line.split()[0]  + " " # Takes the first string as atomic symbol
            x = str(new_atomic_pos[0,0]) + " " # Writes the first atomic position  
            y = str(new_atomic_pos[0,1]) + " " #
            z = str(new_atomic_pos[0,2]) + " " #
            new_input_file.write(atom_name+x+y+z+"\n")
            k.pop(0) # IMPORTANT: first element in list k is removed
            new_atomic_pos = new_atomic_pos[1:,:] # first element in list of atomic positions is removed
        else:
            new_input_file.write(line) # Copies line from template input to new input as it is 
    input_file.close()
    new_input_file.close()
    
def create_nscf(index_image,input_nscf_name,new_atomic_pos): #Creates nscf DFT input for each displacement
    # input_nscf_name = template DFT nscf input file name
    # new_atomic_pos = atomic positions after diplacement
    input_file = open(input_nscf_name,'r') # Reads template file
    lines = input_file.readlines()
    new_input_nscf_name = f'nscf_displacement_{index_image}.inp' # Creates new file
    new_input_file = open(new_input_nscf_name,'w')
    k=[] #Dummy variable list (useful to index lines with atomic positions)
    for i, line in enumerate(lines):
        if "ATOMIC_POSITIONS" in line:
            new_input_file.write("ATOMIC_POSITIONS angstrom\n") # Force switches the units to angstrom
            for j in range(len(new_atomic_pos[:,0])):
                k.append(i+j+1) # Stores indices of all lines with atomic positions
                # For this dummy variable to work properly, there should be no empty line after ATOMIC_POSITIONS
        elif i in k:
            atom_name = line.split()[0]  + " " # Takes the first string as atomic symbol
            x = str(new_atomic_pos[0,0]) + " " # Writes the first atomic position
            y = str(new_atomic_pos[0,1]) + " " #
            z = str(new_atomic_pos[0,2]) + " " #
            new_input_file.write(atom_name+x+y+z+"\n") 
            k.pop(0) # IMPORTANT: first element in list k is removed
            new_atomic_pos = new_atomic_pos[1:,:]  # first element in list of atomic positions is removed
        else:
            new_input_file.write(line) # Copies line from template input to new input as it is
    input_file.close()
    new_input_file.close()
    
def create_yambo_input(level,polarization,path): # Creates Yambo input 
    os.system(f'cp {path}/yambo_{level}.inp yambo_{level}_temp.inp') # Copies existing input file 
    given_input = open(f'yambo_{level}_temp.inp','r')
    new_input =open(f'yambo_{level}.inp','w') # Creates new input file 
    lines = given_input.readlines()
    for i, line in enumerate(lines):
        if level == 'IP':
            if "[Xd] [cc] Electric Field" in line:
                new_input.write(f'{polarization[0]} | {polarization[1]} | {polarization[2]} |\n') # Specifies new polarization
            else:
                new_input.write(line) # Copies line in given input as it is 
        elif level == 'BSE_SEX_TDA':
            if "[BSS] [cc] Electric Field" in line:
                new_input.write(f'{polarization[0]} | {polarization[1]} | {polarization[2]} |\n') # Specifies new polarization
            else:
                new_input.write(line) # Copies line in given input as it is
    given_input.close()
    new_input.close()
    os.system(f'rm yambo_{level}_temp.inp')
            
    

###########################
########  Script   ########
###########################

if UseMPI:
    mpi_string=f'mpirun -np {MPI_processes}'

if UseOMP:
    os.system(f'export OMP_NUM_THREADS={OMP_threads}')

##                          ##
##  PART 1: QE calculatons  ##
##                          ##

# Perform a scf DFT calculation for subsequent DFPT calculation 

if run_qe_ph:
    os.system(f'{mpi_string } pw.x < {pw_input} > {pw_output}')
    print('Running DFT (pw.x) to calculate phonons')

# Now reading the output of scf DFT to get cartesian coordinates of atoms
# (just to get them readymade)

pw_out = open(pw_output,'r')
lines = pw_out.readlines() # Reading the output as lines 

for i, line in enumerate(lines):
    if "number of atoms/cell      =" in line:
        n_atom = int(line.split('=')[-1]) # Get number of atoms
        pos_atom = np.zeros((n_atom,3))   # Create array to store positions of atoms
    elif "lattice parameter (alat)  =" in line:
        alat = float(line.split()[-2]) * bohr_ang # Get the lattice parameter in Angstrom
    elif "positions (alat units)" in line:
        for j in range(n_atom):
            pos_atom[j,0] = float(lines[i+j+1].split()[-4]) * alat # Get the atomic positions
            pos_atom[j,1] = float(lines[i+j+1].split()[-3]) * alat # which are provided in 
            pos_atom[j,2] = float(lines[i+j+1].split()[-2]) * alat # alat or lattice parameter units 

pw_out.close()

# Phonon imput file 
# Do make changes if needed

ph_input_text = f"""phonons
&inputph
              prefix = '{prefix}'
              tr2_ph = 1e-16
           verbosity = 'high'
              fildyn = 'matdyn'
!            fildvscf = 'dvscf'
!     electron_phonon = 'dvscf'
               ldisp = .true.
        nq1 = 1, nq2 = 1, nq3 = 1
/"""

# Wrting the Phonon input if phonons to be computed 
if run_qe_ph:
    ph_inp = open(ph_input, "w")
    ph_inp.write(ph_input_text)
    ph_inp.close()

if run_qe_ph:
    os.system(f"{mpi_string} ph.x < ph.inp > ph.out")
    print('Running DFPT (ph.x) to calculate phonons')

# Dynmat input, this post-processing applies ASR to get correct phonons
dynmat_input_text = """&input
        fildyn = 'matdyn1'
           asr = 'crystal'
        filout = 'matdyn1.modes'
        fileig = 'matdyn1.freq'
        filmol = ''
        filxsf = 'modes.xsf'
!        q(1)=1,q(2)=0,q(3)=0
/
"""
if run_qe_ph: # Creating input file 
    dynmat_inp = open(dynmat_input, "w")
    dynmat_inp.write(dynmat_input_text)
    dynmat_inp.close()

if run_qe_ph:
    os.system(f"dynmat.x < dynmat.inp > dynmat.out")
    print('Calculating phonon eigenvectors with dynmat.x')
    
##                           ##
##  PART 2: Reading phonons  ##
##                           ##
    
# reading the ph.x output 
ph_out = open(ph_output,'r')
lines = ph_out.readlines()

# Checking the number of irreducible representations 
# We will only consider one phonon mode corresponding to 
# each of the irreducible representation 

for i, line in enumerate(lines):
    if "irreducible representations" in line:
        n_ireps=int(line.split()[-3])
        ireps = []
    elif "To be done" in line:
            ireps.append(lines[i+3].split()[2::3]) # Putting ireps in list inside list
            # So the lists inside list have indices of equivalent phonon modes
            

#Turning it into a numpy array of integers 
ireps = np.array(ireps)       
ireps = ireps.astype(int)

ph_out.close()

raman_modes = ireps[1:,0] # Excluding accoustic and taking index of one of the mode belonging to each irrep
displacement = np.zeros((len(raman_modes),n_atom,3)) # For each Raman mode, 3 cart. coords. of displacement of each atom

ph_eig_file = open("matdyn1.freq", "r") # Phonon modes and frequencies stored here 
lines = ph_eig_file.readlines()         # Note that the eigenvectors are normalized to 1 

k=-1 # Dummy index for looping over selected phonon modes
# Reading the eigenvectors now
for i, line in enumerate(lines):
    if "freq (" in line:
        if int(line[11:16]) in raman_modes:
            print(f'Selected phonon mode: {int(line[11:16])}')
            k += 1
            for j in range(n_atom):
                displacement[k,j,0] = float(lines[i+j+1].split()[1])
                displacement[k,j,1] = float(lines[i+j+1].split()[3])
                displacement[k,j,2] = float(lines[i+j+1].split()[5])

ph_eig_file.close()

##                           ##
##  PART 3: The main loop    ##
##                           ##
                
atom_pos_new = np.zeros((len(raman_modes),num_fd_images,n_atom,3)) # Creating an array of req. dims

for j in range(len(raman_modes)): # Loop Raman modes
    os.system(f'mkdir mode_{raman_modes[j]}') # Create a directory for each Raman mode
    for i in range(0,num_fd_images): #Loop displacements
        os.system(f'mkdir mode_{raman_modes[j]}/displacement_{i}') # Create directory for each displacement
        atom_pos_new[j,i,:,:] = pos_atom[:,:] + displacement[j,:,:] * delta_r[i] # Displaced positions
        # Calculation of phonons
        if run_qe_yambo: 
            create_scf(i,pw_input,atom_pos_new[j,i,:,:]) # Create scf input 
            create_nscf(i,pw_nscf_input,atom_pos_new[j,i,:,:]) # Create nscf input
        os.chdir(f'mode_{raman_modes[j]}/displacement_{i}')
        # Calculation of DFT g.s. for Yambo
        if run_qe_yambo:
            os.system(f'mv ../../*scf_displacement_{i}.inp .') # Move the scf/nscf inputs 
            os.system(f'cp ../../*.upf .') # Copy the pseudos
            os.system(f'{mpi_string} pw.x < scf_displacement_{i}.inp > scf_displacement_{i}.out') 
            os.system(f'{mpi_string} pw.x < nscf_displacement_{i}.inp > nscf_displacement_{i}.out')
        os.chdir(f'{prefix}.save')
        # Note: pw.x will overwrite if there is a previous DFT calculation
        # If new DFT calculation, one can create Yambo database
        if yambo_initialize:
            os.system(f'rm -r SAVE r*_setup*')
            os.system('p2y') # Dumping QE output to Yambo database
            os.system('yambo') # Yambo summerizes the database for itself
        if run_yambo:
            os.system('rm SAVE/ndb.em1s*')   # Clear old screening
            for p, pol in enumerate(polarizations): # Loop polarizations
                create_yambo_input(level,pol,inputs_path) # Creates Yambo input 
                pol_string = pol_strings[p] 
                os.system(f'mkdir {level}_{pol_string}')
                os.system(f'rm -r *{level}_{pol_string}*') # Clear old calculations
                os.system(f'mv yambo_{level}.inp yambo_{level}_{pol_string}.inp') # Move created inputs
                os.system(f'{mpi_string} yambo -F yambo_{level}_{pol_string}.inp -J {level}_{pol_string}')
                os.system(f'mv o-{level}_{pol_string}*eps* o-{level}_{pol_string}.eps ')
                if p == 0:
                    os.system(f'mv {level}_{pol_string}/ndb.dipoles* SAVE') # Copy dipoles for reuse
                    os.system(f'mv {level}_{pol_string}/ndb.em1s* SAVE') # Copy static screening for reuse
                
                
        os.chdir('../../..')
        
############ DFT and Yambo calculations finished ##################


        
        

