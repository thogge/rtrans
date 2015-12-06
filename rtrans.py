import numpy as np
import pdb
import matplotlib.pyplot as plt
import glob, os
import os




def rtrans():
    '''
    PURPOSE:
    Does radiative transfer for a collapsing cloud

    '''
    
    size = 11
    mole = 'co'
    massmole = 28*1.67e-24#g
    nc = 1e3#cm^-3
    rc = 1#pc
    abundance = 1e-15#molecules/H2


    nures = 100
    nu_0 = 115e9#Hz (Using CO1-0


    #Set up coordinates
    b = np.tile(np.arange(size), (size,1))
    x = np.transpose(b)
    r = np.sqrt(b**2 + x**2)

    T = Tfunc(r, power = 2)
    den = denfunc(r, nc = nc, rc = rc, abundance = abundance)
    n = den[0]
    N = den[1]

    Rout = len(r)-1
    
    quad = makemap(size, mole, nc, rc, abundance)

    #Flips and copies the tau map
    taumap = np.hstack((quad[:,::-1], quad))

    vel = make_vel()
    velmap = np.hstack((-vel[:,::-1], vel))

    T = Tfunc(r, power = 2)
    Tmap = np.hstack((T[:,::-1], T))

    

    delvel =1

    velrange = np.arange(0, nures, delvel)

    gauss(taumap, Tmap, velmap, nu_0,  massmole, velrange)
    
    #Get tau as a function of nu via velocity
    
    
def gauss(tau, T, v, nu_0, massmole, velrange):
    '''
    PURPOSE:

    Broadens tau into a gaussian and shifts it based on velocity

    INPUTS:

    tau: Array of tau. Not yet a function of nu

    T: Array of temperatures in the cloud

    v: Array of radial velocities in the cloud

    nures: Resolution of the spectra

    nu_0: Rest value of nu

    OUTPUTS:

    3D Map of Tau with velocity (frequency)  dependance on the Z axis

    '''

    k = 1.38e-16#cgs units
    
    pdb.set_trace()
    
    taunu = np.zeros((np.shape(tau)[0], np.shape(tau)[1], len(velrange)))

    vtherm = np.sqrt(2*k*T/massmole) 
    
    for t,index in self.tau.flat:      
        taunu(index) = 1/(vtherm(index)*sqrt(np.pi)) * np.exp()
    
    
    
    

def Inuout(Iin, tau, T, nu):
    '''
    PURPOSE:
    
    Updates the intensity in a cell as a function of nu

    INPUT:
    
    Iin: Incoming intensity as a function of nu

    tau: Optical depth as a function of nu

    T: Kinetic temperature of the gas
    
    nu: Array containing nu

    '''
    
    Iout = Iin*exp(-tau) + Bnu(T,nu)*(1-np.exp(-tau))
    
    return Iout
    

def Bnu(T, nu):

    h = 6.63e-27#erg s
    c = 3e10#cm/s
    k = 1.38e-16#erg K^-1
    
    B = 2*h*nu**3/c**2/(np.exp(h*nu/(k*T))-1)

    return B

def makemap(size, mole, nc, rc, abundance):

    '''
    PURPOSE:

    Returns an array of values of tau for each location in the collapsing
    sphere based on the outputs from radex. 

    Parameters for the model are defined below

    '''

    #Path to jobfiles
    path = '/Users/connorr/Desktop/Radex/bin/jobfiles/'

    #Set up the grid
    b = np.tile(np.arange(size), (size,1))

    x = np.transpose(b)

    r = np.sqrt(b**2 + x**2)


    #Calculate density and temperature
    T = Tfunc(r, power = 2)
    
    den = denfunc(r, nc = nc, rc = rc, abundance = abundance)
    n = den[0]
    N = den[1]

    #Make the job files
    for bind in b[0]:
        for xind in x[:,0]:
            input(mole, T[xind, bind], N[xind, bind], n[xind, bind], bind, xind, path)

    #Run the job files
    run(path = path)

    tau = extract(r,size, mole, path)


    #Assume that tau = 0 outside of sphere    
    sphere = -(r-size+1)>0

    return tau*sphere
    
def vel(r,m_mol,T_in,T_out,gamma):

    # r is in units of cloud radius
    k = 1.38e-16
    v = np.sqrt(-(4*k/m_mol)*(T_in*np.log(r)+((T_in-T_out)/(gamma*(gamma+1)))*(r**gamma-1)))
    return v

def vel_rad(v_center,x,b):

    r = np.sqrt(x**2+b**2)
    v_rad = v_center*x/r
    return v_rad 
    
def make_vel():
    k = 1.38e-16
    m_mol = 28*1.67e-24
    T_in = 100
    T_out = 10
    gamma = 2
    x_range = 25
    y_range = 25
    r_max = np.sqrt(x_range**2+y_range**2)
    vel_arr = np.zeros((y_range,x_range))
    for i in range(0,y_range):
        for j in range(1,x_range):
            r = np.sqrt(i**2+j**2)
            vel_arr[i,j] = vel_rad(vel(r/r_max,m_mol,T_in,T_out,gamma),j,i)
    #plt.pcolor(np.arange(x_range),np.arange(y_range),vel_arr/1e5)
    #plt.colorbar()
    #plt.show()
    return vel_arr

def extract(r,size,mole,path):
    '''
    PURPOSE:
    
    Extracts tau from the jobfiles and then returns an array.
    
    INPUTS:
    
    size: Length of one size of R array

    mole: String containing name of molecule e.g. 'co'

    OPTIONAL INPUTS:
    
    Path to jobfiles

    OUTPUTS:

    Array containing all the information about tau ordered by x and b.

    '''

    fnames = glob.glob(path+'*.out')

    #Set the number of rows to skip
    #Seems to be right number for co, but I do remember something about
    #it changing between molecules so may want to check the .out files
    skip = 13

    #Read in data
    unsorted = []
    tau = np.zeros((size,size))

    #For now only grab the first line
    line = 0
    
    for f in fnames:
        #Grab the impact parameter and x.
        bind = f.split(path+mole+'_')[1].split('_')[0]
        xind = f.split(bind+'_')[1].split('.')[0]

        temptau =  np.loadtxt(f, skiprows = skip, usecols = [7])
        tau[xind,bind] = temptau[line]

        #Remove the jobfiles
        os.system('rm '+f)
        os.system('rm '+f[0:-3]+'inp')

    #Remove all the files that did not run.
    #Note: If there are more than 10000 files this will fail.
    os.system('rm '+path+'*inp')

    return tau
    
def denfunc(r, nc, rc, abundance):

    '''
    PURPOSE:

    Returns a density field given by the Bonnie-Ebert sphere for a given radius

    INPUTS:

    r: array of radii

    OPTIONAL INPUTS:
    
    rhoc: Central number density in cm-3

    rc: Cloud radius in pc

    abundance: Specify the abundance relative to H2.

    OUTPUT:

    Two arrays, one containing the volume density at each location and the other containing 
    the column density at each location.

    '''
    
    pc = 3.09e18#cm
    
    rout = len(r)-1

    n = nc * (rout**2/(rout**2 + r**2))

    l = (rc * pc)/rout

    N = n*l*abundance

    return (n,N)
    

def Tfunc(r, power = 2, Tin=100, Tout=10, constant = False):

    '''
    PURPOSE:

    Returns a value for temperature as a function of radius.

    INPUTS:

    r: array of radii

    OPTIONAL INPUTS:
    
    power: Power law form of the temperature of the cloud.

    Tin: Inner temperature of the cloud. set to 100 by default.
    
    Tout: Outer temperate of the cloud. Set to 10 by default.

    OUTPUT:

    Array containing the temperature at each location.

    '''
    
    rout = len(r)-1
    
    T = Tin - ((Tin-Tout)/(rout**power))*r**power

    return T
    


    
def input(mole,tkin,N,den,b,x,path, trad = 2.73):

    '''
    PURPOSE: 

    Creates input files for radex

    INPUTS:
    
    mole: String containting the molecule you want to model, currently 3 available, 'hco+', 'co', and 'hcn'. 

    tkin: Kinetic temperate of the cell

    N: Column density of the cell

    den: Volume density of the cell

    b: Impact parameter from the center of the cloud. Should be in the form of the index of the resolution element
       Define the line of sight through the center of the cloud as b=0

    x: x coordinate measured from the center of the cloud. Again should be in the form of the index of a resolution element.

    OPTIONAL INPUTS:
    path: Path containing the output job files. 
    trad: Radiation temperature, by default set to 2.73.
    '''
    
    #Makes the jobfiles
    
    linewidth = 4


    #format b with 5 significant digits
    bstr = str(b)
    places = 5    
    bzero = '0'*(places-len(bstr))+bstr

    xstr = str(x)
    xzero = '0'*(places-len(xstr))+xstr

    
    filename = mole+'_'+bzero+'_'+xzero
    
    infile = open(path+filename+'.inp', 'w')

    #Set the range to look for lines for a specific molecule
    if mole == 'hco+':
        minmax = ['50','1400']
    elif mole == 'co':
        minmax = ['10', '1000']
    elif mole == 'hcn':
        minmax = ['1','1400']
    
    Densci = "{:.2E}".format(den)
    Nsci = "{:.2E}".format(N)
        
    infile.write(mole+'.dat\n')
    infile.write(filename+'.out\n')
    infile.write(minmax[0]+' '+minmax[1]+'\n')
    infile.write(str(tkin)+'\n')
    infile.write('1\n')
    infile.write('H2\n')
    infile.write(str(den)+'\n')
    infile.write(str(trad)+'\n')
    infile.write(str(Nsci)+'\n')
    infile.write(str(linewidth)+'\n')
    infile.write(str(0)+'\n')


def run(path):
    #Runs ALL the jobfiles in the path at once
    filelist = glob.glob(path+'*.inp')

    for file in filelist:

        os.system('radex < '+file)
        os.system('mv '+ os.getcwd()+'/'+str.split(file[0:-3],path)[1]+'out '+file[0:-3]+'out')
        


