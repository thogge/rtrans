import numpy as np
import pdb
import matplotlib.pyplot as plt
import glob, os
import os

def makemap():

    '''
    PURPOSE:

    Returns an array of values of tau for each location in the collapsing
    sphere based on the outputs from radex. 

    Parameters for the model are defined below

    '''


    #Set parameters
    size = 11
    mole = 'co'


    #Set up the grid
    b = np.tile(np.arange(size), (size,1))
    x = np.transpose(b)

    r = np.sqrt(b**2 + x**2)


    #Calculate density and temperature
    T = Tfunc(r)

    den = denfunc(r)
    n = den[0]
    N = den[1]

    #Make the job files
    for bind in b[0]:
        for xind in x[:,0]:
            input(mole, T[xind, bind], N[xind, bind], n[xind, bind], bind, xind)

    #Run the job files
    run()

    tau = extract(size, mole)

    pdb.set_trace()
    

def extract(size,mole,path ='/Users/connorr/Desktop/Radex/bin/jobfiles/'):
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

        #unsorted.append([temptau,float(bind),float(xind)])
        
        
        temptau =  np.loadtxt(f, skiprows = skip, usecols = [7])
        tau[xind,bind] = temptau[line]
 #       os.system('rm '+f)
        os.system('rm '+f[0:-3]+'inp')

    #Remove all the files that did not run.
    os.system('rm '+path+'*inp')

    return tau
    
def denfunc(r, nc = 1e2, rc = 3, abundance = 1e-4):

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
    


    
def input(mole,tkin,N,den,b,x,path ='/Users/connorr/Desktop/Radex/bin/jobfiles/', trad = 2.73):

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
    path: Path containing the output job files. By default set to /Users/connorr/Desktop/Radex/bin/jobfiles/
    trad: Radiation temperature, by default set to 2.73.
    '''
    
    #Makes the jobfiles
    
    linewidth = 0


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
    infile.write(str(1.0)+'\n')
    infile.write(str(linewidth)+'\n')


def run(path= '/Users/connorr/Desktop/Radex/bin/jobfiles/'):
    #Runs ALL the jobfiles in the path at once
    filelist = glob.glob(path+'*.inp')

    for file in filelist:
        os.system('radex < '+file)
        os.system('mv '+ os.getcwd()+'/'+str.split(file[0:-3],path)[1]+'out '+file[0:-3]+'out')
        


