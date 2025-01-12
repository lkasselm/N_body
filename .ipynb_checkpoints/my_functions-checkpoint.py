# random collection of useful functions
import numpy as np

def cm(snap):
    # return center of mass from a pynbody SimSnap
    cm_x = np.sum(snap['mass'] * snap['x'])/np.sum(snap['mass'])
    cm_y = np.sum(snap['mass'] * snap['y'])/np.sum(snap['mass'])
    cm_z = np.sum(snap['mass'] * snap['z'])/np.sum(snap['mass'])
    return cm_x, cm_y, cm_z

def recenter_snap(snap):
    # recenter a pynbody SimSnap (so that its center of mass is zero)
    cm_x, cm_y, cm_z = cm(snap)
    snap['x'] -= cm_x
    snap['y'] -= cm_y
    snap['z'] -= cm_z
    
# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

import h5py
def write_IC_hdf5(filename, data):
 def ifset(d,key):
  if key in d.keys():
   result=d[key]
  else:
   result=None
   if key in ['lt','fmt']:  #for plot
    result=''
  return result

 def ifset2(d,key,value):
  if value is None:
   result=ifset(d,key)
  else:
   result=value
  return result

 if isinstance(data, dict):
  data=[data]

 BoxSize = None
 NumPartType = 6
 MassTable = np.zeros(NumPartType, dtype = float)
 NumPart = np.zeros(NumPartType, dtype = int)
 
 for d in data:
  BoxSize = ifset2(d, 'BoxSize', BoxSize)
  i = d['PartType']
  MassTable[i] = d['PartMass']
  NumPart[i] = d['count']
  
 file = h5py.File(filename+'.hdf5','w')

 group = file.create_group("Header")
 if BoxSize is not None:
  group.attrs.create("BoxSize", BoxSize, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 else:
  group.attrs.create("BoxSize", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Flag_Cooling", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_DoublePrecision", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_Feedback", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_IC_Info", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_Metals", 0, shape=None, dtype=h5py.h5t.STD_I32LE) #in makegal ics is 1
 group.attrs.create("Flag_Sfr", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_StellarAge", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("HubbleParam", 1, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("MassTable", MassTable, shape=None,
                                                     dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("NumFilesPerSnapshot", 1, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("NumPart_ThisFile", NumPart, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("NumPart_Total", NumPart, shape=None, dtype=h5py.h5t.STD_U32LE)
 group.attrs.create("NumPart_Total_HighWord", (0,0,0,0,0,0), shape=None, dtype=h5py.h5t.STD_U32LE)
 group.attrs.create("Omega0", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("OmegaLambda", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Redshift", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Time", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)

 ID_offset = 0
 for d in data:
  group = file.create_group("PartType"+str(d['PartType']))
  dataset = group.create_dataset("Coordinates", (d['count'],3), data=d['Coordinates'],
                                                            dtype=h5py.h5t.IEEE_F32LE)
  dataset = group.create_dataset("ParticleIDs", (d['count'],),
                   data=np.array(range(ID_offset,ID_offset+d['count'])), dtype=h5py.h5t.STD_I32BE)
  ID_offset += d['count']
  dataset = group.create_dataset("Velocities", (d['count'],3), data=d['Velocities'],
                                                            dtype=h5py.h5t.IEEE_F32LE)
  if d["PartType"] == 0:
   dataset = group.create_dataset("InternalEnergy", (d['count'],),
                   data=d["InternalEnergy"], dtype=h5py.h5t.IEEE_F32LE)
 file.close()

def read_IC_hdf5(filename):
 data = {}
 file = h5py.File(filename,'r')
 keys = list(file.keys())
 for key in keys:
  if key == 'Header':
   header = dict(file['Header'].attrs.items())
   data['Header'] = header
  else:
   partData = file.get(key)
   partKeys = list(partData.keys())
   data[key] = {}
   for partKey in partKeys:
    data[key][partKey] = np.array(partData.get(partKey))
 file.close()
 return data
    
def append_to_new_line(filename, text):
    # Open the file in append mode
    with open(filename, 'a') as file:
        # Add a newline character before the text if the file isn't empty
        file.write('\n' + text)

def dens_profile(x, y, z, mass, recenter=False):
    import pynbody as pn
    import numpy as np
    from matplotlib import pyplot as plt
    """
    This is a wrapper for the pynbody profile function
    which can be used to plot a 3d density profile from some 
    coordinate and mass arrays.
    """
    N_part = len(x)
    #coords_ = np.transpose(coords)
    data = pn.new(n_particles = N_part)
    data['x'] = x
    data['y'] = y
    data['z'] = z
    data['mass'] = np.ones(N_part) * mass
    if recenter:
        recenter_snap(data)
    p = pn.analysis.profile.Profile(data, ndim=3, type='log')
    plt.loglog(p['rbins'], p['density'])
    plt.loglog(p['rbins'], Hernquist(p['rbins']))
    plt.xlabel('r [kpc]')
    plt.ylabel(r'$\rho$ [M$_\odot$ kpc$^{-3}$]')

def Hernquist(r): # Hernquist density profile 
    return 1/(r * (1+r)**3)

import glob
import os
def get_num_snaps(sim_directory):
    # Get number of snapshots in directory.
    files = glob.glob(os.path.join(sim_directory, 'snapshot*'))
    return len(files)
    
def get_snapname(snap_num):
    # return the filename of a given snapshot
    # assumes that snapshots are in the format 
    # 'snapshot_xxx.hdf5'
    num_str = str(snap_num)
    if len(num_str) == 1:
        num_str = '00' + num_str
    elif len(num_str) == 2:
        num_str = '0' + num_str
    snapname = 'snapshot_'+num_str+'.hdf5'
    return snapname

def get_time(snapname):
    return read_IC_hdf5(snapname)['Header']['Time']

import pynbody as pn
from matplotlib import pyplot as plt

def plot_center_of_mass(sim_directory, num_snaps=None, do_plot=True, species=None):
    # plot center of mass of the system over time
    if num_snaps is None:
        num_snaps = get_num_snaps(sim_directory)
    time = []
    cmx = []
    cmy = []
    cmz = []
    for i in range(num_snaps):
        snapname = get_snapname(i)
        path = sim_directory+snapname
        data = read_IC_hdf5(path)
        time.append(data['Header']['Time'])
        snap = pn.load(path)
        if species == 'gas':
            snap = snap.gas
        elif species == 'dm':
            snap = snap.dm
        cm_ = cm(snap)
        cmx.append(cm_[0])
        cmy.append(cm_[1])
        cmz.append(cm_[2])
    if do_plot:
        plt.plot(time, cmx, label='x')
        plt.plot(time, cmy, label='y')
        plt.plot(time, cmz, label='z')
        plt.xlabel('Time [Gyr]')
        plt.ylabel('Center of mass [kpc]')
        plt.legend()
        plt.show()
    return time, cmx, cmy, cmz

def plot_momentum(sim_directory, num_snaps=None, do_plot=True, species=None):
    # plot total momentum over time
    if num_snaps is None:
        num_snaps = get_num_snaps(sim_directory)
    time = []
    px = []
    py = []
    pz = []
    for i in range(num_snaps):
        snapname = get_snapname(i)
        path = sim_directory+snapname
        data = read_IC_hdf5(path)
        time.append(data['Header']['Time'])
        snap = pn.load(path)
        if species == 'gas':
            snap = snap.gas
        elif species == 'dm':
            snap = snap.dm
        px.append(np.sum(snap['mass']*snap['vx']))
        py.append(np.sum(snap['mass']*snap['vy']))
        pz.append(np.sum(snap['mass']*snap['vz']))
    if do_plot:
        plt.plot(time, px, label='x')
        plt.plot(time, py, label='y')
        plt.plot(time, pz, label='z')
        plt.xlabel('Time [Gyr]')
        plt.ylabel(r'Total momentum [$10^{10}$ M$_\odot$ km/s]')
        plt.legend()
        plt.show()
    return time, px, py, pz

def get_energy(sim_directory):
    # read energy.txt file and return as dictionary
    data = np.loadtxt(sim_directory+'energy.txt')
    data = np.transpose(data)
    output =  {}
    tot = {}
    gas = {}
    dm = {}
    output['time'] = data[0]
    tot['int'], tot['pot'], tot['kin'] = data[1], data[2], data[3]
    gas['int'], gas['pot'], gas['kin'] = data[4], data[5], data[6]
    dm['int'], dm['pot'], dm['kin'] = data[7], data[8], data[9]
    tot['tot'] = tot['int'] + tot['pot'] + tot['kin']
    gas['tot'] = gas['int'] + gas['pot'] + gas['kin']
    dm['tot'] = dm['int'] + dm['pot'] + dm['kin']
    output['tot'] = tot
    output['gas'] = gas
    output['dm'] = dm
    return output

def plot_total_energy(sim_directory):
    # plot total energy over time
    data = get_energy(sim_directory)
    plt.plot(data['time'], data['gas']['tot']/np.abs(data['gas']['tot'][0]), label='Gas')
    plt.plot(data['time'], data['dm']['tot']/np.abs(data['dm']['tot'][0]), label='DM')
    plt.legend()
    plt.xlabel('Time [Gyr]')
    plt.ylabel(r'$E(t)/|E(t_0)|$')
    plt.show()

def plot_energy_components(sim_directory, key):
    data_ = get_energy(sim_directory)
    time = data_['time']
    data = data_[key]
    plt.plot(time, data['int']/np.abs(data['tot']), label=r'$E_{int}/|E_{tot}|$')
    plt.plot(time, data['kin']/np.abs(data['tot']), label=r'$E_{kin}/|E_{tot}|$')
    plt.plot(time, data['pot']/np.abs(data['tot']), label=r'$E_{pot}/|E_{tot}|$')
    plt.legend()
    plt.xlabel('Time [Gyr]')
    plt.show()

keys = {'gas':'PartType0', 'dm':'PartType1'}
def available_keys(datapath, species='gas'):
    data = read_IC_hdf5(datapath)[keys[species]]
    return data.keys()

def center_of_mass(coords, numpart):
    c = np.transpose(coords)
    cx = np.sum(c[0])/numpart
    cy = np.sum(c[1])/numpart
    cz = np.sum(c[2])/numpart
    return np.array([cx, cy, cz])
    
def scatter_plot(datapath, species='gas', 
                 qty='Density', xscale='log',
                 yscale='log', s=1, xlabel='r [kpc]', 
                 ylabel='',
                 recenter=True, N_sample=None, label=None,
                 savepath=None, xlim=None, ylim=None):
    # do a scatter plot of radius vs some quantity
    data = read_IC_hdf5(datapath)
    numpart_list = data['Header']['NumPart_ThisFile']
    if species == 'gas':
        numpart = numpart_list[0]
    elif species == 'dm':
        numpart = numpart_list[1]
    data = data[keys[species]]
    coords = data['Coordinates']
    if recenter:
        cm = center_of_mass(coords, numpart)
        coords = coords - cm
    c = np.transpose(coords)
    r = np.sqrt(c[0]**2 + c[1]**2 + c[2]**2)
    if N_sample is not None:
        # only plot a random subset of datapoints
        # You should always use this if 
        # the number of particles is large 
        # to keep the vector graphic size small
        ind = np.random.choice(len(r), N_sample, 
                               replace=False)
    else:
        ind = len(r)
    plt.scatter(r[ind], data[qty][ind], s=s, label=label)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

import h5py
# scatter plot with lazy loading
def scatter_plot_2(datapath, species='gas', 
                 qty='Density', xscale='log',
                 yscale='log', s=1, xlabel='r [kpc]', 
                 ylabel='',
                 recenter=True, N_sample=None, label=None,
                 savepath=None, xlim=None, ylim=None):
    # do a scatter plot of radius vs some quantity
    file = h5py.File(datapath, 'r')
    header = dict(file['Header'].attrs.items())
    #data = read_IC_hdf5(datapath)
    numpart_list = header['NumPart_ThisFile']
    if species == 'gas':
        numpart = numpart_list[0]
    elif species == 'dm':
        numpart = numpart_list[1]
    #data = data[keys[species]]
    coords = file[keys[species]]['Coordinates'][:]
    #coords = data['Coordinates']
    if recenter:
        cm = center_of_mass(coords, numpart)
        coords = coords - cm
    c = np.transpose(coords)
    r = np.sqrt(c[0]**2 + c[1]**2 + c[2]**2)
    if N_sample is not None:
        # only plot a random subset of datapoints
        # You should always use this if 
        # the number of particles is large 
        # to keep the vector graphic size small
        ind = np.random.choice(len(r), N_sample, 
                               replace=False)
    else:
        ind = len(r)
    plt.scatter(r[ind], data[qty][ind], s=s, label=label)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
def calc_rel_dens_err(datapath, des_profile='Hernquist', 
                      r_s = 0.39, rho_s = 0.26, return_rms=True): 
    if des_profile == 'Hernquist':
        def des_dens(r):
            x = r/r_s
            return rho_s/(x * (1+x)**3)
    data = read_IC_hdf5(datapath)['PartType0']
    coords = data['Coordinates']
    c = np.transpose(coords)
    r = np.sqrt(c[0]**2 + c[1]**2 + c[2]**2)
    des_densities = des_dens(r)
    actual_densities = data['Density']
    rel_dens_err = (actual_densities - des_densities)/des_densities
    # return rms 
    return np.sqrt(np.mean(rel_dens_err**2))

def calc_rel_dens_errs(simpath, 
                       des_profile='Hernquist', 
                       r_s = 0.39, rho_s = 0.26, 
                       return_rms=True, num_snaps=None, mult=100):
    if num_snaps is None:
        num_snaps = get_num_snaps(simpath)
    times = []
    rdes = []
    for i in range(num_snaps):
        if i % mult == 0:
            snapname = simpath + get_snapname(i)
            times.append(get_time(snapname))
            rdes.append(calc_rel_dens_err(snapname, des_profile=des_profile, r_s=r_s, rho_s=rho_s, return_rms=return_rms))
    return times, rdes

def NFW_Mvir_c_to_rs_rhos(Mvir, c, rho_c=130, Delta=100):
    """
    Inputs:
    * Virial Mass of the halo in solar masses
    * concentration parameter of the halo
    * critical density of the universe (optional, default 130 M_solar/kpc^3)
    * Overdensity which defines the virial radius (optional, default: 100)

    Outputs:
    * scale radius in kpc
    * scale density in M_solar/kpc^3
    """
    rho_s = Delta * rho_c * c * (1+c)**3
    r_s = (Mvir/(4*np.pi*rho_s*( np.log(1+c) - c/(1+c)) ))**(1/3)
    return r_s, rho_s

    
    