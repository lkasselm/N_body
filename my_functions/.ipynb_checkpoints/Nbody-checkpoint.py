# functions related to analyzing N-body simulation data
import numpy as np
import pynbody as pn
from matplotlib import pyplot as plt

def read_data(path):
    """
    Read profile data stored in csv file

    Inputs:
    * path (string): path to the file

    Outputs:
    * (Dict): Data as dictionary object
    """
    dict = {}
    dict['header'] = {}
    # get header:
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#') and '=' in line:
                key, value = line.strip('#').split('=')
                key = key.strip()
                value = value.strip()
                if key in ('MVIR','RVIR','z','h0'):
                    dict['header'][key] = float(value)
                else:
                    break
    Mcode_to_physMsun = 1e10/dict['header']['h0']
    Dcode_to_physKpc = 1./dict['header']['h0']/(1.+dict['header']['z'])
    Rhocode_to_phys = Mcode_to_physMsun/Dcode_to_physKpc**3
    # get data:
    data = np.loadtxt(path, delimiter=',', dtype='float', comments='#', usecols=(0,1,2,3,4,5,6,7,8,9,10,11))
    DM_r0,DM_r1,DM_r,DM_rho,Bar_r0,Bar_r1,Bar_r,Bar_rho,Tot_r0,Tot_r1,Tot_r,Tot_rho = data.T
    dict['DM_r'] = DM_r 
    dict['DM_rho'] = DM_rho 
    dict['Bar_r'] = Bar_r 
    dict['Bar_rho'] = Bar_rho 
    dict['Tot_r'] = Tot_r 
    dict['Tot_rho'] = Tot_rho 
    return dict

from scipy.integrate import simpson
def get_mass(r_arr, dens_arr, r):
    """
    Calculate total mass enclosed in r

    Inputs:
    * r_arr (Array, radii)
    * dens_arr (Array, densities)
    * r (float)

    Outputs:
    * M (float)
    """
    # remove all nans:
    r_arr = np.nan_to_num(r_arr, nan=0)
    dens_arr = np.nan_to_num(dens_arr, nan=0)
    M_arr = np.zeros_like(dens_arr)
    dMdr = 4*np.pi*dens_arr*r_arr**2

    # get index of r_arr closest to r:
    idx = np.argmin(np.abs(r_arr-r))
    if idx == 0:
        return 0
    else:
        return simpson(dMdr[:idx], x=r_arr[:idx])

def mean_density(r_arr, rho_arr, r):
    """
    Return the mean density within r
    
    Inputs:
    * r_arr (Array, radii)
    * dens_arr (Array, densities)
    * r (float)
    """
    encl_mass = get_mass(r_arr, rho_arr, r)
    return 3*encl_mass/(4*np.pi*r**3)

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

def dens_profile(x, y, z, mass, recenter=False):
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