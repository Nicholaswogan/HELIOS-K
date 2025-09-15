import numpy as np
from matplotlib import pyplot as plt
import h5py
import os
import numba as nb
from numba import types
from wogan_data import bins as wogan_bins

# @nb.njit()
def remove_duplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

@nb.experimental.jitclass()
class KFile():
    weights : types.double[:]
    xs: types.double[:]
    T : types.double[:]
    P : types.double[:]
    weightsn : types.int32[:]
    def __init__(self, n):
        self.weights = np.empty(n)
        self.xs = np.empty(n)
        self.T = np.empty(n)
        self.P = np.empty(n)
        self.weightsn = np.empty(n,np.int32)

@nb.experimental.jitclass()
class KGrid():
    g : types.double[:]
    ng : types.int32[:]
    P : types.double[:]
    T : types.double[:]
    def __init__(self, k):
        self.g = np.unique(k.weights)
        self.ng = np.unique(k.weightsn)
        self.P = np.unique(k.P)
        self.T = np.unique(k.T)

@nb.njit()
def helper(k, g, i):
    iT = np.argmin(np.abs(k.T[i]-g.T))
    assert np.isclose(k.T[i],g.T[iT],rtol=1e-8,atol=1e-20)

    iP = np.argmin(np.abs(k.P[i]-g.P))
    assert np.isclose(k.P[i],g.P[iP],rtol=1e-8,atol=1e-20)

    ig = np.argmin(np.abs(k.weightsn[i]-g.ng))
    assert np.isclose(k.weightsn[i],g.ng[ig],rtol=1e-8,atol=1e-20)
    return iT, iP, ig

# @nb.njit()
def read_helios_results(folder, sp, nw):

    for l in range(nw):

        print(l,end='\r')

        with open(folder+'/Out_'+sp+'_bin'+str(l).rjust(4,'0')+'.dat','r') as f:
            lines = f.readlines()
        lines = remove_duplicates(lines)

        k = KFile(len(lines))

        for i,line in enumerate(lines):
            tmp = line.split()
            k.weights[i] = float(tmp[0])
            k.xs[i] = float(tmp[1])
            k.T[i] = float(tmp[2])
            # atm * (1.01325 bar / 1 atm)
            k.P[i] = float(tmp[3])*(1.01325/1)
            k.weightsn[i] = float(tmp[4])

        if l == 0:
            g = KGrid(k)
            kcoeff = np.ones((len(g.ng),len(g.P),len(g.T),nw),order='F')*-1

        for i in range(len(k.xs)):
            iT, iP, ig = helper(k, g, i)
            kcoeff[ig,iP,iT,nw-l-1] = k.xs[i]
            
    assert np.all(np.min(kcoeff) >= 0)
            
    return g.g, g.P, g.T, kcoeff

def create_k_dataset(filename, sp, notes, weights, T, log10P, wavelengths, log10k):

    with h5py.File(filename, "w") as f:

        l = "50"
        dset = f.create_dataset("species",(),dtype="S"+l)
        dset[()] = ('{:'+l+'}').format(sp).encode()

        dset = f.create_dataset("notes",(),dtype="S1000")
        dset[()] = '{:1000}'.format(notes).encode()

        dset = f.create_dataset("weights", weights.shape, 'f')
        dset[:] = weights
        
        dset = f.create_dataset("T", T.shape, 'f')
        dset[:] = T
        
        dset = f.create_dataset("log10P", log10P.shape, 'f')
        dset[:] = log10P
        
        dset = f.create_dataset("wavelengths", wavelengths.shape, 'f')
        dset[:] = wavelengths

        dset = f.create_dataset("log10k", log10k.T.shape, 'f')
        dset[:] = log10k.T

def make_notes(commit, date, creator):

    notes = \
f"""These k-coefficients were created with HITEMP and HITRAN
line lists processed with HELIOS-k. HELIOS-k parameters and 
data used are given at this github repository and commit,

repository: https://github.com/Nicholaswogan/HELIOS-K
commit: {commit}

Date: {date}
Creator: {creator}

units:
wavelengths: [um]
log10P: [bar]
T: [K]
log10k: [cm^2/molecule], dimensions (len(weights),len(log10P),len(T),len(wavelengths))
"""
    assert len(notes) < 999
    return notes

def main():
    # Inputs
    species = ['C2H2','C2H6','CH4','CO','CO2','H2O','HCl','N2O','NH3','O2','O3','OCS','SO2']
    output_folder = 'kdistributions/'
    commit = None
    date = '9/12/25'
    creator = 'Nicholas Wogan'

    # Make folder if needed
    if not os.path.isdir(output_folder):
        os.mkdir('output_folder')

    # Make bins file
    sol_wavl = wogan_bins.wavl[:-2].copy()
    ir_wavl = wogan_bins.wavl[9:].copy()
    with h5py.File(os.path.join(output_folder,"bins.h5"), "w") as f:
        dset = f.create_dataset("sol_wavl", sol_wavl.shape, 'f')
        dset[:] = sol_wavl
        dset = f.create_dataset("ir_wavl", ir_wavl.shape, 'f')
        dset[:] = ir_wavl

    # Loop through all species
    nw = len(wogan_bins.wavnum)-1
    for i,sp in enumerate(species):

        # Get results
        g_value, P_grid, T_grid, kcoeff = read_helios_results('./', sp, nw)

        # checks
        assert np.all(np.isclose(P_grid,wogan_bins.P_grid,rtol=1e-10,atol=1e-20))
        # assert np.all(np.isclose(T_grid,bins.T_grid,rtol=1e-10,atol=1e-20))
        tmp = wogan_bins.weights_to_bins(wogan_bins.weights)
        g = (tmp[1:]+tmp[:-1])/2
        assert np.all(np.isclose(g, g_value))

        # clip
        kcoeff = np.clip(kcoeff, a_min=1.0e-60, a_max=np.inf)
        log10k = np.log10(kcoeff) # log10
        log10P = np.log10(P_grid)

        notes = make_notes(commit, date, creator)

        # create file
        outfilename = os.path.join(output_folder,sp+".h5")
        create_k_dataset(outfilename, sp, notes, wogan_bins.weights, T_grid, log10P, wogan_bins.wavl, log10k)

if __name__ == '__main__':
    main()