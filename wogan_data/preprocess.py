import subprocess
import os
import shutil
import hapi
import bz2
import zipfile

MAIN_DIR = os.path.dirname(os.path.realpath(__file__))+'/../'

def get_global_ids(molecule):
    global_ids = []
    for key in hapi.ISO:
        if molecule == hapi.ISO[key][-1]:
            global_ids.append(hapi.ISO[key][0])
    return global_ids

def get_molecule_id(molecule):
    val = None
    for key in hapi.ISO:
        if molecule == hapi.ISO[key][-1]:
            val = key[0]
            break
    return val

def download_isotope_files(molecule, data_dir):
    global_ids = get_global_ids(molecule)

    for val in global_ids:
        # isotope stuff
        url = 'https://hitran.org/data/Q/q'
        cmd = 'wget '+url+str(val)+'.txt'
        subprocess.call(cmd.split())
        os.rename('q'+str(val)+'.txt', os.path.join(data_dir,'q'+str(val)+'.txt'))

    return global_ids

def hapi_download_and_process(molecule, numin=0.0, numax=1000000.0, data_folder='data/', data_name='hitrandata'):

    global_ids = download_isotope_files(molecule, data_folder)
    molecule_id = get_molecule_id(molecule)

    # Download
    hapi.fetch_by_ids(molecule, global_ids, numin, numax)
    os.remove(molecule+'.header')

    # rename
    Species_name = f'{molecule_id:02d}_'+data_name
    filename = Species_name+'.par'
    os.rename(molecule+'.data', os.path.join(data_folder,filename))

    run_hitran_preprocess(data_folder, data_name, molecule_id)

def run_hitran_preprocess(data_folder, data_name, molecule_id):
    
    # copy files to the main directory
    tmp_files = []
    for a in os.listdir(data_folder):
        if data_name in a and a.endswith('.par'):
            shutil.copy(os.path.join(data_folder,a), os.path.join(MAIN_DIR,a))
            tmp_files.append(os.path.join(MAIN_DIR,a))

    # preprocess the files
    cmd = "./hitran -M "+f'{molecule_id:02d}'+" -in "+data_name
    subprocess.run(cmd.split(), cwd=MAIN_DIR)

    # move processesed data files into data dir
    for a in os.listdir(MAIN_DIR):
        if data_name in a and ".bin" in a:
            os.rename(os.path.join(MAIN_DIR,a), os.path.join(data_folder,a))
        if data_name+".param" in a:
            os.rename(os.path.join(MAIN_DIR,a), os.path.join(data_folder,a))
    
    # delete the temporary files
    for tmp in tmp_files:
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass

def write_param_file(filename, name, Species_Name, pathToData, cut, qalphaL, subLorentzianfile, removePlinth):

    param_file = f"""name = {name}
T = TEMPERATURE
P = 1.0
PFile = wogan_data/PFile.txt
Species Name = {Species_Name}
SpeciesFile = -
ciaSystem = -
pathToData = {pathToData}
numin = -
numax = -
dnu = 0.1
Nnu per bin = 100000
cutMode = 0
cut = {cut} 
doResampling = 0
nC = 20
doTransmission = 0
nTr = 1000
dTr = 0.05
doStoreFullK = 2
pathToK =
doStoreSK = 2
nbins = -
binsFile = wogan_data/bins.txt
OutputEdgesFile = wogan_data/edges.txt
kmin = 0.0
qalphaL = {qalphaL}
gammaF = 1.0
doMean = 0
Units = 1
ReplaceFiles = 0
profile = 1
subLorentzianfile = {subLorentzianfile}
removePlinth = {removePlinth}
doTuning = 0
"""

    with open(filename,'w') as f:
        f.write(param_file)

def preprocess_default(molecule, cut=25.0, qalphaL=0.0, subLorentzianfile='-', removePlinth=0):

    data_folder = molecule
    data_name = 'hitrandata'
    molecule_id = get_molecule_id(molecule)
    Species_Name = f'{molecule_id:02d}_'+data_name
    hapi_download_and_process(
        molecule, 
        numin=0.0, 
        numax=1000000.0, 
        data_folder=data_folder, 
        data_name=data_name
    )
    write_param_file(
        filename=molecule+'/param.dat_'+molecule, 
        name=molecule, 
        Species_Name=Species_Name, 
        pathToData='wogan_data/'+molecule+'/', 
        cut=cut, 
        qalphaL=qalphaL, 
        subLorentzianfile=subLorentzianfile, 
        removePlinth=removePlinth
    )

def preprocess_CH4():
    preprocess_default('CH4')

def preprocess_CO():
    preprocess_default('CO', qalphaL=1.0)

def preprocess_CO2():
    preprocess_default('CO2', cut=500.0, qalphaL=1.0, subLorentzianfile='chiPH89__CO2-CO2.dat')

def preprocess_H2O():
    preprocess_default('H2O', removePlinth=1)

def preprocess_O2():
    preprocess_default('O2')

def preprocess_O3():
    preprocess_default('O3')

def main():

    # clean
    species = ['C2H2','C2H6','CH4','CO','CO2','H2O','HCl','N2O','NH3','O2','O3','OCS','SO2']
    for sp in species:
        for a in os.listdir(sp):
            aa = os.path.join(sp,a)
            if a != '.gitignore' and 'param.dat' not in a:
                if os.path.isdir(aa):
                    shutil.rmtree(aa)
                else:
                    os.remove(aa)
             
    preprocess_CH4()
    preprocess_CO()
    preprocess_CO2()
    preprocess_H2O()
    preprocess_O2()
    preprocess_O3()

if __name__ == "__main__":
    main()