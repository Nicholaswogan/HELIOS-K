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

def preprocess_default(molecule):

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
        cut=25.0, 
        qalphaL=0.0, 
        subLorentzianfile='-', 
        removePlinth=0
    )

def preprocess_C2H2():
    preprocess_default('C2H2')

def preprocess_C2H6():
    preprocess_default('C2H6')

def preprocess_CH4():

    molecule = 'CH4'
    molecule_zipname = "06_HITEMP2020.par.bz2"
    molecule_filename = "06_HITEMP2020.par"
    Species_Name = '06_hitrandata'

    data_folder = molecule
    data_name = 'hitrandata'

    global_ids = download_isotope_files(molecule, data_folder)
    molecule_id = get_molecule_id(molecule)

    # # Download
    folder = 'https://hitran.org/files/HITEMP/bzip2format/'
    cmd = 'wget --load-cookies='+os.path.join(MAIN_DIR,'cookies.txt')+' '+folder+molecule_zipname
    subprocess.run(cmd.split())
    os.rename(molecule_zipname, os.path.join(data_folder,molecule_zipname))

    # unzip the HITEMP data
    with open(os.path.join(data_folder,molecule_filename), 'wb') as new_file, bz2.BZ2File(os.path.join(data_folder,molecule_zipname), 'rb') as f:
        for data in iter(lambda : f.read(100 * 1024), b''):
            new_file.write(data)

    # rename the file
    os.rename(os.path.join(data_folder,molecule_filename),os.path.join(data_folder,Species_Name+'.par'))

    # Run preprocess
    run_hitran_preprocess(data_folder, data_name, molecule_id)

    # Param file
    write_param_file(
        filename=molecule+'/param.dat_'+molecule, 
        name=molecule, 
        Species_Name=Species_Name, 
        pathToData='wogan_data/'+molecule+'/', 
        cut=25.0, 
        qalphaL=0.0, 
        subLorentzianfile='-', 
        removePlinth=0
    )

def preprocess_CO():

    molecule = 'CO'
    molecule_zipname = "05_HITEMP2019.par.bz2"
    molecule_filename = "05_HITEMP2019.par"
    Species_Name = '05_hitrandata'

    data_folder = molecule
    data_name = 'hitrandata'

    global_ids = download_isotope_files(molecule, data_folder)
    molecule_id = get_molecule_id(molecule)

    # # Download
    folder = 'https://hitran.org/files/HITEMP/bzip2format/'
    cmd = 'wget --load-cookies='+os.path.join(MAIN_DIR,'cookies.txt')+' '+folder+molecule_zipname
    subprocess.run(cmd.split())
    os.rename(molecule_zipname, os.path.join(data_folder,molecule_zipname))

    # unzip the HITEMP data
    with open(os.path.join(data_folder,molecule_filename), 'wb') as new_file, bz2.BZ2File(os.path.join(data_folder,molecule_zipname), 'rb') as f:
        for data in iter(lambda : f.read(100 * 1024), b''):
            new_file.write(data)

    # rename the file
    os.rename(os.path.join(data_folder,molecule_filename),os.path.join(data_folder,Species_Name+'.par'))

    # Run preprocess
    run_hitran_preprocess(data_folder, data_name, molecule_id)

    # Param file
    write_param_file(
        filename=molecule+'/param.dat_'+molecule, 
        name=molecule, 
        Species_Name=Species_Name, 
        pathToData='wogan_data/'+molecule+'/', 
        cut=25.0, 
        qalphaL=1.0, 
        subLorentzianfile='-', 
        removePlinth=0
    )

def preprocess_CO2():

    molecule = 'CO2'
    data_folder = molecule
    data_name = 'hitrandata'
    Species_Name = '02_hitrandata'

    global_ids = download_isotope_files(molecule, data_folder)
    molecule_id = get_molecule_id(molecule)

    folder = 'https://hitran.org/files/HITEMP/HITEMP-2010/CO2_line_list/'
    files = """
    02_00000-00500_HITEMP2010.zip
    02_00500-00625_HITEMP2010.zip
    02_00625-00750_HITEMP2010.zip
    02_00750-01000_HITEMP2010.zip
    02_01000-01500_HITEMP2010.zip
    02_01500-02000_HITEMP2010.zip
    02_02000-02125_HITEMP2010.zip
    02_02125-02250_HITEMP2010.zip
    02_02250-02500_HITEMP2010.zip
    02_02500-03000_HITEMP2010.zip
    02_03000-03250_HITEMP2010.zip
    02_03250-03500_HITEMP2010.zip
    02_03500-03750_HITEMP2010.zip
    02_03750-04000_HITEMP2010.zip
    02_04000-04500_HITEMP2010.zip
    02_04500-05000_HITEMP2010.zip
    02_05000-05500_HITEMP2010.zip
    02_05500-06000_HITEMP2010.zip
    02_06000-06500_HITEMP2010.zip
    02_06500-12785_HITEMP2010.zip
    """.split()

    # Download HITEMP data
    for ffile in files:
        cmd = 'wget --load-cookies='+os.path.join(MAIN_DIR,'cookies.txt')+' '+folder+ffile
        subprocess.run(cmd.split())
        os.rename(ffile, os.path.join(data_folder,ffile))

        # unzip
        with zipfile.ZipFile(os.path.join(data_folder,ffile), 'r') as zip_ref:
            zip_ref.extractall(data_folder)

    # Rename par files
    for a in os.listdir(data_folder):
        if 'HITEMP2010' in a and a.endswith('.par'):
            aa = a.replace('HITEMP2010',data_name)
            tmp = aa.split('_')
            start = tmp[1].split('-')[0]
            start = start.rjust(5, '0')
            end = tmp[1].split('-')[1]
            end = end.rjust(5, '0')
            tmp[1] = start+'-'+end
            aa = "_".join(tmp)
            os.rename(os.path.join(data_folder,a),os.path.join(data_folder,aa))

    run_hitran_preprocess(data_folder, data_name, molecule_id)
    
    # Param file
    write_param_file(
        filename=molecule+'/param.dat_'+molecule, 
        name=molecule, 
        Species_Name=Species_Name, 
        pathToData='wogan_data/'+molecule+'/', 
        cut=500.0, 
        qalphaL=1.0, 
        subLorentzianfile='chiPH89__CO2-CO2.dat', 
        removePlinth=0
    )

def preprocess_H2O():

    molecule = 'H2O'
    data_folder = molecule
    data_name = 'hitrandata'
    Species_Name = '01_hitrandata'

    global_ids = download_isotope_files(molecule, data_folder)
    molecule_id = get_molecule_id(molecule)

    folder = 'https://hitran.org/files/HITEMP/HITEMP-2010/H2O_line_list/'
    files = """
    01_00000-00050_HITEMP2010.zip
    01_00050-00150_HITEMP2010.zip
    01_00150-00250_HITEMP2010.zip
    01_00250-00350_HITEMP2010.zip
    01_00350-00500_HITEMP2010.zip
    01_00500-00600_HITEMP2010.zip
    01_00600-00700_HITEMP2010.zip
    01_00700-00800_HITEMP2010.zip
    01_00800-00900_HITEMP2010.zip
    01_00900-01000_HITEMP2010.zip
    01_01000-01150_HITEMP2010.zip
    01_01150-01300_HITEMP2010.zip
    01_01300-01500_HITEMP2010.zip
    01_01500-01750_HITEMP2010.zip
    01_01750-02000_HITEMP2010.zip
    01_02000-02250_HITEMP2010.zip
    01_02250-02500_HITEMP2010.zip
    01_02500-02750_HITEMP2010.zip
    01_02750-03000_HITEMP2010.zip
    01_03000-03250_HITEMP2010.zip
    01_03250-03500_HITEMP2010.zip
    01_03500-04150_HITEMP2010.zip
    01_04150-04500_HITEMP2010.zip
    01_04500-05000_HITEMP2010.zip
    01_05000-05500_HITEMP2010.zip
    01_05500-06000_HITEMP2010.zip
    01_06000-06500_HITEMP2010.zip
    01_06500-07000_HITEMP2010.zip
    01_07000-07500_HITEMP2010.zip
    01_07500-08000_HITEMP2010.zip
    01_08000-08500_HITEMP2010.zip
    01_08500-09000_HITEMP2010.zip
    01_09000-11000_HITEMP2010.zip
    01_11000-30000_HITEMP2010.zip
    """.split()

    # Download HITEMP data
    for ffile in files:
        cmd = 'wget --load-cookies='+os.path.join(MAIN_DIR,'cookies.txt')+' '+folder+ffile
        subprocess.run(cmd.split())
        os.rename(ffile, os.path.join(data_folder,ffile))

        # unzip
        with zipfile.ZipFile(os.path.join(data_folder,ffile), 'r') as zip_ref:
            zip_ref.extractall(data_folder)

    # Rename par files
    for a in os.listdir(data_folder):
        if 'HITEMP2010' in a and a.endswith('.par'):
            aa = a.replace('HITEMP2010',data_name)
            tmp = aa.split('_')
            start = tmp[1].split('-')[0]
            start = start.rjust(5, '0')
            end = tmp[1].split('-')[1]
            end = end.rjust(5, '0')
            tmp[1] = start+'-'+end
            aa = "_".join(tmp)
            os.rename(os.path.join(data_folder,a),os.path.join(data_folder,aa))

    # Download
    hapi.fetch_by_ids(molecule, global_ids, 30_000, 42_000)
    os.remove(molecule+'.header')
    filename = f'{molecule_id:02d}_30000-42000_'+data_name+'.par'
    os.rename(molecule+'.data', os.path.join(data_folder,filename))

    run_hitran_preprocess(data_folder, data_name, molecule_id)
    
    # Param file
    write_param_file(
        filename=molecule+'/param.dat_'+molecule, 
        name=molecule, 
        Species_Name=Species_Name, 
        pathToData='wogan_data/'+molecule+'/', 
        cut=25.0, 
        qalphaL=0.0, 
        subLorentzianfile='-', 
        removePlinth=1
    )

def preprocess_HCl():
    preprocess_default('HCl')

def preprocess_N2O():

    molecule = 'N2O'
    molecule_zipname = "04_HITEMP2019.par.bz2"
    molecule_filename = "04_HITEMP2019.par"
    Species_Name = '04_hitrandata'

    data_folder = molecule
    data_name = 'hitrandata'

    global_ids = download_isotope_files(molecule, data_folder)
    molecule_id = get_molecule_id(molecule)

    # # Download
    folder = 'https://hitran.org/files/HITEMP/bzip2format/'
    cmd = 'wget --load-cookies='+os.path.join(MAIN_DIR,'cookies.txt')+' '+folder+molecule_zipname
    subprocess.run(cmd.split())
    os.rename(molecule_zipname, os.path.join(data_folder,molecule_zipname))

    # unzip the HITEMP data
    with open(os.path.join(data_folder,molecule_filename), 'wb') as new_file, bz2.BZ2File(os.path.join(data_folder,molecule_zipname), 'rb') as f:
        for data in iter(lambda : f.read(100 * 1024), b''):
            new_file.write(data)

    # rename the file
    os.rename(os.path.join(data_folder,molecule_filename),os.path.join(data_folder,Species_Name+'.par'))

    # Run preprocess
    run_hitran_preprocess(data_folder, data_name, molecule_id)

    # Param file
    write_param_file(
        filename=molecule+'/param.dat_'+molecule, 
        name=molecule, 
        Species_Name=Species_Name, 
        pathToData='wogan_data/'+molecule+'/', 
        cut=25.0, 
        qalphaL=0.0, 
        subLorentzianfile='-', 
        removePlinth=0
    )

def preprocess_NH3():
    preprocess_default('NH3')

def preprocess_O2():
    preprocess_default('O2')

def preprocess_O3():
    preprocess_default('O3')

def preprocess_OCS():
    preprocess_default('OCS')

def preprocess_SO2():
    preprocess_default('SO2')

def main():

    # # clean
    # species = ['C2H2','C2H6','CH4','CO','CO2','H2O','HCl','N2O','NH3','O2','O3','OCS','SO2']
    # for sp in species:
    #     for a in os.listdir(sp):
    #         aa = os.path.join(sp,a)
    #         if a != '.gitignore' and 'param.dat' not in a:
    #             if os.path.isdir(aa):
    #                 shutil.rmtree(aa)
    #             else:
    #                 os.remove(aa)
             
    preprocess_C2H2()
    preprocess_C2H6()
    preprocess_CH4()
    preprocess_CO()
    preprocess_CO2()
    preprocess_H2O()
    preprocess_HCl()
    preprocess_N2O()
    preprocess_NH3()
    preprocess_O2()
    preprocess_O3()
    preprocess_OCS()
    preprocess_SO2()

if __name__ == "__main__":
    main()