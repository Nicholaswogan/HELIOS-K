# Photochem Opacities with HELIOS-K

This repo is a modified version of [HELIOS-K](https://github.com/exoclime/HELIOS-K) to compute the k-tables for the Photochem code (links [here](https://github.com/Nicholaswogan/photochem) and [here](https://github.com/Nicholaswogan/photochem_clima_data)). The repo also generates [PICASO](https://github.com/natashabatalha/picaso) resampled opacity databases for computing planetary spectra. All of these opacities are targeted at rocky planets that are relatively cold (< 1000 K). Below are the instructions for running this code.

## Step 1: Compile HELIOS-K

You will need to compile HELIOS-K relevant code. This is accomplished by running `make heliosk` and `make hitran`. Please see the HELIOS-K docs for how to do this. The command `make heliosk` will require the nvcc cuda compiler and relevant cuda drivers. The command `make hitran` will require a C++ compiler.

## Step 2: Get cookies for downloading HITEMP data.

The code in the next steps needs to download HITEMP data. To accomplish this, you must make an account on https://hitran.org/. Next, login to the account. Next, you must export the cookies from your browser to a file named "cookies.txt" and place the file in the same directory as this README. I use this chrome extension to export all of my cookies: https://chromewebstore.google.com/detail/get-cookiestxt-locally/cclelndahbckbenkjhflpdbgdldlbecc?pli=1

## Step 3: Compute the opacities

First, create the following conda environment:

```sh
conda create -n test -c conda-forge python numba=0.60.0 h5py=3.11.0 matplotlib=3.9.2 pandas=2.2.3 scipy=1.11.4 requests=2.32.3 wget pip
conda activate test
python -m pip install hitran-api==1.2.2.2
```

Now run the main script. This script will first download all relevant line-by-line data from the internet (see `wogan_data/preprocess.py`). This takes a long time (~1 hour). Next, the script runs HELIOS-K for all molecules, pressures and temperatures. This takes about 48 hours on the GPU that I used.

```sh
python run_heliosk.py
```

Once the HELIOS-K calculation is complete, now you can generate the Photochem k-tables and the PICASO opacity DB with the following:

```sh
python make_photochem_ktable.py
python make_picaso_db.py
```


