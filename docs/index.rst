HELIOS-K
========


| HELIOS-K calculates opacity functions for planetary atmopheres by using opacity
  line lists from different databases. Before the opacity functions can be calculated,
  the line lists need to be downloaded and preprocessed into binary files that can be 
  read from HELIOS-K. 
| HELIOS-K provides tools to automatically download and preprocess files from the
  ExoMol, HITRAN, HITEMP, NIST, Kurucz and VALD3 databases. 
| HELIOS-K is running on GPUs and require a Nvidia GPU with compute capability of 3.0 or higher.


| Developed by Simon Grimm & Kevin Heng.
| University of Bern.

Setup
~~~~~

.. toctree::
  :maxdepth: 1
  :caption: Setup:

  helios_k/compilation.rst


Databases
~~~~~~~~~

.. toctree::
  :maxdepth: 1
  :caption: Databases:

  helios_k/param.rst
  helios_k/exomol.rst
  helios_k/exomolSuperlines.rst
  helios_k/hitran.rst
  helios_k/hitemp.rst
  helios_k/kurucz.rst
  helios_k/nist.rst
  helios_k/vald.rst

Running HELIOS-K
~~~~~~~~~~~~~~~~

.. toctree::
  :maxdepth: 1
  :caption: Running HELIOS-K:

  helios_k/heliosk.rst


OPTIONS
~~~~~~~

.. toctree::
  :maxdepth: 1
  :caption: Options:

  helios_k/profiles.rst
  helios_k/subLorentzian.rst
  helios_k/cut.rst
  helios_k/plinth.rst
  helios_k/bins.rst
  helios_k/resampling.rst
  helios_k/options.rst
  helios_k/transmission.rst
  helios_k/mean.rst


Output Files
~~~~~~~~~~~~

.. toctree::
  :maxdepth: 1
  :caption: Output Files:

  helios_k/output.rst

Atomic Opacities plots
~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
  :maxdepth: 1
  :caption: Atomic Opacitities Plots:

  helios_k/atoms.rst


Citations
~~~~~~~~~

If you make use of HELIOS-K in your work, please cite

`Grimm et Al. (2021) <https://ui.adsabs.harvard.edu/abs/2021ApJS..253...30G/abstract>`_:

.. code-block::

	@ARTICLE{2021ApJS..253...30G,
	       author = {{Grimm}, Simon L. and {Malik}, Matej and {Kitzmann}, Daniel and {Guzm{\'a}n-Mesa}, Andrea and {Hoeijmakers}, H. Jens and {Fisher}, Chloe and {Mendon{\c{c}}a}, Jo{\~a}o M. and {Yurchenko}, Sergey N. and {Tennyson}, Jonathan and {Alesina}, Fabien and {Buchschacher}, Nicolas and {Burnier}, Julien and {Segransan}, Damien and {Kurucz}, Robert L. and {Heng}, Kevin},
		title = "{HELIOS-K 2.0 Opacity Calculator and Open-source Opacity Database for Exoplanetary Atmospheres}",
	      journal = {\apjs},
	     keywords = {Exoplanet atmospheres, 487, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
		 year = 2021,
		month = mar,
	       volume = {253},
	       number = {1},
		  eid = {30},
		pages = {30},
		  doi = {10.3847/1538-4365/abd773},
	archivePrefix = {arXiv},
	       eprint = {2101.02005},
	 primaryClass = {astro-ph.EP},
	       adsurl = {https://ui.adsabs.harvard.edu/abs/2021ApJS..253...30G},
	      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
	}


or 

`Grimm & Heng (2015) <https://ui.adsabs.harvard.edu/abs/2015ApJ...808..182G/abstract>`_:

.. code-block::

    @ARTICLE{2015ApJ...808..182G,
           author = {{Grimm}, Simon L. and {Heng}, Kevin},
            title = "{HELIOS-K: An Ultrafast, Open-source Opacity Calculator for Radiative Transfer}",
          journal = {\apj},
         keywords = {methods: numerical, planets and satellites: atmospheres, radiative transfer, Astrophysics - Earth and Planetary Astrophysics, Physics - Atmospheric and Oceanic Physics},
             year = "2015",
            month = "Aug",
           volume = {808},
           number = {2},
              eid = {182},
            pages = {182},
              doi = {10.1088/0004-637X/808/2/182},
    archivePrefix = {arXiv},
           eprint = {1503.03806},
     primaryClass = {astro-ph.EP},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2015ApJ...808..182G},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

