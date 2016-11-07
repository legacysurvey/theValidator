# theValidator
python scripts for validating Legacy Survey Data

python code that only relies on anaconda installable packages (matplotlib, numpy, scipy, astropy, ...). Everything can be run inside a NERSC ipynb (jupyter.nersc.gov) which has access to /project and $HOME. A few functions are imported from Dustins astrometry and tractor, but no installation of those packages is needed just appending to your PYTHONPATH. thon notebooks via NERSCs Jupyter Hub jupyter.nersc.gov, that can read and analyze [Legacy Survey](http://legacysurvey.org/) data products, such as Tractor Catalogues, concantenating them, matching two sets of catalogues, etc. 

Instructions
============

- login to edison or cori, cd $HOME, git clone https://github.com/legacysurvey/theValidator.git
- open a browser, login to jupyter.nersc.gov, open the ipynb: $HOME/theValidator/compare_two_cats.ipynb
- follow the instructions

Features
========

- Anyone -- if you have a NERSC account, you can run all the validation code 
- Tested -- reproduces results and plots made by different people: Dustin, Enrique, Kaylan
- Fast -- at least for python, memory efficient using Dustins fits_table objects which uses fitsio

Current
========

- cosmos tests -- the plots made by Dustin, Enrique
- DECals vs. MzLS/BASS -- the comparison plots made by Kaylan 

Future
========

- healpy -- healipx maps of depth, seeing, airmass, etc (Ashely has already done some of this separately) 
- psycopg2 -- analysis using the desi PostgreSQL DB at NERSC
- pre-tractor validation
- post-tractor validation 
