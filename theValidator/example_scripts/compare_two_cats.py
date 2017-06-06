from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import fitsio
import glob
import os
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from astropy import units
from astropy.coordinates import SkyCoord

from astrometry.util.fits import fits_table, merge_tables
from tractor.brightness import NanoMaggies


parser = ArgumentParser()
parser.add_argument('--ref_cat_list',action='store',required=True,help="text file listing refernce tractor catalogues, one per line")
parser.add_argument('--obs_cat_list',action='store',required=True,help="text file listing tractor catalogues to compare with, one per line")
parser.add_argument('--ref_name',action='store',default='REF',required=False,help="name for reference catalogues")
parser.add_argument('--obs_name',action='store',default='OBS',required=False,help="name for catalogues comparing to")
args = parser.parse_args()

import catalogues as cat
v2= cat.CatalogueFuncs().stack(args.ref_cat_list)
v3= cat.CatalogueFuncs().stack(args.obs_cat_list)

imatch,imiss,d2d= cat.Matcher().match_within(v2,v3) #,dist=1./3600)
#cat.CatalogueFuncs().set_mags_OldDataModel(v2)
#cat.CatalogueFuncs().set_mags_OldDataModel(v3)
cat.CatalogueFuncs().set_mags(v2)
cat.CatalogueFuncs().set_mags(v3)

#import plots_OldDataModel as plots
import plots
k=plots.Kaylans(v2,v3,imatch,
                ref_name=args.ref_name,obs_name=args.obs_name,savefig=True)
d=plots.Dustins(v2,v3,imatch,d2d,
                ref_name=args.ref_name,obs_name=args.obs_name,savefig=True)
