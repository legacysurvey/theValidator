from __future__ import print_function
import numpy as np
import fitsio
from glob import glob
import os
from astropy import units
from astropy.coordinates import SkyCoord

from astrometry.util.fits import fits_table, merge_tables
from tractor.brightness import NanoMaggies

def read_lines(fn_list):
    fin=open(fn_list,'r')
    lines=fin.readlines()
    fin.close()
    return np.sort(np.array( list(np.char.strip(lines)) ))

class CatalogueFuncs(object):
    '''funcs for using Dustins fits_table objects'''
    def stack(self,fn_list,textfile=True):
        '''concatenates fits tables'''
        if textfile: 
            fns=read_lines(fn_list)
        else:
            fns= fn_list
        if len(fns) < 1: raise ValueError('Error: fns=',fns)
        cats= []
        for i,fn in enumerate(fns):
            print('reading %d/%d' % (i+1,len(fns)))
            try: 
                tab= fits_table(fn) 
                cats.append( tab )
            except IOError:
                print('Fits file does not exist: %s' % fn)
        return merge_tables(cats, columns='fillzero')

    def set_mags(self,cat):
        '''adds columns to fits_table cat'''
        # Remove white spaces
        cat.set('type', np.char.strip(cat.get('type')))
        # AB mags
        # Two kinds, as observed (including dust) and instrinsic (dust removed)
        for whichmag in ['wdust','nodust']:
            # DECam
            for band in ['g','r','z','w1','w2']:
                if whichmag == 'wdust':
                    flux= cat.get('flux_%s' % band)
                    flux_ivar= cat.get('flux_ivar_%s' % band)
                elif whichmag == 'nodust':
                    flux= cat.get('flux_%s' % band)/cat.get('mw_transmission_%s' % band)
                    flux_ivar= cat.get('flux_ivar_%s' % band)*\
                                np.power(cat.get('mw_transmission_%s' % band),2)
                else: raise ValueError()
                mag,mag_err= NanoMaggies.fluxErrorsToMagErrors(flux, flux_ivar)
                cat.set('mag_%s_%s' % (whichmag,band),mag)
                cat.set('mag_ivar_%s_%s' % (whichmag,band),1./np.power(mag_err,2))
        # Instrinsic fluxes
        whichmag='nodust'
        # DECam
        for band in ['g','r','z','w1','w2']:
            flux= cat.get('flux_%s' % band)/cat.get('mw_transmission_%s' % band)
            flux_ivar= cat.get('flux_ivar_%s' % band)*\
                                    np.power(cat.get('mw_transmission_%s' % band),2)
        cat.set('flux_%s_%s' % (whichmag,band),flux)
        cat.set('flux_ivar_%s_%s' % (whichmag,band),flux_ivar)

    def set_mags_OldDataModel(self,cat):
        '''for tractor catalogues with columns like decam_flux.shape(many,6)
        adds columns to fits_table cat'''
        # Remove white spaces
        cat.set('type', np.char.strip(cat.get('type')))
        # AB mags
        # Two kinds, as observed (including dust) and instrinsic (dust removed)
        for whichmag in ['wdust','nodust']:
            # DECam
            shp=cat.get('decam_flux').shape
            mag,mag_err= np.zeros(shp),np.zeros(shp)
            for iband in range(shp[1]):
                if whichmag == 'wdust':
                    flux= cat.get('decam_flux')[:,iband]
                    flux_ivar= cat.get('decam_flux_ivar')[:,iband]
                elif whichmag == 'nodust':
                    flux= cat.get('decam_flux')[:,iband]/cat.get('decam_mw_transmission')[:,iband]
                    flux_ivar= cat.get('decam_flux_ivar')[:,iband]*\
                                np.power(cat.get('decam_mw_transmission')[:,iband],2)
                else: raise ValueError()
                mag[:,iband],mag_err[:,iband]=NanoMaggies.fluxErrorsToMagErrors(flux, flux_ivar)
            cat.set('decam_mag_%s' % whichmag,mag)
            cat.set('decam_mag_ivar_%s' % whichmag,1./np.power(mag_err,2))
            # WISE
            if 'wise_flux' in cat.get_columns(): 
                shp=cat.get('wise_flux').shape
                mag,mag_err= np.zeros(shp),np.zeros(shp)
                for iband in range(shp[1]):
                    if whichmag == 'wdust':
                        flux= cat.get('wise_flux')[:,iband]
                        flux_ivar= cat.get('wise_flux_ivar')[:,iband]
                    elif whichmag == 'nodust':
                        flux= cat.get('wise_flux')[:,iband]/cat.get('wise_mw_transmission')[:,iband]
                        flux_ivar= cat.get('wise_flux_ivar')[:,iband]*\
                                    np.power(cat.get('wise_mw_transmission')[:,iband],2)
                    mag[:,iband],mag_err[:,iband]=NanoMaggies.fluxErrorsToMagErrors(flux, flux_ivar)
                cat.set('wise_mag_%s' % whichmag,mag)
                cat.set('wise_mag_ivar_%s' % whichmag,1./np.power(mag_err,2))
        # Instrinsic fluxes
        whichmag='nodust'
        # DECam
        shp=cat.get('decam_flux').shape
        flux,flux_ivar= np.zeros(shp),np.zeros(shp)
        for iband in range(shp[1]):
            flux[:,iband]= cat.get('decam_flux')[:,iband]/cat.get('decam_mw_transmission')[:,iband]
            flux_ivar[:,iband]= cat.get('decam_flux_ivar')[:,iband]*\
                                    np.power(cat.get('decam_mw_transmission')[:,iband],2)
        cat.set('decam_flux_%s' % whichmag,flux)
        cat.set('decam_flux_ivar_%s' % whichmag,flux_ivar)
        # WISE 
        if 'wise_flux' in cat.get_columns(): 
            shp=cat.get('wise_flux').shape
            flux,flux_err= np.zeros(shp),np.zeros(shp)
            for iband in range(shp[1]):
                flux[:,iband]= cat.get('wise_flux')[:,iband]/cat.get('wise_mw_transmission')[:,iband]
                flux_ivar[:,iband]= cat.get('wise_flux_ivar')[:,iband]*\
                                        np.power(cat.get('wise_mw_transmission')[:,iband],2)
            cat.set('wise_flux_%s' % whichmag,flux)
            cat.set('wise_flux_ivar_%s' % whichmag,flux_ivar)

        
class Cuts4MatchedCats(object):
    '''Cuts for MATCHED cats only'''
    def __init__(self,matched1,matched2):
        self.psf1 = (matched1.get('type') == 'PSF')
        self.psf2 = (matched2.get('type') == 'PSF')
        # Band dependent
        bands='grz'
        self.good={}
        for band,iband in zip(bands,range(6)):
            self.good[band]= ((matched1.get('flux_ivar_%s' % band) > 0) *\
                              (matched2.get('flux_ivar_%s' % band) > 0))      
        

#     def get_mags(self,cat):
#         mag= cat.get('decam_flux')/cat.get('decam_mw_transmission')
#         return 22.5 -2.5*np.log10(mag)

#     def get_mags_ivar(self,cat):
#         return np.power(np.log(10.)/2.5*cat.get('decam_flux'), 2)* \
#                                     cat.get('decam_flux_ivar')

class Matcher(object):
    '''sphere matches two ra,dec lists,
    ref,obs are astronometry.net "fits_tables" objects
    '''
    def __init__(self):
        pass
    def match_within(self,ref,obs,dist=1./3600):
        '''Find obs ra,dec that are within dist of ref ra,dec
        default=1 arcsec'''
        # cat1 --> ref cat
        # cat2 --> cat matching to the ref cat
        cat1 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        cat2 = SkyCoord(ra=obs.get('ra')*units.degree, dec=obs.get('dec')*units.degree)
        idx, d2d, d3d = cat1.match_to_catalog_3d(cat2)
        b= np.array(d2d) <= dist
        # Return 4 index arrays for indices where matches and where missing
        imatch=dict(ref=[],obs=[])
        imiss=dict(ref=[],obs=[])
        # 
        imatch['ref']= np.arange(len(ref))[b]
        imatch['obs']= np.array(idx)[b]
        print("Matched: %d/%d objects" % (imatch['ref'].size,len(ref)))
        imiss['ref'] = np.delete(np.arange(len(ref)), imatch['ref'], axis=0)
        imiss['obs'] = np.delete(np.arange(len(obs)), imatch['obs'], axis=0)
        return imatch,imiss,d2d

    def nearest_neighbor(self,ref):
        '''return indices of nearest nearsest and distances to them'''
        cat1 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        cat2 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        idx, d2d, d3d = cat1.match_to_catalog_3d(cat2,nthneighbor=2)
        #b= np.array(d2d) <= within
        ref_nn= np.arange(len(ref))
        obs_nn= np.array(idx)
        dist= np.array(d2d)
        return ref_nn,obs_nn,dist

    def nearest_neighbors_within(self,ref,obs,within=1./3600,min_nn=1,max_nn=5):
        '''Find obs ra,dec that are within dist of ref ra,dec
        default=1 arcsec'''
        # cat1 --> ref cat
        # cat2 --> cat matching to the ref cat
        cat1 = SkyCoord(ra=ref.get('ra')*units.degree, dec=ref.get('dec')*units.degree)
        cat2 = SkyCoord(ra=obs.get('ra')*units.degree, dec=obs.get('dec')*units.degree)
        ref_nn,obs_nn,dist={},{},{}
        for nn in range(min_nn,max_nn+1):
            idx, d2d, d3d = cat1.match_to_catalog_3d(cat2,nthneighbor=nn)
            b= np.array(d2d) <= within
            ref_nn[str(nn)]= np.arange(len(ref))[b]
            obs_nn[str(nn)]= np.array(idx)[b]
            dist[str(nn)]= np.array(d2d)[b]
            print("within 1arcsec, nn=%d, %d/%d" % (nn,ref_nn[str(nn)].size,len(ref)))
        return ref_nn,obs_nn,dist

class TargetTruth(object):
    '''Build Target Truth catalogues, matching to DR3
    '''
    def __init__(self):
        self.truth_dir='/project/projectdirs/desi/target/analysis/truth'
        self.dr3_dir='/global/project/projectdirs/cosmo/data/legacysurvey/dr3'
        self.save_dir='/project/projectdirs/desi/users/burleigh/desi/target/analysis/truth'
  
    def in_region(self,ra,dec, rlo=0.,rhi=360.,dlo=0.,dhi=30.):
        '''ra,dec are numpy arrays,
        return indices of ra,dec where they are inside specified region'''
        if rlo < rhi:
            return (ra >= rlo) * (ra <= rhi) *\
                   (dec >= dlo) * (dec <= dhi)
        else: # RA wrap
            return np.logical_or(ra >= rlo, ra <= rhi) *\
                   (dec >= dlo) * (dec <= dhi)

    def bricks_in_region(self,rlo=0.,rhi=360.,dlo=0.,dhi=30.):
        print('Region: rlo=%.2f, rhi=%.2f, dlo=%.2f, dhi=%.2f' % (rlo,rhi,dlo,dhi))
        bricks=fits_table(os.path.join(self.dr3_dir,'survey-bricks.fits.gz'))
        i={}
        # Loop over 4 corners of each Brick
        for cnt,(ra,dec) in zip(range(1,5),[('ra1','dec1'),('ra1','dec2'),('ra2','dec1'),('ra2','dec2')]):
            i[str(cnt)]= self.in_region(bricks.get(ra),bricks.get(dec),\
                                        rlo=rlo,rhi=rhi,dlo=dlo,dhi=dhi)
            print('corner=%s, number bricks=%d' % (str(cnt),len(bricks.get('ra')[ i[str(cnt)] ])))
        i= np.any((i['1'],i['2'],i['3'],i['4']),axis=0)
        print('any corner, number bricks=%d' % len(bricks.get('ra')[i]))
        names= bricks.get('brickname')[i] 
        if not len(list(set(names))) == len(names):
            raise ValueError('Repeated brick names')
        return names

    def sweep_corner2radec(self,text='000m005'):
        ra,dec= float(text[:3]),float(text[4:])
        if text[3] == 'm': 
            dec*= -1
        return ra,dec

    def sweeps_in_region(self,rlo=0.,rhi=360.,dlo=0.,dhi=30.):
        print('Region: rlo=%.2f, rhi=%.2f, dlo=%.2f, dhi=%.2f' % (rlo,rhi,dlo,dhi))
        fns=glob(os.path.join(self.dr3_dir,'sweep/3.0/','sweep-*.fits'))
        fns=np.array(fns)
        assert(len(fns) > 0)
        # Loop over sweep regions
        # Ask if any corner of each region is in the region we are interested in
        keep=np.zeros(len(fns)).astype(bool)
        for cnt,fn in enumerate(fns):
            left,right= os.path.basename(fn).split('.')[0].split('-')[1:]
            ra1,dec1= self.sweep_corner2radec(text=left)
            ra2,dec2= self.sweep_corner2radec(text=right)
            # Loop over the 4 corners
            b=False
            for ra,dec in [(ra1,dec1),(ra1,dec2),(ra2,dec1),(ra2,dec2)]:
                b= np.logical_or(b, self.in_region(ra,dec,\
                                                   rlo=rlo,rhi=rhi,dlo=dlo,dhi=dhi)
                                 )
            if b:
                print(ra1,'-->',ra2,'  ',dec1,'-->',dec2)
                keep[cnt]= True
        if not len(fns[keep]) > 0:
            raise ValueError('Something amiss, no sweeps overlap region')
        print('%d/%d sweeps in the region' % (len(fns[keep]),len(fns)))
        return fns[keep] 

    def cosmos_zphot(self):
        # Data
        cosmos=fits_table(os.path.join(self.truth_dir,'cosmos-zphot.fits.gz'))
        # Bricks
        bnames= self.bricks_in_region(rlo=cosmos.get('ra').min(), rhi=cosmos.get('ra').max(),\
                                          dlo=cosmos.get('dec').min(),dhi=cosmos.get('dec').max())
        # Tractor Catalogues --> file list
        catlist= os.path.join(self.save_dir,'cosmos_dr3_bricks.txt')
        if not os.path.exists(catlist):
            fout=open(catlist,'w')
            for b in bnames:
                fn= os.path.join(self.dr3_dir,'tractor/%s/tractor-%s.fits' % (b[:3],b))
                fout.write('%s\n' % fn)
            fout.close()
            print('Wrote %s' % catlist)
        # Match
        fits_funcs= CatalogueFuncs()
        dr3=fits_funcs.stack(os.path.join(self.save_dir,'cosmos_dr3_bricks.txt'))
        mat=Matcher()
        imatch,imiss,d2d= mat.match_within(cosmos,dr3) #,dist=1./3600)
        cosmos.cut(imatch['ref'])
        dr3.cut(imatch['obs'])
        # Save
        cosmos.writeto(os.path.join(self.save_dir,'cosmos-zphot-dr3matched.fits'))
        dr3.writeto(os.path.join(self.save_dir,'dr3-cosmoszphotmatched.fits'))
        print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'cosmos-zphot-dr3matched.fits'),\
                                      os.path.join(self.save_dir,'dr3-cosmoszphotmatched.fits')))

    def vipers(self):
        # Data
        w1=fits_table(os.path.join(self.truth_dir,'vipers-w1.fits.gz'))
        w4=fits_table(os.path.join(self.truth_dir,'vipers-w4.fits.gz'))
        # Bricks
        for data in [w1,w4]:
            data.set('ra',data.get('alpha'))
            data.set('dec',data.get('delta'))
        bnames={}
        for data,key in zip([w1,w4],['w1','w4']):
            bnames[key]= self.bricks_in_region(rlo=data.get('ra').min(), rhi=data.get('ra').max(),\
                                              dlo=data.get('dec').min(),dhi=data.get('dec').max())
        bricks=np.array([])
        for key in bnames.keys():
            bricks=np.concatenate((bricks,bnames[key]))
        # Tractor Catalogues --> file list
        catlist= os.path.join(self.save_dir,'vipers_dr3_bricks.txt')
        if not os.path.exists(catlist):
            fout=open(catlist,'w')
            for b in bricks:
                fn= os.path.join(self.dr3_dir,'tractor/%s/tractor-%s.fits' % (b[:3],b))
                fout.write('%s\n' % fn)
            fout.close()
            print('Wrote %s' % catlist)
        # Merge w1,w4 for matching
        vipers= []
        for fn in [os.path.join(self.truth_dir,'vipers-w1.fits.gz'),\
                   os.path.join(self.truth_dir,'vipers-w4.fits.gz')]:
            vipers.append( fits_table(fn) )
        vipers= merge_tables(vipers, columns='fillzero')
        vipers.set('ra',data.get('alpha'))
        vipers.set('dec',data.get('delta'))
        # Match
        fits_funcs= CatalogueFuncs()
        dr3=fits_funcs.stack(os.path.join(self.save_dir,'vipers_dr3_bricks.txt'))
        mat=Matcher()
        imatch,imiss,d2d= mat.match_within(vipers,dr3) #,dist=1./3600)
        vipers.cut(imatch['ref'])
        dr3.cut(imatch['obs'])
        # Save
        vipers.writeto(os.path.join(self.save_dir,'vipersw1w4-dr3matched.fits'))
        dr3.writeto(os.path.join(self.save_dir,'dr3-vipersw1w4matched.fits'))
        print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'vipersw1w4-dr3matched.fits'),\
                                      os.path.join(self.save_dir,'dr3-vipersw1w4matched.fits')))

    def deep2(self):
        # Data
        deep2={}
        for key in ['1','2','3','4']:
            deep2[key]=fits_table('/project/projectdirs/desi/target/analysis/truth/deep2-field%s.fits.gz' % key)
        # Bricks
        bnames={}
        for key in deep2.keys():
            bnames[key]= self.bricks_in_region(rlo=deep2[key].get('ra').min(), rhi=deep2[key].get('ra').max(),\
                                              dlo=deep2[key].get('dec').min(),dhi=deep2[key].get('dec').max())
            print('Field=%s, Num Bricks=%d, Bricks:' % (key,len(bnames[key])), bnames[key])
        bricks=np.array([])
        for key in bnames.keys():
            bricks=np.concatenate((bricks,bnames[key]))
        # Tractor Catalogues --> file list
        catlist= os.path.join(self.save_dir,'deep2_dr3_bricks.txt')
        if not os.path.exists(catlist):
            fout=open(catlist,'w')
            for b in bricks:
                fn= os.path.join(self.dr3_dir,'tractor/%s/tractor-%s.fits' % (b[:3],b))
                fout.write('%s\n' % fn)
            fout.close()
            print('Wrote %s' % catlist)
        # Merge for matching
        dp2= [deep2['2'],deep2['3'],deep2['4']]
        dp2= merge_tables(dp2, columns='fillzero')
        # Match
        fits_funcs= CatalogueFuncs()
        dr3=fits_funcs.stack(os.path.join(self.save_dir,'deep2_dr3_bricks.txt'))
        mat=Matcher()
        imatch,imiss,d2d= mat.match_within(dp2,dr3) #,dist=1./3600)
        dp2.cut(imatch['ref'])
        dr3.cut(imatch['obs'])
        fits_funcs.set_extra_data(dr3)
        # Save
        dp2.writeto(os.path.join(self.save_dir,'deep2f234-dr3matched.fits'))
        dr3.writeto(os.path.join(self.save_dir,'dr3-deep2f234matched.fits'))
        print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'deep2f234-dr3matched.fits'),\
                                      os.path.join(self.save_dir,'dr3-deep2f234matched.fits')))

    def qso(self):
        # Christophe's qso catalogue
        qso=fits_table(os.path.join(self.save_dir,'CatalogQSO.fits.gz'))
        # DR7 and boss only
        keep={}
        for key in list(set(qso.get('source'))):
            print('%s: %d' % (key,len(qso[qso.get('source') == key])))
            keep[key.strip()]= qso.get('source') == key
        qso.cut( np.logical_or(keep['BOSS'],keep['SDSSDR7QSO']) )
        # Stripe 82
        rlo,rhi= 315., 45.
        dlo,dhi= -1.25, 1.25
        qso.cut( self.in_region(qso.get('ra'),qso.get('dec'), \
                                rlo=rlo,rhi=rhi,dlo=dlo,dhi=dhi) )
        # Bricks
        sweeps= self.sweeps_in_region(rlo=rlo, rhi=rhi,\
                                      dlo=dlo,dhi=dhi)
        print('sweeps= ',sweeps)
        sys.exit('early')
        # Tractor Catalogues --> file list
        #catlist= os.path.join(self.save_dir,'CatalogQSO_dr3_sweeps.txt')
        #if not os.path.exists(catlist):
        #    fout=open(catlist,'w')
        #    for b in bricks:
        #        fn= os.path.join(self.dr3_dir,'tractor/%s/tractor-%s.fits' % (b[:3],b))
        #        fout.write('%s\n' % fn)
        #    fout.close()
        #    print('Wrote %s' % catlist)
        # Match
        fits_funcs= CatalogueFuncs()
        #dr3=fits_funcs.stack(os.path.join(self.save_dir,'CatalogQSO_dr3_bricks.txt'))
        dr3=fits_funcs.stack(sweeps,textfile=False)
        mat=Matcher()
        imatch,imiss,d2d= mat.match_within(qso,dr3) #,dist=1./3600)
        qso.cut(imatch['ref'])
        dr3.cut(imatch['obs'])
        fits_funcs.set_extra_data(dr3)
        # Save
        qso.writeto(os.path.join(self.save_dir,'qso-dr3sweepmatched.fits'))
        dr3.writeto(os.path.join(self.save_dir,'dr3-qsosweepmatched.fits'))
        print('Wrote %s\nWrote %s' % (os.path.join(self.save_dir,'qso-dr3sweepmatched.fits'),\
                                      os.path.join(self.save_dir,'dr3-qsosweepmatched.fits')))





