import os
import sys
import io
import requests
import pandas
import numpy
import urllib
from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt
import sqlalchemy
import pandas.io.sql as sqlio
import glob
import pyarrow.parquet as pq
import h5py

'''

Products

1. DESI QSO list (fits)
2. DESI QSO list (IPAC Table) : qsoToTable()
3. QSO-matched ZTF object list (IPAC Table) : ZTF query through web interface tableToObjects() run separately for each filter
4. DESI QSO MJD of observation (hdf5) : all_DESI_MJDs()
5. ZTF light curves (parquet) : Though overkill fastest to download all light curves rather than database query
6. DESI QSO light curves ([field].hdf5) : Save subset of the light curves to separate file targetidLC()
7. DESI QSO light crves (targetidLCFile.hdf5) : One master hdf5 file with links to other files and a phone book

'''

radius_arcsec = 0.5
radius_degree = radius_arcsec/3600

lc_window = 5.5 # +/- days

secretsFile = "/global/cfs/cdirs/desi/science/td/secrets/desi_pg.txt"
uploadFile = "../data/upload.tbl"
#returnFile = "../data/ztf.ztf_objects_dr18_17743.tbl"
returnFiles = ["../data/ztf.ztf_objects_dr18_4140.tbl","../data/ztf.ztf_objects_dr18_10711.tbl","../data/ztf.ztf_objects_dr18_31886.tbl"]
datesFile = "../data/dates.hdf5"
lcdir = "../data/lc/"
targetidLCFile = "../data/targetidLC.hdf5"
# targetidLCFile = "../data/temp.h5"
ztflc_dir='/pscratch/sd/a/akim/ZTF/irsa.ipac.caltech.edu/data/ZTF/lc/lc_dr18/[01]/'

lcDir = "../data/lc/" # depracated

# cached variables
conn_global = None
ztfdf_global = None

# DESI database
def get_desi_conn():
    global conn_global
    if conn_global is None:
        with open(secretsFile) as f:
            file = f.read()
            db_name, db_user, db_pwd, db_host = file.split()
        conn_global = sqlalchemy.create_engine("postgresql+psycopg2://{}:{}@{}:{}/{}".format(db_user,db_pwd,"decatdb.lbl.gov",5432,"desidb"))
    return conn_global

# ZTF information
def get_ztf_df():
    global ztfdf_global
    if ztfdf_global is None:
        dum=[]
        for returnFile in returnFiles:
            data = ascii.read(returnFile)
            data = data.to_pandas()
            data=data[data['oid'].notna()] # prune out DESI objects not returned by ZTF 
            dum.append(data)

        ztfdf_global = pandas.concat(dum, ignore_index=True, copy=False)
    return ztfdf_global

# Convert DESI fits file to IPAC Table
# IRSA database query accepts IPAC Table
# https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html
def qsoToTable():
    fname = "/global/cfs/cdirs/desi/survey/catalogs/Y1/QSO/iron/QSO_cat_iron_main_dark_healpix_v0.fits"
    hdul = fits.open(fname,memmap=True)
    head = """|    ra      |     dec     |     targetid      |
|  double    |    double   |      long         |
|   deg      |    deg      |                   |
|   null     |    null     |      null         |
"""
    with open(uploadFile, 'w') as f:
        f.write(head)
        for d in hdul['QSO_CAT'].data:
            f.write("{:12.7f} {:13.7f} {:19d}\n".format(d['TARGET_RA'],d['TARGET_DEC'], d['TARGETID']))
    # # coords['Q1']=(298.0025, 29.87147)
    # # coords['Q2']=(269.84158, 45.35492)
    # return hdul['QSO_CAT'].data


# Query the ZTF database
# Have to run the query through the web interface as the API does not accept huge numbers of rows
# Have to run this 3 times for each of the 3 filters
# https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd
def tableToObjects():
    query = 'curl -o output.tbl -F filename=@upload.tbl -F catalog=ztf_objects_dr18 -F spatial=Upload -F uradius=0.5 -F uradunits=arcsec -F one_to_one=1 -F selcols=oid,ra,dec,field,ccdid,qid,filtercode,ngoodobsrel -F outfmt=1 -F constraints=ngoodobsrel\>0 and filtercode=quote zgquote "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"'
    return query

# get the dates this was observed by DESI
def targetid_DESI_MJDs(targetid):
    conn = get_desi_conn()
    lt = ", ".join(targetid.astype('str'))
    # curr = """SELECT f.targetid, f.mjd, f.night
    #          FROM iron.healpix_expfibermap f
    #          WHERE {}""".format(lt)
    curr = """SELECT f.targetid, f.mjd, f.night
             FROM iron.healpix_expfibermap f
             INNER JOIN
             unnest(ARRAY[{}]) as tid
             ON tid=f.targetid""".format(lt)
    desi_df = sqlio.read_sql_query(curr, conn)
    desi_df = desi_df.groupby(['targetid', 'night']).mean().reset_index()
    return desi_df

def all_DESI_MJDs():
    df = get_ztf_df()
    ans = targetid_DESI_MJDs(df['targetid_01'].values)
    store = pandas.HDFStore(datesFile,mode='w')
    store['df']=ans
    store.close()

# make per-field lc files
def targetidLC(overwrite=False):
    if not overwrite:
        print("Doing nothing")
        sys.exit(-1)

    df = get_ztf_df()

    ccds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    qs=[1,2,3,4]

    filteridtostr = {1:'g', 2:'r', 3:'i'}

    # get all of the available fields
    dirs = glob.glob(ztflc_dir+'field*')
    print("number of directories/fields {}".format(len(dirs)))

    counter = 0
    # loop for each field
    # for completeness make an hdf5 file for every field even if it is empty

    for di in dirs:
        dibase  = os.path.basename(di)
        dibase=dibase[5:]
        f_sub = h5py.File(lcdir+"{}.hdf5".format(dibase), 'w')
        print(counter, dibase)

        # try every possible patch of sky
        for ccd in ccds:
            for q in qs:

                files = glob.glob(di+'/ztf_*_*_c{}_q{}_dr18.parquet'.format(str(ccd).zfill(2), q))
                if len(files) ==0: continue # if this patch has no lightcurves skip
               
                # merge all the objects in the different bands into one table
                dum_=[]
                for file in files:
                    df_ = pq.read_table(file).to_pandas()
                    dum_.append(df_)

                lcdf = pandas.concat(dum_, ignore_index=True, copy=False)

                # connect light curves with DESI targetid through an inner join for efficiency
                jdf = pandas.merge(lcdf , df, left_on='objectid', right_on='oid', how='inner')

                # gather LC information on each DESI targetid
                for utid in jdf['targetid_01'].unique():
                    subjdf_ = jdf[jdf['targetid_01']==utid]
                    if str(utid) not in f_sub:
                        grp = f_sub.create_group(str(utid))
                    else:
                        grp = f_sub[str(utid)]
                    
                    # collect for each filter
                    for fid_ in range(1,4):
                        subsubjdf_ = subjdf_[subjdf_['filterid']==fid_]
                        if filteridtostr[fid_] not in grp:
                         grp2 = grp.create_group(filteridtostr[fid_])
                         grp2.attrs.create('nepochs',0)
                        else:
                         grp2 = grp[filteridtostr[fid_]]
                        
                        nep=grp2.attrs.get('nepochs')

                         #grp2 = grp.create_group(filteridtostr[fid_])

                        if nep ==0 and subsubjdf_.shape[0] > 0:
                            #grp2.attrs.create('nepochs',0)
                        #else:
                            for r in subsubjdf_.itertuples():
                                # grp2.attrs.create('nepochs',r.nepochs)
                                grp2.attrs.modify("nepochs", r.nepochs)
                                grp2.create_dataset('hmjd',data=r.hmjd)
                                grp2.create_dataset('mag',data=r.mag)
                                grp2.create_dataset('magerr',data=r.magerr)
                                grp2.create_dataset('catflags',data=r.catflags)
        
        f_sub.close()
        counter +=1

# make uber list file that links to everything else
def linktargetidLCFile():
    f = h5py.File(targetidLCFile, 'w')
    grp=f.create_group('fields')

    dirs= glob.glob(lcdir+"0*.hdf5")

    # dirs = glob.glob(ztflc_dir+'field*')
    tid=[]
    fid=[]
    for di in dirs:
        dibase  = os.path.basename(di)
        dibase=dibase[:-5]
        grp[dibase] = h5py.ExternalLink(lcdir+"{}.hdf5".format(dibase), "/")
        for t in grp[dibase].keys():
            tid.append(t)
            fid.append(dibase)

    tid = numpy.array(tid,dtype='uint')
    fid = numpy.array(fid,dtype='uint')
    asort = numpy.argsort(tid)
    tid=tid[asort]
    fid=tid[asort]
    pgrp=f.create_group('pbook')
    pgrp.create_dataset('targetids', data= tid)
    pgrp.create_dataset('fieldids',data = fid)
    f.close()

# select on dates
# only good photometry kept based on 32768 bitmask
def trimLC(hmjd, mag, magerr, catflags, dates, window):
    index = numpy.full(len(hmjd), False)
    for date in dates.itertuples():
        index = numpy.logical_or(index, numpy.logical_and(numpy.abs(hmjd-date.mjd) < window/2,(catflags & 32768)==0))
    return hmjd[index], magerr[index], mag[index]

def countQSOwithData():
    # lcs 
    lcstore = h5py.File(targetidLCFile,'r')
    fields = lcstore['fields']

    # dates table from DESI
    store = pandas.HDFStore(datesFile,mode='r')
    dates_df = store['df']

    # number of target ids
    ntid=0
    for k in fields.keys():
        ntid += len(fields[k].keys())
    print("number of targetids {}".format(ntid))

    perfilt = {'g': 0, 'r': 0, 'i': 0}
    ans={1:perfilt, 3:dict(perfilt), 5:dict(perfilt), 7:dict(perfilt), 9:dict(perfilt)}

    count=0
    for k in fields.keys():  # loop over fields
          for t in fields[k].keys():  # loop over target ids
                dates = dates_df[dates_df['targetid']==int(t)]
                for b in perfilt.keys():  # loop over filters
                    for window in ans.keys():  # loop over time windows
                         if fields[k][t][b].attrs.get('nepochs') ==0: continue
                         # print(fields[k][t][b].attrs.get('nepochs'))
                         dum = trimLC(numpy.array(fields[k][t][b]['hmjd']), numpy.array(fields[k][t][b]['mag']), numpy.array(fields[k][t][b]['magerr']),numpy.array(fields[k][t][b]['catflags']), dates , window)
                         if (len(dum[0]) !=0):
                             ans[window][b] +=1
          print(count, ans)
          count += 1



def main():
    f = h5py.File(targetidLCFile, 'r')
    print(len(f.keys()))
    tidwithlc=0
    for k in f['/targetids'].keys():
        tidwithlc += len(f[k].keys())
    print(tidwithlc)

if __name__ == '__main__':
    # all_DESI_MJDs()
    # targetidLC(overwrite=True)
    linktargetidLCFile()
    # sys.exit(countQSOwithData())
    # sys.exit(main())
