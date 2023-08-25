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
3. QSO-matched ZTF object list (IPAC Table) : ZTF query through web interface tableToObjects()
4. DESI QSO MJD of observation (h5) : all_DESI_MJDs()
5. ZTF light curves (parquet) : Though overkill fastest to download all light curves rather than database query
3. Get LCs of ZTF objects in matched Table: all_DESI_MJDs

'''

radius_arcsec = 0.5
radius_degree = radius_arcsec/3600

lc_window = 5.5 # +/- days

secretsFile = "/global/cfs/cdirs/desi/science/td/secrets/desi_pg.txt"
uploadFile = "../upload.tbl"
returnFile = "../data/ztf.ztf_objects_dr18_17743.tbl"
datesFile = "../data/dates.h5"
targetidLCFile = "../data/targetidLC.h5"
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
        data = ascii.read(returnFile)
        ztfdf_global = data.to_pandas()
        ztfdf_global = ztfdf_global[ztfdf_global['oid'].notna()]  # prune out DESI objects not returned by ZTF
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
# https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd
def tableToObjects():
    query = 'curl -o output.tbl -F filename=@upload.tbl -F catalog=ztf_objects_dr18 -F spatial=Upload -F uradius=0.5 -F uradunits=arcsec -F one_to_one=1 -F selcols=oid,ra,dec,field,ccdid,qid,ngoodobsrel -F outfmt=1 -Fconstraints=ngoodobsrel\>0 "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"'
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


def targetidLC():
    df = get_ztf_df()

    ccds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    qs=[1,2,3,4]
    dirs = glob.glob(ztflc_dir+'field*')
    print("number of directories {}".format(len(dirs)))
    f = h5py.File(targetidLCFile, 'w')
    counter = 0
    for di in dirs:
        print(counter)
        for ccd in ccds:
            for q in qs:
                files = glob.glob(di+'/ztf_*_*_c{}_q{}_dr18.parquet'.format(str(ccd).zfill(2), q))
                if len(files) ==0: continue
                
                dum_=[]
                for file in files:
                    df_ = pq.read_table(file).to_pandas()
                    dum_.append(df_)
                
                lcdf = pandas.concat(dum_, ignore_index=True, copy=False)
                jdf = pandas.merge(lcdf , df, left_on='objectid', right_on='oid', how='inner')
                for r in jdf.itertuples():
                    grp = f.create_group(str(r.targetid_01))
                    grp.create_dataset('hmjd',data=r.hmjd)
                    grp.create_dataset('mag',data=r.mag)
                    grp.create_dataset('magerr',data=r.magerr)
                    grp.create_dataset('catflags',data=r.catflags)
        counter +=1
    f.close()


# The database query is slow.  Downloaded the entire ZTF light curve data
def targetidLC_(targetid):
    df = get_ztf_df()
    df2 = df[df['targetid_01']==targetid]
    files = glob.glob(ztflc_dir+'field{0}/ztf_{0}_*_c{1}_q{2}_dr18.parquet'.format(df2.field.iloc[0].astype('str').zfill(6),df2.ccdid.iloc[0].astype('str').zfill(2),df2.qid.iloc[0]))
    ans=[]
    for f in files:
        df = pq.read_table(f).to_pandas()
        ans.append(df.loc[df['objectid']==df2.oid.iloc[0]])

    result = pandas.concat(ans, ignore_index=True, copy=False)
    return result

# select on dates
# only good photometry kept based on 32768 bitmask
def trimLC(lc, dates, window):
    index = numpy.full(lc.nepochs, False)
    hmjd = lc.hmjd
    for date in dates.itertuples():
        index = numpy.logical_or(index, numpy.logical_and(numpy.abs(hmjd-date.mjd) < window/2,(lc.catflags & 32768)==0))
    return hmjd[index], lc.magerr[index], lc.mag[index]

def countQSOwithData():
    # object table from ZTF
    df = get_ztf_df()

    # dates table from DESI
    store = pandas.HDFStore(datesFile,mode='r')
    dates_df = store['df']

    # merge the ztf and desi dates
    # mdf = pandas.merge(df, dates_df, left_on='targetid_01', right_on='targetid', how='inner')

    print("number of targetids {}".format(df.shape[0]))
    ans={1:0, 3:0, 5:0, 7:0, 9:0}
    count=0

    ccds = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
    qs=[1,2,3,4]
    dirs = glob.glob(ztflc_dir+'field*')
    print("number of directories {}".format(len(dirs)))
    for di in dirs:
        for ccd in ccds:
            for q in qs:
                files = glob.glob(di+'/ztf_*_*_c{}_q{}_dr18.parquet'.format(str(ccd).zfill(2), q))
                if len(files) ==0: continue
                
                dum_=[]
                for file in files:
                    df_ = pq.read_table(file).to_pandas()
                    dum_.append(df_)
                
                lcdf = pandas.concat(dum_, ignore_index=True, copy=False)
                jdf = pandas.merge(lcdf , df, left_on='objectid', right_on='oid', how='inner')

                # to do merge jdf with dates_df based on targetid

                for r in jdf.itertuples():
                    dates = dates_df[dates_df['targetid'] == r.targetid_01]
                    for k in ans.keys():
                        dum = trimLC(r, dates ,k)
                        if (len(dum[0]) !=0):
                            ans[k] +=1
        print(count, ans)
        count += 1

if __name__ == '__main__':
    sys.exit(targetidLC())
