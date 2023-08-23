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
'''

Steps

1. Convert DESI QSO list to an IPAC Table:      qsoToTable()
2. Get ZTF objects in IPAC Table into xxx asynchronously:  due to query limitations must be done on web interface.  tableToObjects()
3. Get LCs of ZTF objects in matched Table

'''

radius_arcsec = 0.5
radius_degree = radius_arcsec/3600
desistart = 59197   # Dec 14 2020 first night of iron

lc_window = 5.5 # +/- days

uploadFile = "../upload.tbl"
returnFile = "../data/ztf.ztf_objects_dr18_17743.tbl"
datesFile = "../data/dates.h5"

lcDir = "../data/lc/"

# cached variables
conn_global = None
ztfdf_global = None

# DESI database
def get_desi_conn():
    global conn_global
    if conn_global is None:
        with open('/global/cfs/cdirs/desi/science/td/secrets/desi_pg.txt') as f:
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
    lt = ["targetid={} ".format(t) for t in targetid]
    lt = " or ".join(lt)
    curr = """SELECT f.targetid, f.mjd, f.night
             FROM iron.healpix_expfibermap f
             WHERE {}""".format(lt)
    desi_df = sqlio.read_sql_query(curr, conn)
    desi_df = desi_df.groupby(['targetid', 'night']).mean().reset_index()
    return desi_df

def all_DESI_MJDs():
    df = get_ztf_df()
    ans = targetid_DESI_MJDs(df['targetid_01'].values)
    store = pandas.HDFStore(datesFile)
    store['df']=ans
    store.close()
   
# The database query is slow.  Downloaded the entire ZTF light curve data
def objectsToLCs():

    # DESI information
    fname = "/global/cfs/cdirs/desi/survey/catalogs/Y1/QSO/iron/QSO_cat_iron_main_dark_healpix_v0.fits"
    hdul = fits.open(fname,memmap=True)

    with open('/global/cfs/cdirs/desi/science/td/secrets/desi_pg.txt') as f:
        file = f.read()
        db_name, db_user, db_pwd, db_host = file.split()

    # ZTF information
    data = ascii.read(returnFile)
    df = data.to_pandas()
    df = df[df['oid'].notna()]  # prune out DESI objects not returned by ZTF

    # query to IPAC Helpdesk says only one position per query available
    globalcounter=0
    for t in df.itertuples():
        if (globalcounter % 1000 == 0): print(globalcounter)
        
        # get the dates this was observed by DESI
        conn = sqlalchemy.create_engine("postgresql+psycopg2://{}:{}@{}:{}/{}".format(db_user,db_pwd,"decatdb.lbl.gov",5432,"desidb"))

        curr = """SELECT f.mjd, f.night
            FROM iron.healpix_expfibermap f
            WHERE targetid={} AND fiberstatus=0""".format(t.targetid_01)
        try:
            desi_df = sqlio.read_sql_query(curr, conn)
        except:
            with open('error_db.txt','a') as f:
                f.write(str(t.targetid_01))
                wef
            continue
        desi_df = desi_df.groupby(['night']).mean().reset_index()

        for row in desi_df.itertuples():
           # hdul['QSO_CAT'].data['COADD_FIRSTMJD'][hdul['QSO_CAT'].data['TARGETID']== t.targetid_01][0]
            payload = {"ID": t.oid, "FORMAT": "CSV", "BAD_CATFLAGS_MASK": 32768,"TIME":"{} {}".format(row.mjd-lc_window, row.mjd+lc_window)}
            params = urllib.parse.urlencode(payload,quote_via=urllib.parse.quote)
            r = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves', params=params)
            if r.status_code != requests.codes.ok:
                with open('error_query.txt',"a") as f:
                    f.write(str(t.targetid_01)+" "+str(mjd['night']))
                    wef
                continue

            ans = pandas.read_csv(io.StringIO(r.text))
            if (not ans.empty):
                with open("{}/{}_{}.tbl".format(lcDir, t.targetid_01,row.night), "w") as f:
                    f.write(r.text)
        globalcounter=globalcounter+1


if __name__ == '__main__':
    sys.exit(all_DESI_MJDs())
