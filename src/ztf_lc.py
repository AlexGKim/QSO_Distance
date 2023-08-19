import sys
import io
import requests
import pandas
import numpy
import urllib
from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt
import psycopg2

'''

Steps

1. Convert DESI QSO list to an IPAC Table:      qsoToTable()
2. Get ZTF objects in IPAC Table into xxx asynchronously: tableToObjects()
3. Get LCs of ZTF objects in matched Table

'''

radius_arcsec = 0.5
radius_degree = radius_arcsec/3600
desistart = 59197   # Dec 14 2020 first night of iron

lc_window = 5. # +/- days

uploadFile = "../upload.tbl"
returnFile = "../data/output.tbl"
lcDir = "../data/lc/"

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


# Have to run the query through the web interface as the API does not accept huge numbers of rows
# https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd
def tableToObjects():

    query = 'curl -o out.tbl -F filename=@upload.tbl -F catalog=ztf_objects_dr18 -F spatial=Upload -F uradius=0.5 -F uradunits=arcsec -F one_to_one=1 -F selcols=oid,ra,dec,ngoodobsrel -F outfmt=1 -Fconstraints=ngoodobsrel\>0 "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"'



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
        conn = psycopg2.connect(dbname=db_name, user=db_user, password=db_pwd, host=db_host)
        cur = conn.cursor()
        curr = """SELECT f.mjd
            FROM iron.healpix_expfibermap f
            WHERE targetid={} AND fiberstatus=0""".format(t.targetid_01)
        cur.execute(curr)
        mjd_desi = cur.fetchall()
        for mjd in mjd_desi:
#            mjd_desi= hdul['QSO_CAT'].data['COADD_FIRSTMJD'][hdul['QSO_CAT'].data['TARGETID']== t.targetid_01][0]
            payload = {"ID": t.oid, "FORMAT": "CSV", "BAD_CATFLAGS_MASK": 32768,"TIME":"{} {}".format(mjd[0]-lc_window, mjd[0]+lc_window)}
            params = urllib.parse.urlencode(payload,quote_via=urllib.parse.quote)
            r = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves', params=params)
            ans = pandas.read_csv(io.StringIO(r.text))
            if (not ans.empty):
                with open("{}/{}_{}.tbl".format(lcDir, t.targetid_01,mjd[0]), "w") as f:
                    f.write(r.text)
        globalcounter=globalcounter+1


if __name__ == '__main__':
    sys.exit(objectsToLCs())
