import sys
import io
import requests
import pandas
import numpy
import urllib
from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt

'''

Steps

1. Convert DESI QSO list to an IPAC Table:		qsoToTable()
2. Get ZTF objects in IPAC Table into xxx asynchronously: tableToObjects()
3. Get LCs of ZTF objects in matched Table

'''

radius_arcsec = 0.5
radius_degree = radius_arcsec/3600
desistart = 59197   # Dec 14 2020 first night of iron

lc_window = 5. # +/- days

uploadFile = "upload.tbl"
returnFile = "output.tbl"

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

    # ZTF information
	data = ascii.read(returnFile)
	df = data.to_pandas()
	df = df[df['oid'].notna()]	# prune out DESI objects not returned by ZTF

	# query to IPAC Helpdesk says only one position per query available
	counter = 0
	for t in df.itertuples():
        mjd_desi= hdul['QSO_CAT'].data['COADD_FIRSTMJD'][hdul['QSO_CAT'].data['TARGETID']== t.targetid_01][0]
        payload = {"ID": t.oid, "FORMAT": "CSV", "BAD_CATFLAGS_MASK": 32768,"TIME":"{} {}".format(mjd_desi-lc_window, mjd_desi+lc_window)}
        params = urllib.parse.urlencode(payload,quote_via=urllib.parse.quote)
        r = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves', params=params)
        ans = pandas.read_csv(io.StringIO(r.text))
        if (not ans.empty):
			with open("{}.tbl".format(t.targetid_01), "w") as f:
				f.write(r.text)
			count = count+1
		if count ==10:
			break

if __name__ == '__main__':
    sys.exit(objectsToLCs())
