import sys
import io
import requests
import urllib
import pandas
import numpy
from astropy.io import fits
from astropy.io import ascii

'''

Steps

1. Convert DESI QSO list to an IPAC Table:		qsoToTable()
2. Get ZTF objects in IPAC Table into xxx asynchronously: tableToObjects()
3. Get LCs of ZTF objects in matched Table

'''

radius_arcsec = 0.5
radius_degree = radius_arcsec/3600
desistart = 59197   # Dec 14 2020 first night of iron

uploadFile = "upload.tbl"
returnFile = "temp.tbl"

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



def tableToObjects():

	query = 'curl -o out.tbl -F filename=@upload.tbl -F catalog=ztf_objects_dr18 -F spatial=Upload -F uradius=0.5 -F uradunits=arcsec -F one_to_one=1 -F selcols=oid,ra,dec,ngoodobsrel -F outfmt=1 -Fconstraints=ngoodobsrel\>0 "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"'
	# query = """
	# 	SELECT TOP 1 oid, ra, dec, DISTANCE(POINT('ICRS',ra, dec), POINT('ICRS',TAP_UPLOAD.my_table.ra,TAP_UPLOAD.my_table.dec)) as dist
	# 	FROM ztf_objects_dr18
	# 	WHERE CONTAINS(POINT('ICRS',ra, dec), CIRCLE('ICRS',TAP_UPLOAD.my_table.ra,TAP_UPLOAD.my_table.dec,{0}))=1
	# 	ORDER BY dist
	# 	""".format(radius_degree)
	# print(query)
	# payload = {"QUERY": query, "FORMAT":"IPAC_TABLE", "UPLOAD": "my_table,param:table.tbl", "table.tbl":"@upload.tbl"}
	# params = urllib.parse.urlencode(payload, quote_via=urllib.parse.quote) 
	# # r = requests.get('https://irsa.ipac.caltech.edu/TAP/sync', params=payload)
	# q = "https://irsa.ipac.caltech.edu/TAP/sync?UPLOAD=my_table,param:table.tbl&table.tbl=@upload.tbl&FORMAT=IPAC_TABLE&QUERY=SELECT ra, dec FROM ztf_objects_dr18 WHERE CONTAINS(POINT('J2000',ra,dec), CIRCLE('J2000',TAP_UPLOAD.my_table.my_ra, TAP_UPLOAD.my_table.my_dec, 0.01)) =1"
	# r = requests.get(q)
	# print(r.url)
	# print(r.text)
	print(query)


def objectsToLCs():

	data = ascii.read(returnFile)
	df = data.to_pandas()
	df = df[df['oid'].notna()]

	# query to IPAC Helpdesk says only one position per query available
	for t in df.itertuples():
		if (pandas.notna(t.oid)):
			payload = {"ID": t.oid, "FORMAT": "CSV", "BAD_CATFLAGS_MASK": 32768}
			r = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves', params=payload)
			with open("{}.tbl".format(t.targetid_01), "w") as f:
				f.write(r.text)


def matchObjects():

	cols = {"oid": numpy.ulonglong, "ra": numpy.double, "dec": numpy.double, "clon": object, "clat": object, "dist": numpy.double,"angle": numpy.double}
	emptySeries = dict(cols)
	for key in emptySeries.keys():
		emptySeries[key] = None

	# coords = qsoCoords()
	qsodata = qsoData()
	frames=[]
	for d in qsodata:
		k = d['TARGETID'];  coord=(d['TARGET_RA'],d['TARGET_DEC'])
	# for k, coord in coords.items():
		payload = {"catalog": "ztf_objects_dr18", "spatial": "cone", "objstr": "{} {}".format(coord[0], coord[1]), \
			"radius": radius_arcsec, "radunits": "arcsec", "outfmt":1, "selcols": "oid,ra,dec"}
		r = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query', params=payload)
		out = r.text
		out=out.replace("|", "\\")
		df = pandas.read_csv(io.StringIO(out), sep="\s+", comment="\\", names=list(cols.keys()), dtype=cols)
		
		# if there are multiple choose the row closest to the unput ra/dec
		df = df[df.dist == df.dist.min()]

		# # if there are no rows make an empty one
		if df.shape[0] ==0:
			df.append(emptySeries, ignore_index=True)

		df['desi_targetid']=[k] 
		frames.append(df)
		print(df)
		wef

	dfs = pandas.concat(frames)
	print(dfs.shape[0])

if __name__ == '__main__':
    sys.exit(objectsToLCs())
