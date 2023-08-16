import sys
import io
import requests
import urllib
import pandas
import numpy
from astropy.io import fits

radius_arcsec = 0.5
radius_degree = radius_arcsec/3600
desistart = 59197   # Dec 14 2020 first night of iron

def qsoData():
	fname = "../data/QSO_cat_iron_main_dark_healpix_v0.fits"
	hdul = fits.open(fname,memmap=True)

	# coords['Q1']=(298.0025, 29.87147)
	# coords['Q2']=(269.84158, 45.35492)
	return hdul['QSO_CAT'].data


def matchTable():
	query = """
		SELECT TOP 1 oid, ra, dec, DISTANCE(POINT('ICRS',ra, dec), POINT('ICRS', 298.0025,29.87147)) as dist
		FROM ztf_objects_dr18
		WHERE CONTAINS(POINT('ICRS',ra, dec), CIRCLE('ICRS',298.0025,29.87147,{0}))=1
		ORDER BY dist
		""".format(radius_degree)
	print(query)
	payload = {"QUERY": query, "outfmt":1}
	r = requests.get('https://irsa.ipac.caltech.edu/TAP/sync', params=payload)
	print(r.text)



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


def getLCs():
	coords = qsoCoords()
	# query to IPAC Helpdesk says only one position per query available

	for k, coord in coords.items():
		payload = {"POS": "CIRCLE {} {} {}".format(coord[0], coord[1], radius), \
			"BAD_CATFLAGS_MASK": 32768}

		params = urllib.parse.urlencode(payload, quote_via=urllib.parse.quote) 

		r = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves', params=params)

		text_file = open("{}_{}.xml".format(coord[0],coord[1]), "wt")
		n = text_file.write(r.text)
		text_file.close()


def main2():
	from astropy.io.votable import parse
	votable = parse("269.84158_45.35492.xml")
	table = votable.get_first_table()
	array = table.array

if __name__ == '__main__':
    sys.exit(matchTable())