import sys
import requests
import urllib


radius = 1./3600	# 1 arcsec
desistart = 59197   # Dec 14 2020 first night of iron

def main():
	coords=[(298.0025, 29.87147), (269.84158, 45.35492)]

	coords=coords[0:2]

	# query to IPAC Helpdesk says only one position per query available
	for coord in coords:

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
    sys.exit(main())