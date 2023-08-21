import sys
import io
# import requests
import pandas
import numpy
# import urllib
# from astropy.io import fits
from astropy.io import ascii
import matplotlib.pyplot as plt
import sqlalchemy
# import pandas.io.sql as sqlio
import glob
import pyarrow.parquet as pq

conn_global = None
df_global = None
lc_dir='/Users/akim/project/QSO_Distance/data/'


# DESI database
def setup_conn():
    global conn_global
    with open('/global/cfs/cdirs/desi/science/td/secrets/desi_pg.txt') as f:
        file = f.read()
        db_name, db_user, db_pwd, db_host = file.split()
    conn_global = sqlalchemy.create_engine("postgresql+psycopg2://{}:{}@{}:{}/{}".format(db_user,db_pwd,"decatdb.lbl.gov",5432,"desidb"))

# ZTF information
def setup_df():
    global df_global
    returnFile = "../data/ztf.ztf_objects_dr18_17743.tbl"
    data = ascii.read(returnFile)
    df_global = data.to_pandas()
    df_global = df_global[df_global['oid'].notna()]  # prune out DESI objects not returned by ZTF


# The database query is slow.  Downloaded the entire ZTF light curve data
def targetidLC(targetid):

    # df2 = df.query('targetid_01==targetid')
    oid = 297103100005126
    tile=297
    ccd_id = 3
    q_id = 1

    files = glob.glob(lc_dir+'ztf_*{}_*_c*{}_q*{}_dr18.parquet'.format(tile,ccd_id,q_id))

    ans=[]
    for f in files:
        df = pq.read_table(f).to_pandas()
        ans.append(df.loc[df['objectid']==oid])

    result = pandas.concat(ans, ignore_index=True, copy=False)
    return result
        
    # # query to IPAC Helpdesk says only one position per query available
    # globalcounter=0
    # for t in df.itertuples():
    #     if (globalcounter % 1000 == 0): print(globalcounter)
        
    #     # get the dates this was observed by DESI

    #     curr = """SELECT f.mjd, f.night
    #         FROM iron.healpix_expfibermap f
    #         WHERE targetid={} AND fiberstatus=0""".format(t.targetid_01)
    #     try:
    #         desi_df = sqlio.read_sql_query(curr, conn)
    #     except:
    #         with open('error_db.txt','a') as f:
    #             f.write(str(t.targetid_01))
    #             wef
    #         continue
    #     desi_df = desi_df.groupby(['night']).mean().reset_index()

    #     for row in desi_df.itertuples():
    #         payload = {"ID": t.oid, "FORMAT": "CSV", "BAD_CATFLAGS_MASK": 32768,"TIME":"{} {}".format(row.mjd-lc_window, row.mjd+lc_window)}
    #         params = urllib.parse.urlencode(payload,quote_via=urllib.parse.quote)
    #         r = requests.get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves', params=params)
    #         if r.status_code != requests.codes.ok:
    #             with open('error_query.txt',"a") as f:
    #                 f.write(str(t.targetid_01)+" "+str(mjd['night']))
    #                 wef
    #             continue

    #         ans = pandas.read_csv(io.StringIO(r.text))
    #         if (not ans.empty):
    #             with open("{}/{}_{}.tbl".format(lcDir, t.targetid_01,row.night), "w") as f:
    #                 f.write(r.text)
    #     globalcounter=globalcounter+1


if __name__ == '__main__':
    sys.exit(targetidLC(None))
