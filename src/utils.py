import sys
import io
import pandas
import numpy
from astropy.io import ascii
import matplotlib.pyplot as plt
import sqlalchemy
import pandas.io.sql as sqlio
import glob
import pyarrow.parquet as pq
import warnings

conn_global = None
df_global = None
desi_global = None
lc_dir='/pscratch/sd/a/akim/ZTF/irsa.ipac.caltech.edu/data/ZTF/lc/lc_dr18/[01]/'

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
    if df_global is None: setup_df()
    df2 = df_global[df_global['targetid_01']==targetid]
    files = glob.glob(lc_dir+'field{0}/ztf_{0}_*_c{1}_q{2}_dr18.parquet'.format(df2.field.iloc[0].astype('str').zfill(6),df2.ccdid.iloc[0].astype('str').zfill(2),df2.qid.iloc[0]))
    ans=[]
    for f in files:
        df = pq.read_table(f).to_pandas()
        ans.append(df.loc[df['objectid']==df2.oid.iloc[0]])

    result = pandas.concat(ans, ignore_index=True, copy=False)
    return result

# get the dates this was observed by DESI
def targetidDESIMJDs(targetid):
    global desi_global
    if conn_global is None: setup_conn()
    curr = """SELECT f.mjd, f.night
             FROM iron.healpix_expfibermap f
             WHERE targetid={}""".format(targetid)
    desi_df = sqlio.read_sql_query(curr, conn_global)
    desi_df = desi_df.groupby(['night']).mean().reset_index()

    # with warnings.catch_warnings(record=True) as w:
    #     warnings.simplefilter("always")
    #     desi_df = desi_df.groupby(['night']).mean().reset_index()
    #     if len(w):
    #         print(desi_df['night'])
    #         print(targetid)
    #         wef

    return desi_df

# select on dates
# only good photometry kept based on 32768 bitmask
def trimLC(lc, dates, window):
    index = numpy.full(lc.nepochs, False)
    hmjd = lc['hmjd'].values[0]
    for date in dates.itertuples():
        index = numpy.logical_or(index, numpy.logical_and(numpy.abs(hmjd-date.mjd) < window/2,(lc['catflags'].values[0] & 32768)==0))
    return hmjd[index], lc['magerr'].values[0][index], lc['mag'].values[0][index]
 
    # # query to IPAC Helpdesk says only one position per query available
    # globalcounter=0
    # for t in df.itertuples():
    #     if (globalcounter % 1000 == 0): print(globalcounter)
        

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

def main():
    if df_global is None: setup_df()
    print("number of targetids {}".format(df_global.shape[0]))
    ans={1:0, 3:0, 5:0, 7:0, 9:0}
    count=0
    for row in df_global.itertuples():
        if (count % 1000 ==0): print(count, ans)
        # targetid = 39627322701128888
        lc = targetidLC(row.targetid_01)
        dates = targetidDESIMJDs(row.targetid_01)
        for k in ans.keys():
            dum = trimLC(lc,dates,k)
            if (len(dum[0]) !=0):
                ans[k] +=1
        count += 1

    print(ans)

if __name__ == '__main__':
    sys.exit(main())
