

import os
import time
import pandas
import pysradb
import datetime

GSE_ID = "GSE152991"


# Initialize pysradb interface
db = pysradb.sraweb.SRAweb()

# Step 1: Get GSMs for GSE
print(f"Retrieving GSMs for {GSE_ID}...")
gsm_df = db.gse_to_gsm(GSE_ID)
print(gsm_df.head())
print('done')

# Step 2: Get SRRs for each GSM
gsm_ids = gsm_df["experiment_alias"].tolist()
print(f"Mapping {len(gsm_ids)} GSMs to SRRs...")
srr_df = db.gsm_to_srr(gsm_ids)
print(srr_df.head())
print('done')


srr_df.to_csv('info.tsv', sep='\t')



    

