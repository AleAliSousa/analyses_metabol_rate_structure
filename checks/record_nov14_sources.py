"""Record checksums and download locators for the source files used on Nov 14."""
import csv
import hashlib
from pathlib import Path
root=Path(__file__).resolve().parents[1]
rows=[]
def add(path,url,member='',kind='file'):
    p=root/path
    rows.append(dict(local_path=path,url=url,archive_member=member,container=kind,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest(),retrieved_date='2026-11-14'))
add('data_raw/Tushaus_2023/Dataset_EV1.xlsx','https://static-content.springer.com/esm/art%3A10.15252%2Fembj.2023114665/MediaObjects/44318_2023_BFEMBJ2023114665_MOESM3_ESM.xlsx')
for name in ['Bmax.mean.MNI152.sm5.nii.gz','Bmax.sd.MNI152.sm5.nii.gz','Bmax.cov.MNI152.sm5.nii.gz','regional_values.csv']:
    add('data_raw/Johansen_2024_atlas/'+name,'https://xtra.nru.dk/SV2A-atlas/data/SV2A_atlas.zip','SV2A_atlas/'+name,'zip')
for name in ['AAL3v1_1mm.nii.gz','AAL3v1_1mm.xml','AAL3_UserGuide_April2024.pdf']:
    add('data_raw/anatomical_atlases/AAL3v2/'+name,'https://www.gin.cnrs.fr/wp-content/uploads/AAL3v2_for_SPM12.tar.gz','AAL3/'+name,'tar.gz')
add('data_raw/anatomical_atlases/CIT168/STH_MNI2009c.nii.gz','https://neurovault.org/media/images/3145/CIT168toMNI152_prob_atlas_bilat_1mm__%28volume%2016%29.nii.gz')
for name,file_id in [('massp2_maxlabel_MNI2009b.nii.gz','50077914'),('labels.txt','58258639')]:
    add('data_raw/anatomical_atlases/MASSP2/'+name,'https://ndownloader.figshare.com/files/'+file_id)
for name in ['atl-Anatom_space-MNI_dseg.nii','atl-Anatom.tsv']:
    add('data_raw/anatomical_atlases/Diedrichsen/'+name,'https://raw.githubusercontent.com/DiedrichsenLab/cerebellar_atlases/master/Diedrichsen_2009/'+name)
with (root/'metadata/nov14_source_downloads.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(f'Recorded {len(rows)} source files.')
