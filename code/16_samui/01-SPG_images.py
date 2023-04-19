from pathlib import Path
from pyhere import here
import json
import os
import scanpy as sc

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes

spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot
img_channels = [
    'DAPI', 'Abeta', 'pTau', 'GFAP', 'MAP2', 'Lipofuscin', 'segmented_Abeta',
    'segmented_pTau'
]
default_channels = {'blue': 'DAPI', 'red': 'Abeta'}
default_gene = 'SNAP25'

#   Names of continuous features expected to be columns in the observation 
#   data (colData) of the AnnData
spe_cont_features = ['PpTau', 'PAbeta']

#   Diagnosis by brain number (not included in the sample_info sheet)
sample_dx = {
    'Br3854': 'AD', 'Br3873': 'AD', 'Br3880': 'AD', 'Br3874': 'control'
}

sample_info_path = here(
    'raw-data', 'Visium_IF_AD_ITG_MasterExcelSummarySheet.xlsx'
)

spe_path = here('processed-data', '16_samui', 'spe.h5ad')
img_path = here('processed-data', '16_samui', 'combined_tiffs', '{}.tif')
json_path = here(
    'processed-data', 'spaceranger', '{}', 'outs', 'spatial',
    'scalefactors_json.json'
)

out_dir = here('processed-data', '16_samui', '{}')

################################################################################
#   Read in sample info and clean
################################################################################

#   Read in sample info, subset to relevant columns, and clean
sample_info = (pd.read_excel(sample_info_path)
    .query('`Sequenced? ` == "Yes"')
    .filter(["Br####", "Slide SN #", "Array #", "Sample #"])
    #   Clean up column names
    .rename(
        columns = {
            "Br####": "br_num",
            "Slide SN #": "sample_id",
            "Array #": "array_num",
            "Sample #": "sample_num"
        }
    )
)

#   Prepend "Br" tp brain number and make a string
sample_info['br_num'] = "Br" + sample_info['br_num'].astype(str)

#   Add diagnosis using brain number
sample_info['diagnosis'] = sample_info['br_num'].replace(sample_dx)

#   Fix the experiment number column (use strings of integers)
sample_info['experiment_num'] = ((sample_info['sample_num'] - 1) // 4  + 1).astype(str)

#   Different forms of sample IDs appear to be used for spaceranger outputs
#   and raw images
sample_info = (sample_info
    .assign(
        spaceranger_id = sample_info['sample_id'].transform(lambda x: x.replace('-', '')) +
            '_' + sample_info['array_num'] + '_' + sample_info['br_num'],
        image_id = 'VIFAD' + sample_info['experiment_num'] + '_' + sample_info['sample_id'] + '_' + sample_info['array_num']
    )
)

#   Subset all types of IDs to this sample only
sample_id_spaceranger = sample_info['spaceranger_id'].iloc[int(os.environ['SGE_TASK_ID']) - 1]
sample_id_image = sample_info['image_id'].iloc[int(os.environ['SGE_TASK_ID']) - 1]
sample_id_samui = sample_id_spaceranger + '_' + sample_info['diagnosis'].iloc[int(os.environ['SGE_TASK_ID']) - 1]

#   Update paths for this sample ID
out_dir = Path(str(out_dir).format(sample_id_samui))
json_path = Path(str(json_path).format(sample_id_spaceranger))
img_path = Path(str(img_path).format(sample_id_image))

out_dir.mkdir(exist_ok = True)

#   All paths should exist
assert all([x.exists() for x in [out_dir, json_path, img_path]])

################################################################################
#   Read in scale-factors info
################################################################################

#   Read in the spaceranger JSON to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################

#   Read in AnnData and subset to this sample
spe = sc.read(spe_path)
spe = spe[spe.obs['sample_id'] == sample_id_spaceranger, :]
spe.obs.index.name = "barcode"

#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
gene_df = pd.DataFrame(
    spe.X.toarray(),
    index = spe.obs.index,
    columns = spe.var['gene_name']
)

#   Some gene symbols are actually duplicated. Just take the first column in
#   any duplicated cases
gene_df = gene_df.loc[: , ~gene_df.columns.duplicated()].copy()

#   Samui seems to break when using > ~ 5,000 genes. Take just the genes where
#   at least 10% of spots have nonzero counts
gene_df = gene_df.loc[:, np.sum(gene_df > 0, axis = 0) > (gene_df.shape[0] * 0.1)].copy()

assert default_gene in gene_df.columns, "Default gene not in AnnData"

################################################################################
#   Use the Samui API to create the importable directory for this sample
################################################################################

this_sample = Sample(name = sample_id_samui, path = out_dir)

this_sample.add_coords(
    spe.obsm['spatial'].rename(
        columns = {'pxl_col_in_fullres': 'x', 'pxl_row_in_fullres': 'y'}
    ),
    name = "coords", mPerPx = m_per_px, size = spot_diameter_m
)

#   Add the IF image for this sample
this_sample.add_image(
    tiff = img_path, channels = img_channels, scale = m_per_px,
    defaultChannels = default_channels
)

#   Add gene expression results (multiple columns) as a feature
this_sample.add_csv_feature(
    gene_df, name = "Genes", coordName = "coords", dataType = "quantitative"
)

#   Add additional requested observational columns (colData columns)
this_sample.add_csv_feature(
    spe.obs[spe_cont_features], name = "Spot Coverage", coordName = "coords",
    dataType = "quantitative"
)

this_sample.set_default_feature(group = "Genes", feature = default_gene)

this_sample.write()
