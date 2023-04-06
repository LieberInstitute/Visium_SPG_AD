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
img_channels = ['Abeta', 'DAPI', 'GFAP', 'Lipofuscin', 'MAP2', 'pTau']
default_channels = {'blue': 'DAPI', 'red': 'Abeta'}

sample_info_path = here(
    'raw-data', 'Visium_SPG_AD_ITG_MasterExcelSummarySheet.xlsx'
)

img_path = here(
    'processed-data', 'Images', 'VistoSeg', 'Capture_Areas', '{}.tif'
)
json_path = here(
    'processed-data', 'spaceranger', '{}', 'outs', 'spatial',
    'scalefactors_json.json'
)
tissue_path = here(
    'processed-data', 'spaceranger', '{}', 'outs', 'spatial',
    'tissue_positions_list.csv'
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

#   Subset both types of IDs to this sample only
sample_id_spaceranger = sample_info['spaceranger_id'].iloc[int(os.environ['SGE_TASK_ID']) - 1]
sample_id_image = sample_info['image_id'].iloc[int(os.environ['SGE_TASK_ID']) - 1]

#   Update paths for this sample ID
out_dir = Path(str(out_dir).format(sample_id_spaceranger))
json_path = Path(str(json_path).format(sample_id_spaceranger))
img_path = Path(str(img_path).format(sample_id_image))
tissue_path = Path(str(tissue_path).format(sample_id_spaceranger))

out_dir.mkdir(exist_ok = True)

#   All paths should exist
assert all([x.exists() for x in [out_dir, json_path, img_path, tissue_path]])

################################################################################
#   Read in spatial coordinates and scale-factors info
################################################################################

#   Read in the tissue positions file to get spatial coordinates. Index by
#   barcode
tissue_positions = pd.read_csv(
    tissue_path,
    header = None,
    names = ["in_tissue", "row", "col", "y", "x"], # Note the switch of x and y
    index_col = 0
)
tissue_positions.index.name = "barcode"

#   Read in the spaceranger JSON, ultimately to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

################################################################################
#   Use the Samui API to create the importable directory for this sample
################################################################################

this_sample = Sample(name = sample_id_spaceranger, path = out_dir)

this_sample.add_coords(
    tissue_positions, name = "coords", mPerPx = m_per_px, size = spot_diameter_m
)

#   Add the IF image for this sample
this_sample.add_image(
    tiff = img_path, channels = img_channels, scale = m_per_px,
    defaultChannels = default_channels
)

#this_sample.set_default_feature(group = "Genes", feature = "SNAP25")

this_sample.write()
