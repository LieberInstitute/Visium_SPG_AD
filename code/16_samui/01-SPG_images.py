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
img_channels = ['Lipofuscin', 'DAPI', 'GFAP', 'NeuN', 'OLIG2', 'TMEM119']
default_channels = {'blue': 'DAPI', 'red': 'NeuN'}

sample_info_path = here(
    'raw-data', 'Visium_IF_AD_ITG_MasterExcelSummarySheet.xlsx'
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
    .filter(["Br####", "Slide SN #", "Array #"])
    #   Clean up column names
    .rename(
        columns = {
            "Br####": "br_num",
            "Slide SN #": "sample_id",
            "Array #": "array_num"
        }
    )
)

#   Prepend "Br" tp brain number and make a string
sample_info['br_num'] = "Br" + sample_info['br_num'].astype(str)

#   Different forms of sample IDs appear to be used for spaceranger outputs
#   and raw images
sample_info = (sample_info
    .assign(
        spaceranger_id = sample_info['sample_id'].transform(lambda x: x.replace('-', '')) +
            '_' + sample_info['array_num'] + '_' + sample_info['br_num'],
        image_id = 'VIFAD1_' + sample_info['sample_id'] + '_' + sample_info['array_num']
    )
)
