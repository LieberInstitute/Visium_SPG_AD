import tifffile
from PIL import Image
from pyhere import here
from pathlib import Path
import pandas as pd
import numpy as np
import session_info
import os

#   Images are too large to open with default PIL settings
Image.MAX_IMAGE_PIXELS = None

#   Diagnosis by brain number (not included in the sample_info sheet)
sample_dx = {
    'Br3854': 'AD', 'Br3873': 'AD', 'Br3880': 'AD', 'Br3874': 'control'
}

sample_info_path = here(
    'raw-data', 'Visium_IF_AD_ITG_MasterExcelSummarySheet.xlsx'
)

raw_img_path = here(
    'processed-data', 'Images', 'VistoSeg', 'Capture_Areas', '{}.tif'
)

segmented_img_path = here(
    'processed-data', 'Images', 'VistoSeg', 'Segmentations',
    '{}_segmentation.tif'
)

out_path = here('processed-data', '16_samui', 'combined_tiffs', '{}.tif')

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

################################################################################
#   Subset variables to this sample and ensure paths exist
################################################################################

#   Subset to this sample only
sample_id_image = sample_info['image_id'].iloc[int(os.environ['SGE_TASK_ID']) - 1]
# sample_id_image = sample_info['image_id'].iloc[5] # for testing

#   Update paths for this sample ID
raw_img_path = Path(str(raw_img_path).format(sample_id_image))
segmented_img_path = Path(str(segmented_img_path).format(sample_id_image))
out_path = Path(str(out_path).format(sample_id_image))

out_path.parent.mkdir(exist_ok = True)

assert raw_img_path.exists() and segmented_img_path.exists()

################################################################################
#   Read in images and merge
################################################################################

#   Read in raw and segmented images
raw_img = tifffile.imread(raw_img_path)
segmented_img = Image.open(segmented_img_path)

#   Grab just the Abeta channel
segmented_img.seek(1)
abeta = np.array(
        segmented_img, dtype = np.uint8
    ).reshape((1, segmented_img.size[1], segmented_img.size[0])) * 255

#   Grab just the pTau channel
segmented_img.seek(2)
ptau = np.array(
        segmented_img, dtype = np.uint8
    ).reshape((1, segmented_img.size[1], segmented_img.size[0])) * 255

#   Now stack the raw images and segmented Abeta and pTau to form a single
#   8-channel image. Then write to disk
full_img = np.concatenate((raw_img, abeta, ptau), axis = 0)

with tifffile.TiffWriter(out_path) as tiff:
    tiff.write(full_img)

session_info.show()
