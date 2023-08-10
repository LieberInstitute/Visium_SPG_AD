import scanpy as sc
import numpy as np
from pyhere import here
from pathlib import Path
from PIL import Image
import json
import session_info

################################################################################
#   Variable definitions
################################################################################

processed_dir = Path(here('processed-data', '21_spot_deconvo'))
plot_dir = Path(here('plots', '21_spot_deconvo'))

plot_dir.mkdir(parents = True, exist_ok = True)
processed_dir.mkdir(parents = True, exist_ok = True)

#   Directory containing hires image and a JSON containing scale factors and
#   spot size for a given sample. Here '{}' will be replaced by a single
#   sample name
spaceranger_dir = here(
    'processed-data', 'spaceranger', '{}', 'outs', 'spatial'
)

#   Naming conventions used for different columns in the AnnDatas
sample_id_var = 'sample_id'          # in spatial object only
ensembl_id_var = 'gene_id'           # in both spatial and single-cell objects
gene_symbol_var = 'gene_name'        # in both spatial and single-cell objects

cell_type_var = "broad.cell.type"

spatial_coords_names = ['pxl_col_in_fullres', 'pxl_row_in_fullres']

################################################################################
#   Preprocessing
################################################################################

#  Load AnnDatas
print('Loading AnnDatas...')
adata_vis = sc.read_h5ad(processed_dir / 'adata_spatial.h5ad')
adata_ref = sc.read_h5ad(processed_dir / 'adata_mathys.h5ad')

adata_vis.obs['sample'] = adata_vis.obs[sample_id_var]

# rename genes to ENSEMBL
adata_vis.var['SYMBOL'] = adata_vis.var[gene_symbol_var]
adata_vis.var_names = adata_vis.var[ensembl_id_var]
adata_vis.var_names.name = None

#   Spatial AnnData needs unique indices. Rather than using barcode (repeated
#   for every sample), use "key" (barcode + sample ID)
adata_vis.obs_names = adata_vis.obs['key']
adata_vis.obs_names.name = None

# Use ENSEMBL as gene IDs to make sure IDs are unique and correctly matched
adata_ref.var['SYMBOL'] = adata_ref.var[gene_symbol_var]
adata_ref.var.index = adata_ref.var[ensembl_id_var]
adata_ref.var_names = adata_ref.var[ensembl_id_var]
adata_ref.var.index.name = None

#   Subset to marker genes
with open(processed_dir / 'markers.txt', 'r') as f:
    selected = f.read().splitlines()

adata_ref = adata_ref[:, selected].copy()

#-------------------------------------------------------------------------------
#   Attach hi-res images and scaleFactors to spatial AnnData
#-------------------------------------------------------------------------------

adata_vis.uns['spatial'] = {}

for sample_id in adata_vis.obs['sample'].cat.categories:
    #   Path to JSON from spaceranger including spot size for this sample
    json_path = here(
        str(spaceranger_dir).format(sample_id), 'scalefactors_json.json'
    )

    with open(json_path) as f: 
        json_data = json.load(f)

    #   Read in high-res image as numpy array with values in [0, 1] rather than
    #   [0, 255], then attach to AnnData object
    img_path = str(
        here(
            str(spaceranger_dir).format(sample_id), 'tissue_hires_image.png'
        )
    )
    img_arr = np.array(Image.open(img_path), dtype = np.float32) / 256

    #   Store image and scalefactors in AnnData as squidpy expects
    adata_vis.uns['spatial'][sample_id] = {
        'scalefactors': json_data,
        'images' : { 'hires' : img_arr }
    }

#-------------------------------------------------------------------------------
#   Attach spatialCoords to spatial AnnData
#-------------------------------------------------------------------------------

#   Correct how spatialCoords are stored. Currently, they are a pandas
#   DataFrame, with the columns potentially in the wrong order (depending on the
#   version of SpatialExperiment used in R). We need them as a numpy array.
adata_vis.obsm['spatial'] = np.array(
    adata_vis.obsm['spatial'][spatial_coords_names]
)

#-------------------------------------------------------------------------------
#   Save AnnDatas
#-------------------------------------------------------------------------------

adata_vis.write_h5ad(processed_dir / 'adata_spatial_orig.h5ad')
adata_ref.write_h5ad(processed_dir / 'adata_mathys_orig.h5ad')

session_info.show(html=False)
