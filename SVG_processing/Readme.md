# Multi-Sample Spatial Variable Genes Analysis with Patch Extraction

A comprehensive pipeline for analyzing spatially variable genes (SVGs) across multiple spatial transcriptomics samples and extracting corresponding histological patches from whole slide images (WSIs).

## Overview

This pipeline processes multiple spatial transcriptomics samples to:
1. Identify spatially variable genes using Moran's I statistic
2. Find common genes across all samples
3. Extract histological patches around each spot
4. Save organized data for downstream machine learning applications

## Mathematical Foundation

### Moran's I Statistic

The pipeline uses Moran's I to measure spatial autocorrelation for each gene:

```
I = (n/W) * (Σᵢ Σⱼ wᵢⱼ(xᵢ - x̄)(xⱼ - x̄)) / (Σᵢ(xᵢ - x̄)²)
```

Where:
- `n` = number of spatial locations (spots)
- `W` = sum of all spatial weights
- `wᵢⱼ` = spatial weight between locations i and j (binary: 1 if neighbors, 0 otherwise)
- `xᵢ` = gene expression at location i
- `x̄` = mean gene expression across all locations

**Interpretation:**
- `I > 0`: Positive spatial autocorrelation (similar values cluster together)
- `I = 0`: No spatial autocorrelation (random distribution)
- `I < 0`: Negative spatial autocorrelation (dissimilar values cluster together)

### Spatial Neighbors

The pipeline uses k-nearest neighbors (default k=6) to define spatial relationships:
- For each spot, find the 6 nearest spots based on Euclidean distance
- Create binary adjacency matrix where `wᵢⱼ = 1` if spots are neighbors

## Installation and Dependencies

### Required Libraries

```bash
# Core dependencies
pip install squidpy scanpy pandas numpy scipy matplotlib scikit-learn torch pillow

# For WSI support (recommended)
pip install openslide-python tifffile

# Alternative for some systems
conda install -c conda-forge openslide-python
```

### Optional Dependencies

```bash
# For additional image format support
pip install tifffile

# For SVS file support on Windows
# Download and install OpenSlide binaries from: https://openslide.org/download/
```

## File Structure Requirements

### Input Structure
```
your_project/
├── h5ad_files/
│   ├── sample1.h5ad
│   ├── sample2.h5ad
│   └── sample3.h5ad
└── wsi_files/
    ├── sample1/
    │   └── sample1.svs (or .tif, .tiff, .ndpi, etc.)
    ├── sample2/
    │   └── sample2.svs
    └── sample3/
        └── sample3.svs
```

### Output Structure
```
output_directory/
├── top_svg_genes.csv
├── processing_summary.csv
├── sample1/
│   ├── patches.npy
│   ├── patches_png/
│   │   ├── patch_0000.png
│   │   ├── patch_0001.png
│   │   └── ...
│   ├── svg_expression.pt
│   └── metadata.csv
├── sample2/
│   └── ...
└── sample3/
    └── ...
```

## Function Documentation

### Core Analysis Functions

#### `calculate_spatial_autocorr_manual(adata, n_neighbors=6)`
**Purpose**: Calculates Moran's I spatial autocorrelation for all genes in a sample.

**Parameters**:
- `adata`: AnnData object with spatial coordinates in `obsm['spatial']`
- `n_neighbors`: Number of nearest neighbors to consider (default: 6)

**Returns**: Array of Moran's I scores for each gene

**Mathematical Process**:
1. Build k-nearest neighbor graph from spatial coordinates
2. For each gene, calculate mean expression and total deviation
3. Compute spatial autocorrelation using Moran's I formula
4. Return scores for all genes

#### `process_single_sample(h5ad_file, sample_name)`
**Purpose**: Processes a single spatial transcriptomics sample.

**Parameters**:
- `h5ad_file`: Path to H5AD file
- `sample_name`: Identifier for the sample

**Returns**: 
- `results_df`: DataFrame with genes, Moran's I scores, and sample name
- `adata`: Processed AnnData object

**Processing Steps**:
1. Load H5AD file
2. Calculate QC metrics
3. Normalize to 10,000 counts per spot
4. Log-transform expression values
5. Calculate spatial autocorrelation
6. Return results

#### `find_common_genes_and_top_svg(file_pairs, top_n=100)`
**Purpose**: Identifies common genes across samples and ranks them by average spatial autocorrelation.

**Parameters**:
- `file_pairs`: List of (h5ad_file, wsi_file, sample_name) tuples
- `top_n`: Number of top SVGs to return (default: 100)

**Returns**:
- `top_svg_genes`: List of top spatially variable genes
- `sample_data`: Dictionary with processed data for each sample
- `gene_avg_moran`: Dictionary with average Moran's I scores

**Algorithm**:
1. Process all samples individually
2. Find intersection of genes across all samples
3. Calculate average Moran's I for each common gene
4. Rank genes by average spatial autocorrelation
5. Return top N genes

### Image Processing Functions

#### `load_wsi_image(wsi_file)`
**Purpose**: Loads whole slide images using appropriate libraries based on file format.

**Parameters**:
- `wsi_file`: Path to WSI file

**Returns**:
- `img_pil`: PIL Image object
- `scale_factor`: Scaling factor if lower resolution level was used

**Supported Formats**:
- `.svs`: Aperio ScanScope Virtual Slide
- `.ndpi`: Hamamatsu NanoZoomer
- `.tif/.tiff`: TIFF images
- `.vms/.vmu`: Hamamatsu Virtual Microscopy
- `.scn`: Leica SCN
- `.mrxs`: 3DHistech MIRAX

**Loading Strategy**:
1. Try OpenSlide for whole slide formats
2. Use lower resolution level for very large images (>50,000 pixels)
3. Fall back to PIL for standard formats
4. Use tifffile as last resort for TIFF files

#### `extract_patches_around_spots(adata, wsi_file, patch_size=224, mpp=0.5)`
**Purpose**: Extracts square patches around each spatial transcriptomics spot.

**Parameters**:
- `adata`: AnnData object with spatial coordinates
- `wsi_file`: Path to corresponding WSI file
- `patch_size`: Size of extracted patches in pixels (default: 224)
- `mpp`: Microns per pixel (currently unused, for future resolution matching)

**Returns**:
- `patches`: NumPy array of extracted patches
- `valid_indices`: Indices of spots with successfully extracted patches

**Extraction Process**:
1. Load WSI image
2. Adjust coordinates for scale factor if needed
3. For each spot coordinate:
   - Calculate patch boundaries (centered on spot)
   - Extract patch ensuring minimum 80% of target size
   - Resize to exact patch size if needed
   - Convert to RGB format
4. Return valid patches and corresponding indices

### Data Management Functions

#### `save_sample_data(sample_name, adata, top_svg_genes, patches, valid_indices, output_dir)`
**Purpose**: Saves processed data for a single sample in organized format.

**Parameters**:
- `sample_name`: Identifier for the sample
- `adata`: AnnData object
- `top_svg_genes`: List of top spatially variable genes
- `patches`: NumPy array of extracted patches
- `valid_indices`: Indices of valid spots
- `output_dir`: Base output directory

**Saves**:
- `patches.npy`: NumPy array of all patches
- `patches_png/`: Individual PNG files for each patch
- `svg_expression.pt`: PyTorch tensor with expression data
- `metadata.csv`: Sample metadata

#### `load_sample_data(sample_path)`
**Purpose**: Loads previously saved sample data.

**Parameters**:
- `sample_path`: Path to sample directory

**Returns**: Dictionary with patches, expression data, genes, coordinates, and metadata

### Utility Functions

#### `normalize_name(name)`
**Purpose**: Normalizes filenames by removing hyphens and underscores for flexible matching.

#### `find_h5ad_and_wsi_files(h5ad_directory, wsi_directory)`
**Purpose**: Finds matching pairs of H5AD and WSI files using normalized name matching.

#### `visualize_sample_patches(sample_path, n_patches=9)`
**Purpose**: Visualizes a grid of patches from a saved sample for quality control.

## Usage Instructions

### Basic Usage

```python
from your_script import main_analysis

# Set your directories
h5ad_directory = '/path/to/h5ad_files/'
wsi_directory = '/path/to/wsi_files/'
output_directory = '/path/to/output/'

# Run analysis
top_svg_genes, summary = main_analysis(
    h5ad_directory=h5ad_directory,
    wsi_directory=wsi_directory,
    output_directory=output_directory,
    top_n_svg=100,        # Number of top SVGs to extract
    patch_size=224,       # Size of patches in pixels
    mpp=0.5              # Microns per pixel (future use)
)

print(f"Top 10 SVG genes: {top_svg_genes[:10]}")
```

### Advanced Usage

#### Loading Saved Data

```python
from your_script import load_sample_data, visualize_sample_patches

# Load a specific sample
sample_data = load_sample_data('/path/to/output/sample1/')

# Access components
patches = sample_data['patches']           # NumPy array (n_patches, 224, 224, 3)
expression = sample_data['expression']     # PyTorch tensor (n_patches, n_genes)
genes = sample_data['genes']               # List of gene names
coordinates = sample_data['coordinates']   # Spot coordinates
metadata = sample_data['metadata']         # Sample metadata

# Visualize patches
visualize_sample_patches('/path/to/output/sample1/', n_patches=9)
```

#### Batch Processing

```python
import pandas as pd
from pathlib import Path

# Load processing summary
summary = pd.read_csv('/path/to/output/processing_summary.csv')

# Process only successful samples
successful_samples = summary[summary['success'] == True]

for _, row in successful_samples.iterrows():
    sample_name = row['sample']
    sample_path = Path('/path/to/output') / sample_name
    
    # Load and process each sample
    data = load_sample_data(sample_path)
    print(f"Sample {sample_name}: {len(data['patches'])} patches")
```

### Parameter Tuning

#### Spatial Neighbors
```python
# Adjust number of neighbors for spatial autocorrelation
# More neighbors = smoother spatial relationships
# Fewer neighbors = more local spatial patterns
moran_scores = calculate_spatial_autocorr_manual(adata, n_neighbors=8)
```

#### Patch Size
```python
# Adjust patch size based on your analysis needs
# Larger patches = more tissue context
# Smaller patches = more focused on spot location
patches, indices = extract_patches_around_spots(adata, wsi_file, patch_size=512)
```

## Output Files Description

### `top_svg_genes.csv`
Contains ranked list of spatially variable genes:
- `gene`: Gene name
- `avg_moran_i`: Average Moran's I score across samples
- `rank`: Rank by spatial autocorrelation

### `processing_summary.csv`
Processing statistics for each sample:
- `sample`: Sample name
- `total_spots`: Total number of spots in sample
- `valid_patches`: Number of successfully extracted patches
- `svg_genes`: Number of SVG genes with expression data
- `success`: Whether processing was successful

### Per-Sample Files

#### `patches.npy`
NumPy array with shape `(n_patches, patch_size, patch_size, 3)` containing RGB patch data.

#### `patches_png/`
Individual PNG files for each patch, named `patch_XXXX.png` where XXXX is the spot index.

#### `svg_expression.pt`
PyTorch tensor file containing:
- `expression`: Tensor with shape `(n_patches, n_svg_genes)`
- `genes`: List of gene names
- `spot_coordinates`: Original spatial coordinates
- `spot_indices`: Indices of spots with valid patches

#### `metadata.csv`
Sample metadata including counts and gene lists.

## Troubleshooting

### Common Issues

#### OpenSlide Installation
```bash
# Ubuntu/Debian
sudo apt-get install openslide-tools python3-openslide

# macOS
brew install openslide

# Windows
# Download from https://openslide.org/download/
```

#### Memory Issues
- Reduce patch size for large datasets
- Process samples sequentially rather than loading all at once
- Use lower resolution WSI levels for very large images

#### File Format Issues
- Ensure WSI files are in supported formats
- Check file permissions and paths
- Verify H5AD files contain 'spatial' coordinates in `obsm`

