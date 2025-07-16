# Visium Data Processing Pipeline

This repository contains two main scripts for processing Visium spatial transcriptomics data with whole slide images (WSI), specifically designed for kidney tissue analysis.

## Overview

The pipeline processes Visium spatial transcriptomics data by:
1. Extracting tissue patches from whole slide images at spot locations
2. Associating cell type information with spatial coordinates
3. Creating a structured dataset for downstream analysis

## Scripts

### 1. `visium_processor.py` (Main Processing Script)

**Purpose**: Optimized batch processing of Visium data to extract tissue patches and associate them with cell type information.

**Key Features**:
- **Multi-format support**: Handles both SVS and TIFF whole slide images
- **Batch processing**: Processes multiple spots simultaneously for efficiency
- **Parallel processing**: Uses multiprocessing for faster execution
- **Memory optimization**: Efficient memory management for large datasets
- **Flexible patch extraction**: Configurable patch size and target microns per pixel (MPP)

**Main Components**:

#### `OptimizedWSIProcessor` Class
- Initializes and manages whole slide image reading
- Calculates appropriate scale factors for target resolution
- Extracts patches from specified coordinates
- Handles both TIFF and SVS file formats

#### `process_visium_data_optimized` Function
- Main orchestration function for the entire pipeline
- Manages batch processing and parallel execution
- Saves extracted patches and cell type vectors
- Generates processing summary statistics

**Input Files**:
- `h5ad_file`: Scanpy AnnData file containing spatial coordinates and gene expression data
- `wsi_file`: Whole slide image file (.svs or .tiff)
- `csv_file`: Cell type counts per spot
- `output_dir`: Directory for saving results

**Output Structure**:
```
output_dir/
├── patches/           # PNG images of tissue patches
├── arrays/           # NumPy arrays of patches (optional)
├── cell_sub_types/   # PyTorch tensors of cell type vectors
└── summary.txt       # Processing statistics
```

**Key Parameters**:
- `patch_size`: Size of extracted patches (default: 224x224)
- `target_mpp`: Target microns per pixel (default: 0.5)
- `batch_size`: Number of spots processed per batch (default: 100)
- `n_workers`: Number of parallel workers

### 2. `cell_type_dictionary.py` (Cell Type Mapping Script)

**Purpose**: Creates a mapping dictionary from main cell types to their subtypes.

**Functionality**:
- Reads a CSV file containing cell type hierarchies
- Parses main types and their corresponding subtypes
- Creates a dictionary mapping for easy lookup

**Input Format** (CSV):
```
Main_Types,Sub_Types
Epithelial,Proximal_Tubule.Distal_Tubule.Collecting_Duct
Immune,T_Cell.B_Cell.Macrophage
Stromal,Fibroblast.Endothelial
```

**Output**: Python dictionary with main types as keys and subtype lists as values


## Usage

### Basic Usage

```python
from visium_processor import process_visium_data_optimized

# Define file paths
h5ad_file = "path/to/spatial_data.h5ad"
wsi_file = "path/to/whole_slide_image.svs"
csv_file = "path/to/cell_types.csv"
output_dir = "path/to/output"

# Process the data
summary = process_visium_data_optimized(
    h5ad_file=h5ad_file,
    wsi_file=wsi_file,
    csv_file=csv_file,
    output_dir=output_dir,
    patch_size=224,
    target_mpp=0.5,
    batch_size=100,
    n_workers=8
)
```

### Cell Type Dictionary Creation

```python
from cell_type_dictionary import create_cell_type_dictionary

# Create cell type mapping
cell_type_dict = create_cell_type_dictionary("path/to/cell_types.csv")
print(cell_type_dict)
```

## Data Flow

1. **Load Spatial Data**: Read H5AD file containing spot coordinates and gene expression
2. **Initialize WSI Processor**: Set up whole slide image reading with appropriate scaling
3. **Batch Processing**: Process spots in batches for memory efficiency
4. **Patch Extraction**: Extract tissue patches at each spot location
5. **Cell Type Association**: Link cell type information to each spot
6. **Save Results**: Store patches, arrays, and cell type vectors

## Performance Optimization

### Memory Management
- Uses batch processing to limit memory usage
- Implements proper cleanup of slide resources
- Configurable batch sizes for different system capabilities

### Parallel Processing
- Multi-worker processing for CPU-bound tasks
- Optimized for systems with multiple cores
- Adjustable worker count based on available resources

### File I/O Optimization
- Efficient patch extraction algorithms
- Option for HDF5 consolidation for faster access
- Compressed storage formats where appropriate

## File Structure Requirements

Your data should be organized as follows:
```
project/
├── spatial_data/
│   ├── sample1.h5ad
│   ├── sample2.h5ad
│   └── ...
├── wsi_images/
│   ├── sample1.svs
│   ├── sample2.svs
│   └── ...
├── cell_type_data/
│   ├── sample1_counts.csv
│   ├── sample2_counts.csv
│   └── Cell_SubTypes_Grouped.csv
└── output/
    ├── sample1/
    ├── sample2/
    └── ...
```

## Troubleshooting

### Common Issues

1. **Memory Errors**: Reduce batch_size or n_workers
2. **File Not Found**: Check file paths and permissions
3. **OpenSlide Errors**: Ensure OpenSlide is properly installed
4. **Slow Processing**: Increase batch_size and n_workers based on system capabilities

### Performance Tuning

- **For large datasets**: Increase batch_size (200-500)
- **For memory-limited systems**: Decrease batch_size (50-100)
- **For faster processing**: Increase n_workers (up to CPU cores - 1)

## Output Validation

The pipeline generates a summary file with key statistics:
- Total spots processed
- Successful extractions
- Processing time
- Spots per second
