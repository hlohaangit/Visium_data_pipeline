# Multi-Sample Spatial Variable Genes Analysis with Patch Extraction
import squidpy as sq
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import os
from pathlib import Path
import torch
from PIL import Image, ImageFile
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings('ignore')

# Try to import openslide for SVS support
try:
    import openslide
    from openslide import OpenSlide
    OPENSLIDE_AVAILABLE = True
    print("OpenSlide available - SVS files supported")
except ImportError:
    OPENSLIDE_AVAILABLE = False
    print("Warning: OpenSlide not available - SVS files may not work")
    print("Install with: pip install openslide-python")

# Try to import alternative libraries
try:
    import tifffile
    TIFFFILE_AVAILABLE = True
except ImportError:
    TIFFFILE_AVAILABLE = False

# Configuration
Image.MAX_IMAGE_PIXELS = None
ImageFile.LOAD_TRUNCATED_IMAGES = True

def normalize_name(name):
    """
    Normalize filename by removing hyphens and underscores for comparison
    """
    return name.replace("-", "").replace("_", "").lower()

def find_h5ad_and_wsi_files(h5ad_directory, wsi_directory):
    """
    Find matching h5ad and WSI files (ignoring hyphens and underscores in names)
    """
    h5ad_dir = Path(h5ad_directory)
    wsi_dir = Path(wsi_directory)
    
    # Find all h5ad files
    h5ad_files = list(h5ad_dir.glob("*.h5ad"))
    
    # Get all WSI subdirectories
    wsi_subdirs = [d for d in wsi_dir.iterdir() if d.is_dir()]
    
    # Create mapping of normalized names to actual directory names
    wsi_name_mapping = {}
    for wsi_subdir in wsi_subdirs:
        normalized = normalize_name(wsi_subdir.name)
        wsi_name_mapping[normalized] = wsi_subdir
    
    file_pairs = []
    for h5ad_file in h5ad_files:
        sample_name = h5ad_file.stem  # filename without extension
        normalized_sample = normalize_name(sample_name)
        
        # Look for corresponding WSI directory using normalized names
        matching_wsi_dir = None
        actual_wsi_name = None
        
        if normalized_sample in wsi_name_mapping:
            matching_wsi_dir = wsi_name_mapping[normalized_sample]
            actual_wsi_name = matching_wsi_dir.name
        
        if matching_wsi_dir and matching_wsi_dir.exists():
            # Find WSI file in the subdirectory (looking for common image formats)
            wsi_extensions = ['*.tif', '*.tiff', '*.svs', '*.ndpi', '*.vms', '*.vmu', '*.scn']
            wsi_file = None
            
            for ext in wsi_extensions:
                wsi_files = list(matching_wsi_dir.glob(ext))
                if wsi_files:
                    wsi_file = wsi_files[0]  # Take the first matching file
                    break
            
            if wsi_file:
                file_pairs.append((str(h5ad_file), str(wsi_file), sample_name))
                print(f"Found pair: {sample_name} <-> {actual_wsi_name}")
            else:
                print(f"Warning: No WSI file found in {actual_wsi_name}")
        else:
            print(f"Warning: No matching WSI directory found for {sample_name}")
            # Show available directories for debugging
            available_dirs = [d.name for d in wsi_subdirs]
            print(f"  Available WSI directories: {available_dirs}")
    
    return file_pairs

def calculate_spatial_autocorr_manual(adata, n_neighbors=6):
    """
    Manual calculation of spatial autocorrelation using Moran's I
    """
    coordinates = adata.obsm['spatial']
    
    # Build k-nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=n_neighbors+1).fit(coordinates)
    distances, indices = nbrs.kneighbors(coordinates)
    
    # Get expression matrix
    if sparse.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X.copy()
    
    n_obs, n_vars = X.shape
    moran_scores = []
    
    for gene_idx in range(n_vars):
        gene_expr = X[:, gene_idx]
        
        # Skip genes with no expression or no variance
        if np.sum(gene_expr) == 0 or np.var(gene_expr) == 0:
            moran_scores.append(0)
            continue
        
        # Calculate Moran's I
        mean_expr = np.mean(gene_expr)
        total_deviation = np.sum((gene_expr - mean_expr) ** 2)
        
        if total_deviation == 0:
            moran_scores.append(0)
            continue
        
        # Calculate spatial autocorrelation
        numerator = 0
        weight_sum = 0
        
        for i in range(n_obs):
            neighbors = indices[i, 1:]  # Skip self
            
            for neighbor_idx in neighbors:
                weight = 1  # Binary weights
                numerator += weight * (gene_expr[i] - mean_expr) * (gene_expr[neighbor_idx] - mean_expr)
                weight_sum += weight
        
        if weight_sum == 0:
            moran_scores.append(0)
        else:
            moran_i = (n_obs / weight_sum) * (numerator / total_deviation)
            moran_scores.append(moran_i)
    
    return np.array(moran_scores)

def process_single_sample(h5ad_file, sample_name):
    """
    Process a single sample and return Moran's I scores
    """
    print(f"Processing {sample_name}...")
    
    # Load data
    adata = sc.read_h5ad(h5ad_file)
    
    # Preprocessing
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Calculate spatial autocorrelation
    moran_scores = calculate_spatial_autocorr_manual(adata, n_neighbors=6)
    
    # Create results dataframe
    results_df = pd.DataFrame({
        'gene': adata.var.index,
        'moran_i': moran_scores,
        'sample': sample_name
    })
    
    return results_df, adata

def find_common_genes_and_top_svg(file_pairs, top_n=100):
    """
    Find common genes across all samples and identify top SVGs
    """
    print("Processing all samples...")
    
    all_results = []
    sample_data = {}
    all_genes_sets = []
    
    # Process each sample
    for h5ad_file, wsi_file, sample_name in file_pairs:
        try:
            results_df, adata = process_single_sample(h5ad_file, sample_name)
            all_results.append(results_df)
            sample_data[sample_name] = {
                'adata': adata,
                'wsi_file': wsi_file,
                'results': results_df
            }
            all_genes_sets.append(set(results_df['gene']))
            print(f"  - {sample_name}: {len(results_df)} genes")
        except Exception as e:
            print(f"Error processing {sample_name}: {e}")
            continue
    
    if not all_results:
        raise ValueError("No samples processed successfully")
    
    # Find common genes across all samples
    common_genes = set.intersection(*all_genes_sets)
    print(f"\nFound {len(common_genes)} common genes across {len(all_results)} samples")
    
    # Calculate average Moran's I across samples for common genes
    print("Calculating average spatial autocorrelation across samples...")
    
    gene_avg_moran = {}
    for gene in common_genes:
        moran_values = []
        for results_df in all_results:
            gene_moran = results_df[results_df['gene'] == gene]['moran_i'].iloc[0]
            moran_values.append(gene_moran)
        gene_avg_moran[gene] = np.mean(moran_values)
    
    # Get top SVGs
    sorted_genes = sorted(gene_avg_moran.items(), key=lambda x: x[1], reverse=True)
    top_svg_genes = [gene for gene, score in sorted_genes[:top_n]]
    
    print(f"\nTop {top_n} spatially variable genes:")
    for i, (gene, score) in enumerate(sorted_genes[:20]):  # Show top 20
        print(f"{i+1:2d}. {gene}: {score:.4f}")
    
    return top_svg_genes, sample_data, gene_avg_moran

def load_wsi_image(wsi_file):
    """
    Load WSI image using appropriate library based on file format
    """
    wsi_path = Path(wsi_file)
    file_extension = wsi_path.suffix.lower()
    
    print(f"  Loading WSI: {wsi_path.name} (format: {file_extension})")
    
    # Try OpenSlide first for SVS and other whole slide formats
    if OPENSLIDE_AVAILABLE and file_extension in ['.svs', '.ndpi', '.vms', '.vmu', '.scn', '.mrxs', '.tiff', '.tif']:
        try:
            print("  Using OpenSlide...")
            slide = OpenSlide(str(wsi_file))
            
            # Get the best level for processing (usually level 0 = highest resolution)
            level = 0
            dimensions = slide.level_dimensions[level]
            print(f"  Image dimensions: {dimensions[0]} x {dimensions[1]}")
            
            # For very large images, use a lower resolution level
            if dimensions[0] > 50000 or dimensions[1] > 50000:
                # Find a more manageable level
                for i, dim in enumerate(slide.level_dimensions):
                    if dim[0] <= 50000 and dim[1] <= 50000:
                        level = i
                        dimensions = dim
                        print(f"  Using level {level} with dimensions: {dimensions[0]} x {dimensions[1]}")
                        break
            
            # Read the image at the selected level
            img_pil = slide.read_region((0, 0), level, dimensions)
            
            # Convert RGBA to RGB if necessary
            if img_pil.mode == 'RGBA':
                # Create white background
                background = Image.new('RGB', img_pil.size, (255, 255, 255))
                background.paste(img_pil, mask=img_pil.split()[-1])  # Use alpha channel as mask
                img_pil = background
            elif img_pil.mode != 'RGB':
                img_pil = img_pil.convert('RGB')
            
            slide.close()
            
            # Calculate the scaling factor if we used a different level
            scale_factor = slide.level_downsamples[level] if level > 0 else 1.0
            
            return img_pil, scale_factor
            
        except Exception as e:
            print(f"  OpenSlide failed: {e}")
            # Fall through to try other methods
    
    # Try PIL for standard formats
    if file_extension in ['.tif', '.tiff', '.jpg', '.jpeg', '.png']:
        try:
            print("  Using PIL...")
            img_pil = Image.open(wsi_file)
            if img_pil.mode != 'RGB':
                img_pil = img_pil.convert('RGB')
            return img_pil, 1.0
        except Exception as e:
            print(f"  PIL failed: {e}")
    
    # Try tifffile for TIFF files
    if TIFFFILE_AVAILABLE and file_extension in ['.tif', '.tiff']:
        try:
            print("  Using tifffile...")
            img_array = tifffile.imread(wsi_file)
            
            # Handle different array shapes
            if len(img_array.shape) == 3:
                if img_array.shape[2] == 3:  # RGB
                    img_pil = Image.fromarray(img_array)
                elif img_array.shape[2] == 4:  # RGBA
                    img_pil = Image.fromarray(img_array[:, :, :3])
                else:
                    img_pil = Image.fromarray(img_array[:, :, 0])  # Take first channel
                    img_pil = img_pil.convert('RGB')
            else:  # Grayscale
                img_pil = Image.fromarray(img_array)
                img_pil = img_pil.convert('RGB')
            
            return img_pil, 1.0
        except Exception as e:
            print(f"  Tifffile failed: {e}")
    
    # If all methods fail
    raise ValueError(f"Cannot load image file: {wsi_file}. Tried OpenSlide, PIL, and tifffile.")

def extract_patches_around_spots(adata, wsi_file, patch_size=224, mpp=0.5):
    """
    Extract patches around each spot at specified resolution
    """
    try:
        # Load WSI using appropriate method
        img, scale_factor = load_wsi_image(wsi_file)
        print(f"  Successfully loaded image with scale factor: {scale_factor}")
        
    except Exception as e:
        print(f"  Failed to load WSI: {e}")
        return np.array([]), []
    
    # Get spot coordinates
    spatial_coords = adata.obsm['spatial']
    
    # Adjust coordinates for scale factor if using lower resolution level
    if scale_factor != 1.0:
        print(f"  Adjusting coordinates for scale factor: {scale_factor}")
        spatial_coords = spatial_coords / scale_factor
    
    patches = []
    valid_indices = []
    
    half_patch = patch_size // 2
    
    print(f"  Extracting patches from {len(spatial_coords)} spots...")
    
    for i, (x, y) in enumerate(spatial_coords):
        try:
            # Convert coordinates to integers
            x_int, y_int = int(x), int(y)
            
            # Calculate patch boundaries
            left = max(0, x_int - half_patch)
            top = max(0, y_int - half_patch)
            right = min(img.width, x_int + half_patch)
            bottom = min(img.height, y_int + half_patch)
            
            # Check if patch is large enough (at least 80% of target size)
            patch_width = right - left
            patch_height = bottom - top
            min_size = int(patch_size * 0.8)
            
            if patch_width >= min_size and patch_height >= min_size:
                # Extract patch
                patch = img.crop((left, top, right, bottom))
                
                # Resize to exact patch size if needed
                if patch.size != (patch_size, patch_size):
                    patch = patch.resize((patch_size, patch_size), Image.LANCZOS)
                
                # Convert to array
                patch_array = np.array(patch)
                
                # Ensure it's RGB
                if len(patch_array.shape) == 3 and patch_array.shape[2] == 3:
                    patches.append(patch_array)
                    valid_indices.append(i)
                
        except Exception as e:
            if i < 10:  # Only print first 10 errors to avoid spam
                print(f"    Error extracting patch for spot {i}: {e}")
            continue
    
    print(f"  Successfully extracted {len(patches)} patches out of {len(spatial_coords)} spots")
    
    return np.array(patches), valid_indices

def save_sample_data(sample_name, adata, top_svg_genes, patches, valid_indices, output_dir):
    """
    Save patches and expression data for a sample
    """
    sample_output_dir = Path(output_dir) / sample_name
    sample_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create patches subdirectory for PNG files
    patches_png_dir = sample_output_dir / 'patches_png'
    patches_png_dir.mkdir(exist_ok=True)
    
    # Filter adata to valid spots and top SVG genes
    adata_filtered = adata[valid_indices, :].copy()
    
    # Get indices of top SVG genes
    svg_gene_indices = [adata_filtered.var.index.get_loc(gene) for gene in top_svg_genes 
                       if gene in adata_filtered.var.index]
    
    # Extract expression data for top SVG genes
    if sparse.issparse(adata_filtered.X):
        svg_expression = adata_filtered.X[:, svg_gene_indices].toarray()
    else:
        svg_expression = adata_filtered.X[:, svg_gene_indices]
    
    # Save patches as numpy array
    np.save(sample_output_dir / 'patches.npy', patches)
    
    # Save individual patches as PNG files
    print(f"  Saving {len(patches)} individual PNG patches...")
    for i, patch in enumerate(patches):
        # Convert patch to PIL Image
        if patch.dtype != np.uint8:
            # Convert to uint8 if needed
            if patch.max() <= 1.0:
                patch_uint8 = (patch * 255).astype(np.uint8)
            else:
                patch_uint8 = patch.astype(np.uint8)
        else:
            patch_uint8 = patch
        
        # Create PIL Image
        patch_img = Image.fromarray(patch_uint8)
        
        # Save as PNG with spot index in filename
        spot_idx = valid_indices[i]
        png_filename = f"patch_{spot_idx:04d}.png"
        patch_img.save(patches_png_dir / png_filename)
    
    # Save expression data as PyTorch tensor
    svg_tensor = torch.tensor(svg_expression, dtype=torch.float32)
    torch.save({
        'expression': svg_tensor,
        'genes': [top_svg_genes[i] for i in range(len(svg_gene_indices))],
        'spot_coordinates': adata_filtered.obsm['spatial'],
        'spot_indices': valid_indices
    }, sample_output_dir / 'svg_expression.pt')
    
    # Save metadata
    metadata = {
        'sample_name': sample_name,
        'n_spots': len(valid_indices),
        'n_patches': len(patches),
        'n_svg_genes': len(svg_gene_indices),
        'patch_size': patches.shape[1:3] if len(patches) > 0 else (0, 0),
        'svg_genes': [top_svg_genes[i] for i in range(len(svg_gene_indices))]
    }
    
    pd.DataFrame([metadata]).to_csv(sample_output_dir / 'metadata.csv', index=False)
    
    print(f"  Saved {len(patches)} patches (numpy + PNG) and expression data for {len(svg_gene_indices)} SVG genes")
    
    return len(patches), len(svg_gene_indices)

def main_analysis(h5ad_directory, wsi_directory, output_directory, top_n_svg=100, patch_size=224, mpp=0.5):
    """
    Main analysis pipeline
    """
    print("="*60)
    print("MULTI-SAMPLE SPATIAL VARIABLE GENES ANALYSIS")
    print("="*60)
    
    # Create output directory
    output_dir = Path(output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find file pairs
    print("Finding h5ad and WSI file pairs...")
    file_pairs = find_h5ad_and_wsi_files(h5ad_directory, wsi_directory)
    
    if len(file_pairs) == 0:
        raise ValueError("No matching h5ad and WSI file pairs found")
    
    print(f"Found {len(file_pairs)} file pairs")
    
    # Find common genes and top SVGs
    top_svg_genes, sample_data, gene_avg_moran = find_common_genes_and_top_svg(file_pairs, top_n_svg)
    
    # Save overall results
    svg_results = pd.DataFrame({
        'gene': top_svg_genes,
        'avg_moran_i': [gene_avg_moran[gene] for gene in top_svg_genes],
        'rank': range(1, len(top_svg_genes) + 1)
    })
    svg_results.to_csv(output_dir / 'top_svg_genes.csv', index=False)
    
    # Process each sample
    print(f"\nExtracting patches and saving data...")
    summary_stats = []
    
    for sample_name, data in sample_data.items():
        print(f"Processing {sample_name}...")
        adata = data['adata']
        wsi_file = data['wsi_file']
        
        try:
            # Extract patches
            patches, valid_indices = extract_patches_around_spots(
                adata, wsi_file, patch_size, mpp
            )
            
            if len(patches) > 0:
                # Save data
                n_patches, n_genes = save_sample_data(
                    sample_name, adata, top_svg_genes, patches, valid_indices, output_dir
                )
                
                summary_stats.append({
                    'sample': sample_name,
                    'total_spots': adata.n_obs,
                    'valid_patches': n_patches,
                    'svg_genes': n_genes,
                    'success': True
                })
            else:
                print(f"  Warning: No valid patches extracted for {sample_name}")
                summary_stats.append({
                    'sample': sample_name,
                    'total_spots': adata.n_obs,
                    'valid_patches': 0,
                    'svg_genes': 0,
                    'success': False
                })
                
        except Exception as e:
            print(f"  Error processing {sample_name}: {e}")
            summary_stats.append({
                'sample': sample_name,
                'total_spots': adata.n_obs if 'adata' in locals() else 0,
                'valid_patches': 0,
                'svg_genes': 0,
                'success': False
            })
    
    # Save summary
    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(output_dir / 'processing_summary.csv', index=False)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"Output directory: {output_dir}")
    print(f"Top {top_n_svg} SVG genes saved to: top_svg_genes.csv")
    print(f"Processing summary saved to: processing_summary.csv")
    print(f"Individual sample data saved in subdirectories")
    
    # Print summary statistics
    total_samples = len(summary_df)
    successful_samples = summary_df['success'].sum()
    total_patches = summary_df['valid_patches'].sum()
    
    print(f"\nSummary Statistics:")
    print(f"  Total samples: {total_samples}")
    print(f"  Successful samples: {successful_samples}")
    print(f"  Total patches extracted: {total_patches}")
    print(f"  Average patches per successful sample: {total_patches/max(successful_samples, 1):.1f}")
    
    return top_svg_genes, summary_df

# Example usage
if __name__ == "__main__":
    # Set your directories here
    h5ad_directory = '/blue/pinaki.sarder/j.fermin/CellAtlas/data/Kidney/H5AD_files/'
    wsi_directory = '/blue/pinaki.sarder/j.fermin/Annotations/Data/'

    output_directory = "/orange/pinaki.sarder/h.lohaan/Hari_data_pipeline/"
    
    # Run analysis
    top_svg_genes, summary = main_analysis(
        h5ad_directory=h5ad_directory,
        wsi_directory=wsi_directory,
        output_directory=output_directory,
        top_n_svg=100,
        patch_size=224,
        mpp=0.5
    )
    
    print(f"\nTop 10 SVG genes: {top_svg_genes[:10]}")

# Additional utility functions for loading saved data
def load_sample_data(sample_path):
    """
    Load saved sample data
    """
    sample_path = Path(sample_path)
    
    # Load patches
    patches = np.load(sample_path / 'patches.npy')
    
    # Load expression data
    expression_data = torch.load(sample_path / 'svg_expression.pt')
    
    # Load metadata
    metadata = pd.read_csv(sample_path / 'metadata.csv').iloc[0].to_dict()
    
    return {
        'patches': patches,
        'expression': expression_data['expression'],
        'genes': expression_data['genes'],
        'coordinates': expression_data['spot_coordinates'],
        'spot_indices': expression_data['spot_indices'],
        'metadata': metadata
    }

def visualize_sample_patches(sample_path, n_patches=9):
    """
    Visualize patches from a sample
    """
    data = load_sample_data(sample_path)
    patches = data['patches']
    
    if len(patches) == 0:
        print("No patches to visualize")
        return
    
    n_patches = min(n_patches, len(patches))
    n_cols = 3
    n_rows = (n_patches + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    for i in range(n_patches):
        axes[i].imshow(patches[i])
        axes[i].set_title(f'Patch {i+1}')
        axes[i].axis('off')
    
    # Hide unused subplots
    for i in range(n_patches, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    plt.show()
    
    print(f"Sample: {data['metadata']['sample_name']}")
    print(f"Total patches: {len(patches)}")
    print(f"Patch size: {patches[0].shape}")
    print(f"SVG genes: {len(data['genes'])}")