# Utils Module Documentation

This module provides shared utility functions for the Dormant Site Activation Pipeline.

## Module Overview

### `constants.py`
Defines all constants used throughout the pipeline:
- **Directory paths**: Base directories for data, results, figures, logs
- **Chromosome lists**: Standard human chromosomes (with and without 'chr' prefix)
- **DNA constants**: Bases, complement mapping
- **Motif tiers**: Classification levels for motif matches
- **gnomAD fields**: Standard VCF field names and descriptions
- **AlphaGenome tracks**: Epigenomic features for prediction
- **File extensions**: Standard extensions for various file types
- **Default parameters**: Sensible defaults for pipeline parameters
- **URLs**: Database URLs for JASPAR and HOCOMOCO
- **Color schemes**: Consistent colors for visualization

### `logging_utils.py`
Logging functionality for pipeline operations:

#### Functions
- **`setup_logger(name, log_file, level, format_string)`**: Configure logger with file and console output
- **`log_step(logger, step_name, start)`**: Log pipeline step boundaries with separators
- **`log_parameters(logger, params)`**: Pretty-print parameter dictionaries
- **`log_file_info(logger, file_path, description)`**: Log file existence and size
- **`create_timestamped_log(base_name, log_dir)`**: Generate timestamped log filenames

#### Classes
- **`ProgressLogger`**: Track and log progress through large iterations
  - `update(n)`: Update counter and log periodic progress with ETA

### `io.py`
File I/O operations for various formats:

#### Configuration
- **`load_config(config_path)`**: Load YAML configuration
- **`save_config(config, output_path)`**: Save configuration to YAML

#### FASTA
- **`load_fasta(fasta_path)`**: Load FASTA into dictionary (handles gzipped)
- **`write_fasta(sequences, output_path)`**: Write sequences to FASTA

#### BED
- **`load_bed(bed_path, columns)`**: Load BED file into DataFrame
- **`write_bed(df, output_path, compress)`**: Write DataFrame to BED

#### TSV
- **`load_tsv(tsv_path, **kwargs)`**: Load TSV into DataFrame
- **`write_tsv(df, output_path, compress, **kwargs)`**: Write DataFrame to TSV

#### JSON
- **`load_json(json_path)`**: Load JSON file
- **`write_json(data, output_path, indent)`**: Write dictionary to JSON

#### Utilities
- **`ensure_dir(dir_path)`**: Create directory if needed
- **`get_file_handle(file_path, mode)`**: Get handle (auto-handles gzip)
- **`reverse_complement(seq)`**: Return reverse complement of DNA sequence

## Usage Examples

### Setting up logging
```python
from utils.logging_utils import setup_logger, log_step

logger = setup_logger(
    name=__name__,
    log_file="logs/my_analysis.log",
    level="INFO"
)

log_step(logger, "Genome scanning", start=True)
# ... do work ...
log_step(logger, "Genome scanning", start=False)
```

### Progress tracking
```python
from utils.logging_utils import ProgressLogger

progress = ProgressLogger(logger, total=1000000, step=10000)

for i in range(1000000):
    # ... process item ...
    progress.update(1)
```

### Loading configuration
```python
from utils.io import load_config

config = load_config("pipeline_config.yaml")
tf_name = config['tf_name']
motif_id = config['motif_id']
```

### Working with FASTA
```python
from utils.io import load_fasta, write_fasta

# Load reference genome
genome = load_fasta("data/reference/GRCh38.fa")
chr1_seq = genome['chr1']

# Write variant sequences
variant_seqs = {
    "variant_1": "ATCGATCG...",
    "variant_2": "GCTAGCTA..."
}
write_fasta(variant_seqs, "output/variants.fa")
```

### Working with BED files
```python
from utils.io import load_bed, write_bed

# Load motif hits
motif_hits = load_bed("results/motif_scan/AP1/motif_hits.bed")

# Filter and write
strong_hits = motif_hits[motif_hits['score'] > 0.95]
write_bed(strong_hits, "results/motif_scan/AP1/strong_hits.bed")
```

### Constants usage
```python
from utils.constants import STANDARD_CHROMOSOMES, COMPLEMENT, MOTIF_TIERS

# Process only standard chromosomes
for chrom in STANDARD_CHROMOSOMES:
    process_chromosome(chrom)

# Get motif tier description
tier_name = MOTIF_TIERS[0]  # "strong"

# Reverse complement
def reverse_complement(seq):
    return ''.join(COMPLEMENT[b] for b in reversed(seq))
```

## Best Practices

1. **Always use the logging utilities** for consistent log formatting
2. **Use constants** instead of hardcoding values
3. **Use the io functions** to handle gzipped files transparently
4. **Track progress** for long-running operations with ProgressLogger
5. **Log parameters** at the start of each module for reproducibility
6. **Create timestamped logs** for each pipeline run

## Dependencies

- `pyyaml`: YAML configuration file handling
- `pandas`: Tabular data (BED, TSV)
- `biopython`: FASTA sequence handling
- Standard library: `logging`, `json`, `gzip`, `pathlib`
