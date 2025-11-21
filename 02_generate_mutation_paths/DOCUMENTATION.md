# Module 02: Mutation Path Generation

This module generates minimal mutation paths from motif sites to consensus sequences.

## Scripts

### `consensus_from_pwm.py`
Infers the consensus sequence from a PWM by taking the maximum probability base at each position.

**Output:** `data/motifs/consensus/{TF}_consensus.txt`

### `compute_distance.py`
Calculates Hamming distance between each motif site and the consensus sequence.

**Output:** `results/mutation_paths/{TF}/distances.tsv`

### `enumerate_paths.py`
Generates all minimal mutation paths (≤3 steps) from each motif site to consensus.

**Output:** `results/mutation_paths/{TF}/paths.tsv`

## Usage

```bash
# Generate consensus
python 02_generate_mutation_paths/consensus_from_pwm.py --config pipeline_config.yaml

# Compute distances
python 02_generate_mutation_paths/compute_distance.py --config pipeline_config.yaml

# Enumerate paths
python 02_generate_mutation_paths/enumerate_paths.py --config pipeline_config.yaml
```
