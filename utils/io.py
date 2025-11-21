"""
I/O utilities for Dormant Site Activation Pipeline
"""

import gzip
import yaml
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def load_config(config_path: str) -> dict:
    """
    Load YAML configuration file.
    
    Parameters
    ----------
    config_path : str
        Path to YAML config file
        
    Returns
    -------
    dict
        Configuration dictionary
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def save_config(config: dict, output_path: str):
    """
    Save configuration to YAML file.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary
    output_path : str
        Path to output YAML file
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)


def load_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Load FASTA file into a dictionary.
    
    Parameters
    ----------
    fasta_path : str
        Path to FASTA file (can be gzipped)
        
    Returns
    -------
    dict
        Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    
    if fasta_path.endswith('.gz'):
        with gzip.open(fasta_path, 'rt') as f:
            for record in SeqIO.parse(f, 'fasta'):
                sequences[record.id] = str(record.seq)
    else:
        with open(fasta_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                sequences[record.id] = str(record.seq)
    
    return sequences


def write_fasta(sequences: Union[Dict[str, str], List[SeqRecord]], output_path: str):
    """
    Write sequences to FASTA file.
    
    Parameters
    ----------
    sequences : dict or list
        Dictionary mapping IDs to sequences, or list of SeqRecord objects
    output_path : str
        Path to output FASTA file
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    if isinstance(sequences, dict):
        records = [SeqRecord(Seq(seq), id=seq_id, description="") 
                   for seq_id, seq in sequences.items()]
    else:
        records = sequences
    
    with open(output_path, 'w') as f:
        SeqIO.write(records, f, 'fasta')


def load_bed(bed_path: str, columns: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Load BED file into a DataFrame.
    
    Parameters
    ----------
    bed_path : str
        Path to BED file (can be gzipped)
    columns : list, optional
        Column names. If None, uses standard BED columns.
        
    Returns
    -------
    pd.DataFrame
        BED file contents
    """
    if columns is None:
        columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    
    # Determine number of columns in file
    if bed_path.endswith('.gz'):
        with gzip.open(bed_path, 'rt') as f:
            first_line = f.readline().strip()
    else:
        with open(bed_path, 'r') as f:
            first_line = f.readline().strip()
    
    n_cols = len(first_line.split('\t'))
    columns = columns[:n_cols]
    
    df = pd.read_csv(bed_path, sep='\t', header=None, names=columns, comment='#')
    return df


def write_bed(df: pd.DataFrame, output_path: str, compress: bool = False):
    """
    Write DataFrame to BED file.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with BED columns
    output_path : str
        Path to output BED file
    compress : bool
        Whether to gzip compress the output
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    if compress and not output_path.endswith('.gz'):
        output_path += '.gz'
    
    df.to_csv(output_path, sep='\t', header=False, index=False)


def load_tsv(tsv_path: str, **kwargs) -> pd.DataFrame:
    """
    Load TSV file into DataFrame.
    
    Parameters
    ----------
    tsv_path : str
        Path to TSV file (can be gzipped)
    **kwargs
        Additional arguments for pd.read_csv
        
    Returns
    -------
    pd.DataFrame
        TSV contents
    """
    return pd.read_csv(tsv_path, sep='\t', **kwargs)


def write_tsv(df: pd.DataFrame, output_path: str, compress: bool = False, **kwargs):
    """
    Write DataFrame to TSV file.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to write
    output_path : str
        Path to output TSV file
    compress : bool
        Whether to gzip compress the output
    **kwargs
        Additional arguments for df.to_csv
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    if compress and not output_path.endswith('.gz'):
        output_path += '.gz'
    
    df.to_csv(output_path, sep='\t', index=False, **kwargs)


def load_json(json_path: str) -> dict:
    """
    Load JSON file.
    
    Parameters
    ----------
    json_path : str
        Path to JSON file
        
    Returns
    -------
    dict
        JSON contents
    """
    with open(json_path, 'r') as f:
        data = json.load(f)
    return data


def write_json(data: dict, output_path: str, indent: int = 2):
    """
    Write dictionary to JSON file.
    
    Parameters
    ----------
    data : dict
        Data to write
    output_path : str
        Path to output JSON file
    indent : int
        JSON indentation level
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=indent)


def ensure_dir(dir_path: str):
    """
    Create directory if it doesn't exist.
    
    Parameters
    ----------
    dir_path : str
        Path to directory
    """
    Path(dir_path).mkdir(parents=True, exist_ok=True)


def get_file_handle(file_path: str, mode: str = 'r'):
    """
    Get file handle, automatically handling gzipped files.
    
    Parameters
    ----------
    file_path : str
        Path to file
    mode : str
        File mode ('r', 'w', etc.)
        
    Returns
    -------
    file handle
        Open file handle
    """
    if file_path.endswith('.gz'):
        if 'b' not in mode:
            mode = mode.replace('r', 'rt').replace('w', 'wt')
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of DNA sequence.
    
    Parameters
    ----------
    seq : str
        DNA sequence
        
    Returns
    -------
    str
        Reverse complement
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))
