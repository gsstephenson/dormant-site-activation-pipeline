#!/usr/bin/env python3
"""
Dormant Site Activation Pipeline - Automated Runner

This script automates execution of pipeline modules in the correct order.
Supports selective module execution via the -m flag.

Usage:
    # Run all modules (0-5)
    python run_pipeline.py --config pipeline_config.yaml

    # Run specific modules
    python run_pipeline.py --config pipeline_config.yaml -m 2,3,4,5

    # Smoke test (validates setup without running full pipeline)
    python run_pipeline.py --config pipeline_config.yaml --smoke-test -m 2,3,4,5

    # Dry run (shows what would be executed)
    python run_pipeline.py --config pipeline_config.yaml -m 2,3,4,5 --dry-run

    # Enable logging (saves to logs/pipeline_run_TIMESTAMP.log)
    python run_pipeline.py --config pipeline_config.yaml -m 2,3,4,5 --log

    # Custom log file
    python run_pipeline.py --config pipeline_config.yaml -m 2,3,4,5 --log-file my_run.log

Author: George Stephenson
Date: November 2025
"""

import argparse
import subprocess
import sys
import os
import time
import io
from pathlib import Path
from datetime import datetime, timedelta
import yaml


# =============================================================================
# LOGGING UTILITY - Tee output to both terminal and file
# =============================================================================

class TeeLogger:
    """
    A class that duplicates stdout/stderr to both terminal and a log file.
    
    This allows real-time viewing of output while also capturing everything
    to a file for later debugging.
    """
    
    def __init__(self, log_file: Path):
        self.log_file = log_file
        self.terminal = sys.stdout
        self.log = open(log_file, 'w', buffering=1)  # Line buffered
        
        # Write header
        self.log.write(f"{'='*80}\n")
        self.log.write(f"PIPELINE LOG - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        self.log.write(f"Log file: {log_file}\n")
        self.log.write(f"{'='*80}\n\n")
        self.log.flush()
    
    def write(self, message):
        self.terminal.write(message)
        self.terminal.flush()
        self.log.write(message)
        self.log.flush()
    
    def flush(self):
        self.terminal.flush()
        self.log.flush()
    
    def close(self):
        self.log.write(f"\n{'='*80}\n")
        self.log.write(f"LOG ENDED - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        self.log.write(f"{'='*80}\n")
        self.log.close()


def setup_logging(log_dir: Path, modules: list) -> Path:
    """
    Setup logging to a timestamped file.
    
    Parameters
    ----------
    log_dir : Path
        Directory to store log files
    modules : list
        List of modules being run (for filename)
        
    Returns
    -------
    Path
        Path to the log file
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # Create timestamped log filename
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    modules_str = '-'.join(map(str, modules))
    log_file = log_dir / f"pipeline_run_{timestamp}_modules_{modules_str}.log"
    
    return log_file


# =============================================================================
# MODULE DEFINITIONS
# =============================================================================

MODULES = {
    0: {
        'name': 'Data Fetching',
        'description': 'Download reference genome, motifs, verify gnomAD',
        'scripts': [
            # These are typically run once manually, but included for completeness
            # ('bash', '00_fetch_data/download_reference.sh'),
            # ('python', '00_fetch_data/download_motifs.py'),
        ],
        'estimated_time': '~5 minutes (if data already present)',
        'skip_by_default': True,  # Usually already done
    },
    1: {
        'name': 'Motif Scanning',
        'description': 'Scan genome for AP1 motif matches and tier by strength',
        'scripts': [
            ('python', '01_scan_motifs/scan_genome_fimo.py', ['--config', '{config}']),
            ('python', '01_scan_motifs/tier_sites.py', ['--config', '{config}']),
        ],
        'estimated_time': '~2-4 hours',
        'skip_by_default': True,  # Usually already done
    },
    2: {
        'name': 'Mutation Path Enumeration',
        'description': 'Generate mutation paths from dormant sites to consensus (with strand-aware alleles)',
        'scripts': [
            ('python', '02_generate_mutation_paths/consensus_from_pwm.py', ['--config', '{config}']),
            ('python', '02_generate_mutation_paths/compute_distance.py', ['--config', '{config}']),
            ('python', '02_generate_mutation_paths/enumerate_paths.py', ['--config', '{config}']),
        ],
        'estimated_time': '~30-60 minutes',
        'skip_by_default': False,
    },
    3: {
        'name': 'gnomAD Intersection',
        'description': 'Query gnomAD for population variants matching mutation paths',
        'scripts': [
            ('python', '03_intersect_gnomad/query_gnomad_vcfs.py', ['--config', '{config}']),
        ],
        'estimated_time': '~3-4 hours',
        'skip_by_default': False,
    },
    4: {
        'name': 'AlphaGenome Scoring',
        'description': 'Score functional impact using AlphaGenome API',
        'scripts': [
            ('python', '04_run_alphagenome/prepare_unique_variants.py', []),
            ('python', '04_run_alphagenome/run_alphagenome_scoring.py', ['--config', '{config}']),
        ],
        'estimated_time': '~5-6 hours',
        'skip_by_default': False,
        'requires_env': 'ALPHA_GENOME_KEY',
    },
    5: {
        'name': 'Activation Landscape',
        'description': 'Compute and visualize the 2D activation landscape',
        'scripts': [
            ('python', '05_compute_activation_landscape/compute_activation_landscape.py', [
                '--predictions', 'results/alphagenome/AP1/predictions.parquet',
                '--gnomad', 'results/gnomad_intersection/AP1/paths_with_gnomad.tsv',
                '--output', 'results/landscape/AP1'
            ]),
            ('python', '05_compute_activation_landscape/plot_activation_landscape.py', [
                '--input', 'results/landscape/AP1/AP1_activation_landscape.tsv',
                '--output-dir', 'figures/landscape'
            ]),
        ],
        'estimated_time': '~1-2 minutes',
        'skip_by_default': False,
    },
    6: {
        'name': 'Disease Overlap Analysis',
        'description': 'ClinVar intersection and GWAS proximity enrichment - validates novel disease candidates',
        'scripts': [
            ('python', '06_disease_overlap/disease_overlap.py', [
                '--ap1-landscape', 'results/landscape/{tf_name}/{tf_name}_activation_landscape.tsv',
                '--clinvar-vcf', 'data/clinvar/raw/clinvar.vcf.gz',
                '--gwas-catalog', 'data/gwas/raw/gwas_catalog_associations.tsv',
                '--output-dir', 'results/disease_overlap/{tf_name}',
                '--gwas-window', '1000'
            ]),
        ],
        'estimated_time': '~1-2 minutes',
        'skip_by_default': False,
    },
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_config(config_path: Path) -> dict:
    """Load pipeline configuration."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def format_duration(seconds: float) -> str:
    """Format duration in human-readable format."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f}m"
    else:
        return f"{seconds/3600:.1f}h"


def run_command(cmd: list, cwd: Path, dry_run: bool = False, env: dict = None, log_file: Path = None) -> bool:
    """
    Run a command and return success status.
    
    Parameters
    ----------
    cmd : list
        Command and arguments
    cwd : Path
        Working directory
    dry_run : bool
        If True, print command but don't execute
    env : dict
        Environment variables to add
    log_file : Path
        If provided, tee output to this log file
        
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    cmd_str = ' '.join(cmd)
    
    if dry_run:
        print(f"  [DRY RUN] Would execute: {cmd_str}")
        return True
    
    print(f"  Executing: {cmd_str}")
    sys.stdout.flush()
    
    # Prepare environment
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    
    try:
        if log_file:
            # Use subprocess with real-time output capture via Popen
            # This streams to both terminal and log file
            process = subprocess.Popen(
                cmd,
                cwd=cwd,
                env=run_env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1  # Line buffered
            )
            
            # Stream output line by line
            with open(log_file, 'a') as log:
                for line in process.stdout:
                    print(line, end='')  # Print to terminal
                    sys.stdout.flush()
                    log.write(line)  # Write to log
                    log.flush()
            
            process.wait()
            return process.returncode == 0
        else:
            # Standard execution without logging
            result = subprocess.run(
                cmd,
                cwd=cwd,
                env=run_env,
                capture_output=False,  # Let output stream to terminal
                text=True
            )
            return result.returncode == 0
    except Exception as e:
        print(f"  ERROR: {e}")
        return False


def check_prerequisites(module_num: int, base_dir: Path, config: dict) -> tuple:
    """
    Check if prerequisites for a module are met.
    
    Returns
    -------
    tuple
        (success: bool, message: str)
    """
    tf_name = config.get('tf_name', 'AP1')
    
    prerequisites = {
        2: [
            ('results/motif_scan/{tf}/motif_hits_tiered.tsv', 'Module 01 output'),
        ],
        3: [
            ('results/mutation_paths/{tf}/paths.tsv', 'Module 02 output'),
        ],
        4: [
            ('results/gnomad_intersection/{tf}/paths_with_gnomad.tsv', 'Module 03 output'),
        ],
        5: [
            ('results/alphagenome/{tf}/predictions.parquet', 'Module 04 output'),
            ('results/gnomad_intersection/{tf}/paths_with_gnomad.tsv', 'Module 03 output'),
        ],
        6: [
            ('results/landscape/{tf}/{tf}_activation_landscape.tsv', 'Module 05 output'),
            ('data/clinvar/raw/clinvar.vcf.gz', 'ClinVar VCF'),
            ('data/gwas/raw/gwas_catalog_associations.tsv', 'GWAS Catalog'),
        ],
    }
    
    if module_num not in prerequisites:
        return True, "No prerequisites"
    
    missing = []
    for path_template, description in prerequisites[module_num]:
        path = base_dir / path_template.format(tf=tf_name)
        if not path.exists():
            missing.append(f"{description}: {path}")
    
    if missing:
        return False, "Missing prerequisites:\n    " + "\n    ".join(missing)
    
    return True, "All prerequisites met"


def check_environment(module_num: int) -> tuple:
    """
    Check if required environment variables are set.
    
    Returns
    -------
    tuple
        (success: bool, message: str)
    """
    module = MODULES.get(module_num, {})
    required_env = module.get('requires_env')
    
    if not required_env:
        return True, "No environment variables required"
    
    if os.environ.get(required_env):
        return True, f"{required_env} is set"
    else:
        return False, f"Missing environment variable: {required_env}"


# =============================================================================
# SMOKE TEST
# =============================================================================

def run_smoke_test(base_dir: Path, config_path: Path, modules: list) -> bool:
    """
    Run smoke tests to validate pipeline setup.
    
    Returns
    -------
    bool
        True if all tests pass
    """
    print("\n" + "=" * 80)
    print("SMOKE TEST - Validating Pipeline Setup")
    print("=" * 80)
    
    config = load_config(config_path)
    tf_name = config.get('tf_name', 'AP1')
    all_passed = True
    
    # Test 1: Config file
    print("\n[1/7] Checking configuration file...")
    if config_path.exists():
        print(f"  ✓ Config found: {config_path}")
        print(f"  ✓ TF name: {tf_name}")
    else:
        print(f"  ✗ Config not found: {config_path}")
        all_passed = False
    
    # Test 2: Python environment
    print("\n[2/7] Checking Python environment...")
    try:
        import pandas
        import numpy
        import yaml
        print(f"  ✓ pandas: {pandas.__version__}")
        print(f"  ✓ numpy: {numpy.__version__}")
        print(f"  ✓ yaml: available")
    except ImportError as e:
        print(f"  ✗ Missing package: {e}")
        all_passed = False
    
    # Test 3: AlphaGenome (if Module 4 selected)
    if 4 in modules:
        print("\n[3/7] Checking AlphaGenome setup...")
        if os.environ.get('ALPHA_GENOME_KEY'):
            print("  ✓ ALPHA_GENOME_KEY is set")
            try:
                from alphagenome.models import dna_client
                print("  ✓ alphagenome package importable")
            except ImportError:
                print("  ✗ alphagenome package not found")
                print("    Hint: conda activate alphagenome-env")
                all_passed = False
        else:
            print("  ✗ ALPHA_GENOME_KEY not set")
            print("    Hint: export ALPHA_GENOME_KEY='your_key'")
            all_passed = False
    else:
        print("\n[3/7] Skipping AlphaGenome check (Module 4 not selected)")
    
    # Test 4: Required scripts exist
    print("\n[4/7] Checking required scripts...")
    missing_scripts = []
    for mod_num in modules:
        module = MODULES.get(mod_num, {})
        for script_info in module.get('scripts', []):
            script_path = base_dir / script_info[1]
            if not script_path.exists():
                missing_scripts.append(str(script_path))
    
    if missing_scripts:
        print(f"  ✗ Missing scripts:")
        for s in missing_scripts:
            print(f"      {s}")
        all_passed = False
    else:
        print(f"  ✓ All {sum(len(MODULES[m].get('scripts', [])) for m in modules)} scripts found")
    
    # Test 5: Prerequisites for selected modules
    print("\n[5/7] Checking module prerequisites...")
    for mod_num in modules:
        success, msg = check_prerequisites(mod_num, base_dir, config)
        if success:
            print(f"  ✓ Module {mod_num}: {msg}")
        else:
            print(f"  ✗ Module {mod_num}: {msg}")
            all_passed = False
    
    # Test 6: Disk space
    print("\n[6/7] Checking disk space...")
    import shutil
    total, used, free = shutil.disk_usage(base_dir)
    free_gb = free / (1024**3)
    if free_gb > 10:
        print(f"  ✓ Free disk space: {free_gb:.1f} GB")
    else:
        print(f"  ⚠ Low disk space: {free_gb:.1f} GB (recommend >10 GB)")
    
    # Test 7: Strand fix validation
    print("\n[7/7] Validating strand-aware fix in Module 02...")
    enumerate_path = base_dir / '02_generate_mutation_paths/enumerate_paths.py'
    if enumerate_path.exists():
        content = enumerate_path.read_text()
        if 'reverse_complement' in content and 'ref_genomic' in content:
            print("  ✓ Strand-aware allele handling is implemented")
        else:
            print("  ✗ Strand fix NOT found in enumerate_paths.py")
            print("    The fix for Bug 1 may not be applied!")
            all_passed = False
    else:
        print(f"  ✗ Script not found: {enumerate_path}")
        all_passed = False
    
    # Summary
    print("\n" + "=" * 80)
    if all_passed:
        print("✓ ALL SMOKE TESTS PASSED")
        print("=" * 80)
        print("\nYou can now run the pipeline with:")
        modules_str = ','.join(map(str, modules))
        print(f"  python run_pipeline.py --config pipeline_config.yaml -m {modules_str}")
    else:
        print("✗ SOME TESTS FAILED - Please fix issues before running")
        print("=" * 80)
    
    return all_passed


# =============================================================================
# MAIN PIPELINE RUNNER
# =============================================================================

def run_pipeline(
    config_path: Path,
    modules: list,
    dry_run: bool = False,
    base_dir: Path = None,
    log_file: Path = None
) -> bool:
    """
    Run selected pipeline modules in order.
    
    Parameters
    ----------
    config_path : Path
        Path to pipeline config
    modules : list
        List of module numbers to run
    dry_run : bool
        If True, print commands without executing
    base_dir : Path
        Pipeline base directory
    log_file : Path
        If provided, capture all output to this file
        
    Returns
    -------
    bool
        True if all modules succeeded
    """
    if base_dir is None:
        base_dir = Path(__file__).parent
    
    config = load_config(config_path)
    
    # Sort modules
    modules = sorted(modules)
    
    print("\n" + "=" * 80)
    print("DORMANT SITE ACTIVATION PIPELINE")
    print("=" * 80)
    print(f"\nConfiguration: {config_path}")
    print(f"TF: {config.get('tf_name', 'AP1')}")
    print(f"Modules to run: {modules}")
    print(f"Mode: {'DRY RUN' if dry_run else 'EXECUTE'}")
    
    # Estimate total time
    total_estimate = ""
    estimates = []
    for mod_num in modules:
        module = MODULES.get(mod_num, {})
        estimates.append(module.get('estimated_time', 'unknown'))
    print(f"Estimated time: {' + '.join(estimates)}")
    
    print("\n" + "-" * 80)
    
    overall_start = time.time()
    failed_modules = []
    
    for mod_num in modules:
        module = MODULES.get(mod_num)
        if not module:
            print(f"\n⚠ Unknown module: {mod_num}, skipping")
            continue
        
        print(f"\n{'='*80}")
        print(f"MODULE {mod_num}: {module['name']}")
        print(f"{'='*80}")
        print(f"Description: {module['description']}")
        print(f"Estimated time: {module['estimated_time']}")
        
        # Check prerequisites
        success, msg = check_prerequisites(mod_num, base_dir, config)
        if not success:
            print(f"\n⚠ Prerequisites not met:\n  {msg}")
            if not dry_run:
                print("  Stopping pipeline.")
                failed_modules.append(mod_num)
                break
        
        # Check environment
        success, msg = check_environment(mod_num)
        if not success:
            print(f"\n⚠ Environment issue: {msg}")
            if not dry_run:
                print("  Stopping pipeline.")
                failed_modules.append(mod_num)
                break
        
        # Run scripts
        module_start = time.time()
        scripts = module.get('scripts', [])
        
        if not scripts:
            print(f"\n  No scripts to run for this module.")
            continue
        
        for i, script_info in enumerate(scripts, 1):
            interpreter = script_info[0]
            script = script_info[1]
            args = script_info[2] if len(script_info) > 2 else []
            
            # Replace {config} and {tf_name} placeholders
            tf_name = config.get('tf_name', 'AP1')
            args = [a.format(config=str(config_path), tf_name=tf_name) for a in args]
            
            print(f"\n  [{i}/{len(scripts)}] {script}")
            
            cmd = [interpreter, script] + args
            success = run_command(cmd, base_dir, dry_run=dry_run, log_file=log_file)
            
            if not success and not dry_run:
                print(f"\n  ✗ Script failed: {script}")
                failed_modules.append(mod_num)
                break
        
        if mod_num in failed_modules:
            break
        
        module_elapsed = time.time() - module_start
        print(f"\n  ✓ Module {mod_num} completed in {format_duration(module_elapsed)}")
    
    # Summary
    overall_elapsed = time.time() - overall_start
    
    print("\n" + "=" * 80)
    print("PIPELINE SUMMARY")
    print("=" * 80)
    print(f"Total time: {format_duration(overall_elapsed)}")
    
    if failed_modules:
        print(f"Failed modules: {failed_modules}")
        print("Status: ✗ FAILED")
        return False
    else:
        print(f"Completed modules: {modules}")
        print("Status: ✓ SUCCESS")
        return True


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Dormant Site Activation Pipeline Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run smoke test first
  python run_pipeline.py --config pipeline_config.yaml --smoke-test -m 2,3,4,5

  # Dry run to see what would execute
  python run_pipeline.py --config pipeline_config.yaml -m 2,3,4,5 --dry-run

  # Run modules 2-5 (after fixing strand bug)
  python run_pipeline.py --config pipeline_config.yaml -m 2,3,4,5

  # Run only landscape computation
  python run_pipeline.py --config pipeline_config.yaml -m 5
        """
    )
    
    parser.add_argument(
        '--config', '-c',
        type=Path,
        default=Path('pipeline_config.yaml'),
        help='Path to pipeline configuration file'
    )
    
    parser.add_argument(
        '--modules', '-m',
        type=str,
        default=None,
        help='Comma-separated list of modules to run (e.g., "2,3,4,5")'
    )
    
    parser.add_argument(
        '--smoke-test',
        action='store_true',
        help='Run smoke tests to validate setup before execution'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print commands without executing them'
    )
    
    parser.add_argument(
        '--log', '-l',
        action='store_true',
        help='Enable logging to a timestamped file in logs/ directory'
    )
    
    parser.add_argument(
        '--log-file',
        type=Path,
        default=None,
        help='Specify custom log file path (implies --log)'
    )
    
    parser.add_argument(
        '--list-modules',
        action='store_true',
        help='List available modules and exit'
    )
    
    args = parser.parse_args()
    
    base_dir = Path(__file__).parent
    
    # List modules
    if args.list_modules:
        print("\nAvailable Pipeline Modules:")
        print("-" * 60)
        for num, module in sorted(MODULES.items()):
            skip = " [skip by default]" if module.get('skip_by_default') else ""
            print(f"  {num}: {module['name']}{skip}")
            print(f"      {module['description']}")
            print(f"      Time: {module['estimated_time']}")
            print()
        return
    
    # Parse modules
    if args.modules:
        try:
            modules = [int(m.strip()) for m in args.modules.split(',')]
        except ValueError:
            print("Error: Modules must be comma-separated integers (e.g., '2,3,4,5')")
            sys.exit(1)
    else:
        # Default: modules 2-5 (skip 0 and 1 which are usually pre-computed)
        modules = [2, 3, 4, 5]
    
    # Validate module numbers
    for m in modules:
        if m not in MODULES:
            print(f"Error: Unknown module {m}. Valid modules: {list(MODULES.keys())}")
            sys.exit(1)
    
    # Check config exists
    config_path = base_dir / args.config if not args.config.is_absolute() else args.config
    if not config_path.exists():
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)
    
    # Run smoke test
    if args.smoke_test:
        success = run_smoke_test(base_dir, config_path, modules)
        if not success:
            sys.exit(1)
        if not args.dry_run:
            print("\nSmoke test passed. Add --dry-run to preview, or remove --smoke-test to execute.")
            return
    
    # Setup logging if requested
    log_file = None
    tee_logger = None
    
    if args.log_file:
        log_file = args.log_file
        args.log_file.parent.mkdir(parents=True, exist_ok=True)
    elif args.log:
        log_file = setup_logging(base_dir / 'logs', modules)
    
    if log_file and not args.dry_run:
        # Setup tee logging for print statements
        tee_logger = TeeLogger(log_file)
        sys.stdout = tee_logger
        print(f"\n📝 Logging enabled: {log_file}")
    
    try:
        # Run pipeline
        success = run_pipeline(
            config_path=config_path,
            modules=modules,
            dry_run=args.dry_run,
            base_dir=base_dir,
            log_file=log_file
        )
    finally:
        # Restore stdout and close log
        if tee_logger:
            sys.stdout = tee_logger.terminal
            tee_logger.close()
            print(f"\n📝 Log saved to: {log_file}")
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
