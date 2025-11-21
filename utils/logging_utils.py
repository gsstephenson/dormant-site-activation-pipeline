"""
Logging utilities for Dormant Site Activation Pipeline
"""

import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional


def setup_logger(
    name: str,
    log_file: Optional[str] = None,
    level: str = "INFO",
    format_string: Optional[str] = None
) -> logging.Logger:
    """
    Set up a logger with file and console handlers.
    
    Parameters
    ----------
    name : str
        Name of the logger (usually __name__)
    log_file : str, optional
        Path to log file. If None, only console logging is enabled.
    level : str
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    format_string : str, optional
        Custom format string for log messages
        
    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))
    
    # Remove existing handlers
    logger.handlers = []
    
    # Default format
    if format_string is None:
        format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    
    formatter = logging.Formatter(format_string)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def log_step(logger: logging.Logger, step_name: str, start: bool = True):
    """
    Log the start or end of a pipeline step with a separator.
    
    Parameters
    ----------
    logger : logging.Logger
        Logger instance
    step_name : str
        Name of the pipeline step
    start : bool
        If True, log step start. If False, log step completion.
    """
    separator = "=" * 80
    if start:
        logger.info(separator)
        logger.info(f"Starting: {step_name}")
        logger.info(separator)
    else:
        logger.info(f"Completed: {step_name}")
        logger.info(separator)


def log_parameters(logger: logging.Logger, params: dict):
    """
    Log a dictionary of parameters.
    
    Parameters
    ----------
    logger : logging.Logger
        Logger instance
    params : dict
        Dictionary of parameter names and values
    """
    logger.info("Parameters:")
    for key, value in params.items():
        logger.info(f"  {key}: {value}")


def log_file_info(logger: logging.Logger, file_path: str, description: str = "File"):
    """
    Log information about a file (existence, size).
    
    Parameters
    ----------
    logger : logging.Logger
        Logger instance
    file_path : str
        Path to file
    description : str
        Description of the file
    """
    path = Path(file_path)
    if path.exists():
        size_mb = path.stat().st_size / (1024 * 1024)
        logger.info(f"{description}: {file_path} ({size_mb:.2f} MB)")
    else:
        logger.warning(f"{description}: {file_path} (NOT FOUND)")


def create_timestamped_log(base_name: str, log_dir: str = "logs") -> str:
    """
    Create a timestamped log file path.
    
    Parameters
    ----------
    base_name : str
        Base name for the log file
    log_dir : str
        Directory for log files
        
    Returns
    -------
    str
        Path to timestamped log file
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir_path = Path(log_dir)
    log_dir_path.mkdir(parents=True, exist_ok=True)
    
    log_file = log_dir_path / f"{base_name}_{timestamp}.log"
    return str(log_file)


class ProgressLogger:
    """
    Log progress through a large iteration with periodic updates.
    """
    
    def __init__(self, logger: logging.Logger, total: int, step: int = 1000):
        """
        Parameters
        ----------
        logger : logging.Logger
            Logger instance
        total : int
            Total number of items to process
        step : int
            Log every N items
        """
        self.logger = logger
        self.total = total
        self.step = step
        self.current = 0
        self.start_time = datetime.now()
    
    def update(self, n: int = 1):
        """
        Update progress counter.
        
        Parameters
        ----------
        n : int
            Number of items processed
        """
        self.current += n
        
        if self.current % self.step == 0 or self.current == self.total:
            elapsed = (datetime.now() - self.start_time).total_seconds()
            rate = self.current / elapsed if elapsed > 0 else 0
            remaining = (self.total - self.current) / rate if rate > 0 else 0
            
            percent = 100 * self.current / self.total
            self.logger.info(
                f"Progress: {self.current:,}/{self.total:,} ({percent:.1f}%) | "
                f"Rate: {rate:.1f} items/sec | "
                f"Remaining: {remaining/60:.1f} min"
            )
