"""
__init__.py
Marcus Viscardi,    March 04, 2025

Initialization file for the andromeda package.
"""
from andromeda.logger import log
import shutil

# Check if the user has the necessary programs installed
# First umi_tools:
if not shutil.which("umi_tools"):
    log.warning(
        "umi_tools not found in PATH, please install it to use the UMI grouping module."
    )
if not shutil.which("samtools"):
    log.critical(
        f"samtools not found in PATH, most (all?) modules will be non-functional."
    )