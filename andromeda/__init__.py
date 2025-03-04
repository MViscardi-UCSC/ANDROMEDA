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
for program in ["samtools", "bcftools"]:
    if not shutil.which(program):
        log.warning(
            f"{program} not found in PATH, please install it to use the BAM tools module."
        )