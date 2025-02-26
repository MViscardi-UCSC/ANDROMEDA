from loguru import logger
import sys

# Configure Loguru
# logger.remove()  # Remove the default configuration
# logger.add(sys.stderr, level="INFO", format="{time} {level} {message}")

# Example: Add a file handler if needed
# logger.add("andromeda.log", level="DEBUG", rotation="1 MB", retention="7 days")

# This makes the logger importable from other modules
log = logger
