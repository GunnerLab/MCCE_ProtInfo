#!/usr/bin/env python

from pathlib import Path
import shutil
import warnings

#................................................................................
APP_NAME = "mcce_protinfo"

NO_MCCE_MSG = """The mcce executable was not found.
The ProtInfo report will not include any information or diagnostics from MCCE step1.py."""

USER_MCCE = shutil.which("mcce")
if USER_MCCE is None:
    warnings.warn(NO_MCCE_MSG)
