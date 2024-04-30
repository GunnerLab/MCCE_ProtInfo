#!/usr/bin/env python

import shutil
import warnings


NO_MCCE_MSG = """The mcce executable was not found.
The ProtInfo report will not include any information or diagnostics
from MCCE step1.py."""

USER_MCCE = shutil.which("mcce")
if USER_MCCE is None:
    warnings.warn(NO_MCCE_MSG)
