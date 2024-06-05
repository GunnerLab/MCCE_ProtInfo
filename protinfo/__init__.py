#!/usr/bin/env python

from protinfo import _version
import logging
import shutil


# Config for root logger:
DT_FMT = "%Y-%m-%d %H:%M:%S"
BODY = "[%(levelname)s]: %(name)s, %(funcName)s:\n\t%(message)s"
logging.basicConfig(
    level=logging.INFO,
    format=BODY,
    datefmt=DT_FMT,
    filename="protinfo.log",
    encoding="utf-8",
)
logger = logging.getLogger("ProtInfo")
logger.info(f"Version :{_version.version_tuple}")


NO_MCCE_MSG = """The mcce executable was not found.
The ProtInfo report will not include any information or diagnostics
from MCCE step1.py."""

USER_MCCE = shutil.which("mcce")
if USER_MCCE is None:
    print(NO_MCCE_MSG)
    logger.warning(NO_MCCE_MSG)
