#!/usr/bin/env python

import logging
from protinfo._version import version_tuple
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
logger.info(f"Version :{version_tuple}")


NO_MCCE_MSG = """The mcce executable was not found.
The ProtInfo report will not include any information or diagnostics
from MCCE step1.py."""


USER_MCCE = shutil.which("mcce")
if USER_MCCE is None:
    print(NO_MCCE_MSG)
    logger.warning(NO_MCCE_MSG)


class Opts:
    def __init__(self, cli_name: str = "ProtInfo", **kwargs):
        """Class Opts makes command line options available to all objects."""

        self.cli_name = cli_name
        self.all = kwargs

    def __str__(self):
        out = f"{self.cli_name} - User options: "
        if not self.all:
            return out + "(not set)"

        for k in self.all:
            out = out + f"{k}: {self.all[k]!r}; "
        return out + "\n"

    def __repr__(self):
        return self.__str__()


cli_opts = Opts()
