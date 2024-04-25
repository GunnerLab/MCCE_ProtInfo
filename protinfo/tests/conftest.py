#!/usr/bin/env python

import pytest

def pytest_tempdir_basename(request):
    return "pinfo_tempdir"
