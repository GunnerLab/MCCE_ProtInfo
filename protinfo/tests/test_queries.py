#!/usr/bin/env python

from protinfo.queries import get_rcsb_pdb
import os
from pathlib import Path
import pytest
import tempfile


def test_create_file(tmp_path):
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "hello.txt"
    p.write_text(CONTENT)
    assert p.read_text() == CONTENT
    assert len(list(tmp_path.iterdir())) == 1
    assert 0



class TestGetRcsbPdb:

    def __init__(self, tmp_path):
        self.tmp_path = tmp_path
        self.here = Path.cwd()
        
    def test_valid_pdbid(self):
        """Given a valid pdbid, the function should download the pdb
        file containing the biological assembly from rcsb.org and
        return a Path object."""

        pdbid = "1abc"
        d = self.tmp_path/pdbid
        os.chdir(d)
        result = get_rcsb_pdb(pdbid)
        os.chdir(self.here)

        assert isinstance(result, Path)

    def test_invalid_pdbid(self):
        """Given an invalid pdbid, the function should return None."""

        pdbid = "xyz"
        d = self.tmp_path/pdbid
        os.chdir(d)
        result = get_rcsb_pdb(pdbid)
        os.chdir(self.here)

        assert result is None

    def test_overwrite_existing_file(self):
        """Given a pdbid that already exists in the current directory,
        the function should overwrite the existing file with the new download.
        """

        pdbid = "1abc"
        d = self.tmp_path/pdbid
        os.chdir(d)
        existing_file = Path("1abc.pdb").resolve()
        existing_file.touch()

        result = get_rcsb_pdb(pdbid)
        os.chdir(self.here)

        assert isinstance(result, Path)
        assert result == existing_file
        assert existing_file.exists()
