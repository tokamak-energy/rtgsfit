import os
import shutil
import subprocess

import numpy as np
import pytest

from rtgsfit_vs_gsfit import cnst

def test_rtgsfit_mds_nodegen(rtgsfit_mds_nodegen):
    """
    Test the RTGSFIT MDSplus node generation.
    """
    assert 1 == 1