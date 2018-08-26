"""
Copyright (C) 2017 Andrew Riha

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import os
import shutil

import pytest

from lineage import Lineage


def del_output_dir_helper():
    if os.path.exists('output'):
        shutil.rmtree('output')


@pytest.fixture(scope='module')
def l():
    return Lineage()


@pytest.fixture(autouse=True)
def del_output_dir():
    """ Delete output directory if it exists during setup / teardown. """
    del_output_dir_helper()
    yield
    del_output_dir_helper()
