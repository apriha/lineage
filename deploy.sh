#!/bin/sh

# Travis CI deployment script - timestamp version for test PyPI deployment

# Copyright (C) 2017 Andrew Riha
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# if TRAVIS_TAG environment variable is empty
# https://stackoverflow.com/a/3061064
if [ -z "$TRAVIS_TAG" ]; then
  # get current version
  # https://stackoverflow.com/a/2439587
  ver=$(head -n 1 lineage/VERSION);
  dev=".dev";

  # append .devNOW to version so that version is unique
  # https://stackoverflow.com/a/4181721
  echo ${ver}${dev}$(date +%Y%m%d%H%M%S) > lineage/VERSION;
fi
