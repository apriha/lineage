""" Class for representing individuals within the `lineage` framework. """

"""
Copyright (C) 2016 Andrew Riha

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

from snps import SNPsCollection
from snps.utils import clean_str


class Individual(SNPsCollection):
    """ Object used to represent and interact with an individual.

    The ``Individual`` object maintains information about an individual. The object provides
    methods for loading an individual's genetic data (SNPs) and normalizing it for use with the
    `lineage` framework.

    """

    def __init__(self, name, raw_data=None, output_dir="output", **kwargs):
        """ Initialize an ``Individual`` object.

        Parameters
        ----------
        name : str
            name of the individual
        raw_data : list or str
            path(s) to file(s) with raw genotype data
        output_dir : str
            path to output directory
        """
        self._name = name
        super().__init__(name=name, raw_data=raw_data, output_dir=output_dir, **kwargs)

    def __repr__(self):
        return "Individual({!r})".format(self._name)

    @property
    def name(self):
        """ Get this ``Individual``'s name.

        Returns
        -------
        str
        """
        return self._name

    def get_var_name(self):
        return clean_str(self.name)
