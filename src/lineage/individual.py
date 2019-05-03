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

from lineage.snps import SNPsCollection


class Individual(SNPsCollection):
    """ Object used to represent and interact with an individual.

    The ``Individual`` object maintains information about an individual. The object provides
    methods for loading an individual's genetic data (SNPs) and normalizing it for use with the
    `lineage` framework.

    """

    def __init__(self, name, raw_data=None, output_dir="output"):
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
        super().__init__(name=name, raw_data=raw_data, output_dir=output_dir)

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
        return self.get_var_repr(self.name)

    def remap_snps(self, target_assembly, complement_bases=True):
        """ Remap the SNP coordinates of this ``Individual`` from one assembly to another.

        This method is a wrapper for `remap_snps` in the ``Lineage`` class.

        This method uses the assembly map endpoint of the Ensembl REST API service to convert SNP
        coordinates / positions from one assembly to another. After remapping, the coordinates /
        positions for the ``Individual``'s SNPs will be that of the target assembly.

        If the SNPs are already mapped relative to the target assembly, remapping will not be
        performed.

        Parameters
        ----------
        target_assembly : {'NCBI36', 'GRCh37', 'GRCh38', 36, 37, 38}
            assembly to remap to
        complement_bases : bool
            complement bases when remapping SNPs to the minus strand

        Returns
        -------
        chromosomes_remapped : list of str
            chromosomes remapped; empty if None
        chromosomes_not_remapped : list of str
            chromosomes not remapped; empty if None

        Notes
        -----
        An assembly is also know as a "build." For example:

        Assembly NCBI36 = Build 36
        Assembly GRCh37 = Build 37
        Assembly GRCh38 = Build 38

        See https://www.ncbi.nlm.nih.gov/assembly for more information about assemblies and
        remapping.

        References
        ----------
        ..[1] Ensembl, Assembly Map Endpoint,
          http://rest.ensembl.org/documentation/info/assembly_map
        """
        from lineage import Lineage

        l = Lineage()
        return l.remap_snps(self, target_assembly, complement_bases)

    def _set_snps(self, snps, build=37):
        """ Set `_snps` and `_build` properties of this ``Individual``.

        Notes
        -----
        Intended to be used internally to `lineage`.

        Parameters
        ----------
        snps : pandas.DataFrame
            individual's genetic data normalized for use with `lineage`
        build : int
            build of this ``Individual``'s SNPs
        """
        self._snps = snps
        self.sort_snps()
        self._build = build
