""" Class for representing individuals within the `lineage` framework. """

"""
MIT License

Copyright (c) 2016 Andrew Riha

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
import inspect

from snps import SNPs
from snps.utils import clean_str


class Individual(SNPs):
    """Object used to represent and interact with an individual.

    The ``Individual`` object maintains information about an individual. The object provides
    methods for loading an individual's genetic data (SNPs) and normalizing it for use with the
    `lineage` framework.

    ``Individual`` inherits from ``snps.SNPs``. See here for details about the ``SNPs`` object:
    https://snps.readthedocs.io/en/latest/snps.html
    """

    def __init__(self, name, raw_data=(), **kwargs):
        """Initialize an ``Individual`` object.

        Parameters
        ----------
        name : str
            name of the individual
        raw_data : str, bytes, ``SNPs`` (or list or tuple thereof)
            path(s) to file(s), bytes, or ``SNPs`` object(s) with raw genotype data
        **kwargs
            parameters to ``snps.SNPs`` and/or ``snps.SNPs.merge``
        """
        self._name = name

        init_args = self._get_defined_kwargs(SNPs, kwargs)
        merge_args = self._get_defined_kwargs(SNPs.merge, kwargs)

        super().__init__(**init_args)

        # load raw data by merging `SNPs` objects into this object
        if not isinstance(raw_data, list) and not isinstance(raw_data, tuple):
            s = (
                SNPs(raw_data, **init_args)
                if not isinstance(raw_data, SNPs)
                else raw_data
            )
            self.merge([s], **merge_args)
        else:
            for file in raw_data:
                s = file
                if not isinstance(file, SNPs):
                    s = SNPs(file, **init_args)

                self.merge([s], **merge_args)

    def _get_defined_kwargs(self, callable, kwargs):
        sig = inspect.signature(callable)
        return {k: kwargs[k] for k in kwargs if k in sig.parameters}

    def __repr__(self):
        return "Individual({!r})".format(self._name)

    @property
    def name(self):
        """Get this ``Individual``'s name.

        Returns
        -------
        str
        """
        return self._name

    def get_var_name(self):
        return clean_str(self.name)
