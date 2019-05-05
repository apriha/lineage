from multiprocessing import Pool
import os


class Parallelizer:
    def __init__(self, parallelize=True, processes=os.cpu_count()):
        """ Initialize a `Parallelizer`.

        Parameters
        ----------
        parallelize : bool
            utilize multiprocessing to speedup calculations
        processes : int
            processes to launch if multiprocessing
        """
        self._parallelize = parallelize
        self._processes = processes

    def __call__(self, f, tasks):
        """ Optionally parallelize execution of a function.

        Parameters
        ----------
        f : func
            function to execute
        tasks : list of dict
            tasks to pass to `f`

        Returns
        -------
        list
            results of each call to `f`
        """
        if self._parallelize:
            with Pool(self._processes) as p:
                return p.map(f, tasks)
        else:
            return map(f, tasks)
