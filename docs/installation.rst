Installation
============

``lineage`` is `available <https://pypi.org/project/lineage/>`_ on the
`Python Package Index <https://pypi.org>`_. Install ``lineage`` (and its required
Python dependencies) via ``pip``::

    $ pip install lineage

Installation and Usage on a Raspberry Pi
----------------------------------------
The instructions below provide the steps to install ``lineage`` on a
`Raspberry Pi <https://www.raspberrypi.org>`_ (tested with
"`Raspberry Pi OS <https://www.raspberrypi.org/downloads/raspberry-pi-os/>`_ (32-bit) Lite",
release date 2020-08-20). For more details about Python on the Raspberry Pi, see
`here <https://www.raspberrypi.org/documentation/linux/software/python.md>`_.

.. note:: Text after a prompt (e.g., ``$``) is the command to type at the command line. The
          instructions assume a fresh install of Raspberry Pi OS and that after logging in as
          the ``pi`` user, the current working directory is ``/home/pi``.

1. Install ``pip`` for Python 3::

    pi@raspberrypi:~ $ sudo apt install python3-pip

   Press "y" followed by "enter" to continue. This enables us to install packages from the
   Python Package Index.

2. Install the ``venv`` module::

    pi@raspberrypi:~ $ sudo apt install python3-venv

   Press "y" followed by "enter" to continue. This enables us to create a
   `virtual environment <https://docs.python.org/3/library/venv.html>`_ to isolate the ``lineage``
   installation from other system Python packages.

3. `Install ATLAS <https://github.com/Kitt-AI/snowboy/issues/262#issuecomment-324997127>`_::

    pi@raspberrypi:~ $ sudo apt install libatlas-base-dev

   Press "y" followed by "enter" to continue. This is required for `NumPy <https://numpy.org>`_, a
   dependency of ``lineage``.

4. `Install Pillow dependencies <https://www.piwheels.org/project/Pillow/>`_::

    pi@raspberrypi:~ $ sudo apt install libjbig0 liblcms2-2 libopenjp2-7 libtiff5 libwebp6 libwebpdemux2 libwebpmux3

   Press "y" followed by "enter" to continue. This is required for
   `Matplotlib <https://matplotlib.org>`_, a dependency of ``lineage``.

5. Create a directory for ``lineage`` and change working directory::

    pi@raspberrypi:~ $ mkdir lineage
    pi@raspberrypi:~ $ cd lineage

6. Create a virtual environment for ``lineage``::

    pi@raspberrypi:~/lineage $ python3 -m venv .venv

   The virtual environment is located at ``/home/pi/lineage/.venv``.

7. Activate the virtual environment::

    pi@raspberrypi:~/lineage $ source .venv/bin/activate

   Now when you invoke Python or ``pip``, the virtual environment's version will be used (as
   indicated by the ``(.venv)`` before the prompt). This can be verified as follows::

    (.venv) pi@raspberrypi:~/lineage $ which python
    /home/pi/lineage/.venv/bin/python

8. Install ``lineage``::

    (.venv) pi@raspberrypi:~/lineage $ pip install lineage

9. Start Python::

    (.venv) pi@raspberrypi:~/lineage $ python
    Python 3.7.3 (default, Jul 25 2020, 13:03:44)
    [GCC 8.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>>

10. Use ``lineage``; examples shown in the README should now work.

11. At completion of usage, the virtual environment can be deactivated::

     (.venv) pi@raspberrypi:~/lineage $ deactivate
     pi@raspberrypi:~/lineage $

Installation on Linux
---------------------
On Linux systems, the following system-level installs may also be required::

    $ sudo apt install python3-tk
    $ sudo apt install gfortran
    $ sudo apt install python-dev
    $ sudo apt install python-devel
    $ sudo apt install python3.X-dev # (where X == Python minor version)
