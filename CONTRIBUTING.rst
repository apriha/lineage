Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

Bug reports
-----------

When `reporting a bug <https://github.com/apriha/lineage/issues>`_ please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Documentation improvements
--------------------------

``lineage`` could always use more documentation, whether as part of the official ``lineage``
docs, in docstrings, or even on the web in blog posts, articles, and such.

Feature requests and feedback
-----------------------------

The best way to send feedback is to file an issue at https://github.com/apriha/lineage/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions are welcome :)

Development
-----------

To set up ``lineage`` for local development:

1. Fork `lineage <https://github.com/apriha/lineage>`_
   (look for the "Fork" button).
2. Clone your fork locally::

    git clone git@github.com:your_name_here/lineage.git

3. Create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

4. Setup a development environment::

    pip install pipenv
    pipenv install --dev

5. When you're done making changes, run all the tests with::

    pipenv run pytest --cov=lineage tests

6. Commit your changes and push your branch to GitHub::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull request guidelines
```````````````````````

If you need some code review or feedback while you're developing the code, just make the pull
request.

For merging, you should:

1. Ensure tests pass.
2. Update documentation when there's new API, functionality, etc.
3. Add yourself to ``CONTRIBUTORS.rst``.
