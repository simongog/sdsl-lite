Biicode C/C++ dependency manager
=================================

New with biicode? Check the [Getting Started Guide](http://docs.biicode.com/c++/gettingstarted.html).

How to build it?
------------------
Building it is too easy:

    $ git clone git@github.com:simongog/sdsl-lite.git
    $ cd sdsl-lite
    $ bii init -l && bii work

It creates all the necessary structure to build sdsl-lite with biicode build system. Then, you can build its default samples:

    $ bii build
    $ ./bin/any_executable

By default, the first use applies all the changes to the repository, if you want to revert these ones, set the `BII_SDSL_LITE_REVERT_CHANGES` environment variable to `True` and run `bii work` to keep your original code and undo the biicode changes.

    $ export BII_SDSL_LITE_REVERT_CHANGES=True
    $ bii work



How to use it in other external projects?
-------------------------------------------

Take a look at any example from [the examples/sdsl-lite block](https://www.biicode.com/examples/sdsl-lite) and try your own using *sdsl-lite library* with biicode.

