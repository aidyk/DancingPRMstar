V-REP Stubs generator
=====================

This utility is used to generate boilerplate code for V-REP Lua callbacks.
It reads an XML file containing a description of the callbacks, script
functions, and enums, and it produces a pair of C++ source/header files.

What you need:
    - Python interpreter (2.7 or greater)
    - lxml package for Python (available from pip)
    - tempita package for Python (available from pip)

Usage:

    > python -m v_repStubsGen -H stubs.h -C stubs.cpp callbacks.xml

    generates cpp/h in one shot.
