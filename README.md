# DancingPRMstar

An implementation of DancingPRMstar. Current version requires a lot of effort to refine the overall performance / readability also
might have lots of bugs. Please enjoy it as an underdeveoped implementation.

How to run
Installation tips:
- OMPL(Open Motion Planning Library, https://ompl.kavrakilab.org/, tested on 1.1.0)
- V-REP(http://www.coppeliarobotics.com/, tested on 3.2.2 EDU)
- It is based on Lazy PRM\*(Lazy Collision Checking in Asymptotically-Optimal Motion Planning, Kris Hauser, 2015) with my own implemenetation of DSPT(Fully dynamic algorithms for maintaining shortest paths trees, D. Frigioni, Journal of Algoritms, 2000).
- Trajectory optimizer is not provided in this version. Please refer to the corresponding source code
for the implementation tip. You can use uBLAS or Eigen library for linear algebra computation; this version
is Eigen-friendly though.
- Routines for experiment and debug might remain incomplete/unused.

1. Copy "DancingPRMstar.cpp" and "DancingPRMstar.h" into *$(OMPL_PATH)/src/ompl/geometric/planners/prm/*

Successful compilation would generate a library file(e.g., libompl.\*).
Copy it into */usr/local/lib* so that V-REP can load the ompl library file successfully.

2. Replace v_repExtOMPL directory in V-REP directory(*$(VREP_PATH)/programming/v_repExtOMPL*) with 
downloaded v_repExtOMPL.

3. Compile v_repExtOMPL and copy libv.repExtOMPL.\* into V-REP directory(*vrep.app/Contents/MacOS/*)

4. Run V-REP and check the console ouptut to make sure that all the library files are successfully loaded.

5. Open example_scene.ttt in V-REP and Run!
![img](https://github.com/aidyk/DancingPRMstar/blob/master/example_screenshot.png)

If you have any problem with installation or following the above steps, feel free to mail me or leave a comment.
