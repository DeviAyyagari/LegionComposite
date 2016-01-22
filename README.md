## Legion Compositor (Simpler Version)
A Legion-based Image Compositor that performs interactive rendering in an HPC environment.

## What's Different in this Version?
This version of the compositor is significantly stripped down compared to the master branch. Notably, all QT components have been removed, as well as isosurfacing and Optix rendering. Further, in order to run better on complex hybrid systems, all task-level recursion has been eliminated from the compositor and redone in a serialized format. However, for the moment all PhaseBarrier concurrency has also been removed.

## Dependencies
 * GNU C Compiler 4.7 or greater (Tested with 4.9.2) or equivalent compiler
 * Stanford Legion (https://github.com/StanfordLegion/legion) master branch
 * If using GPU rendering, NVIDIA CUDA (Only tested with 7.0)

## Compilation Instructions
Compilation of the compositor itself simply requires running the Makefile in the main project directory. Some components may require you to manually edit the Makefile itself to specify paths.
In order to specify where Legion is, you need to set the environment variable 'LG_RT_DIR' to point to the runtime directory of Legion.

 * Compile the main program with 'make' in the main legioncomposite folder. Make sure to set GASNET_ROOT.

Compiling the viewer requires you to be on a 64-bit system for the moment, and can be built in the CompositeViewer folder using the command
```bash
g++ -O2 -o viewer viewer.cc `Magick++-config --cppflags --cxxflags --ldflags --libs`
```
