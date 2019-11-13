#CSCI 580 Final Project - Curly Hair Simulation via Mass-Springs

[Original assignment](http://run.usc.edu/cs520-s19/assign1/)

##Controls
- `b` - Show/Hide Bend Springs
- `e` - Reset viewmode to Springs
- `h` - Show/Hide Shear Springs 
- `p` - Pause
- `s` - Show/Hide Structural Springs
- `v` - Switch viewmode between Springs and Jello
- `Right Click` - Camera Rotation
- `Space Bar` - Take 300 screenshots and save as picxxxx.PPM files (xxxx represents frame number). Program will terminate afterwards

##Relevant Files
###physics.h/cpp
- Majority of the implementation is here. 
- Responsible for calculating all spring, force-field, and collision forces.

###input.h/cpp
- Read/Write World files. World files contain information on the integrator type, timestep size, spring coefficients, force fields, and other information. See detailed info [here](http://run.usc.edu/cs520-s19/assign1/world.html).
- Responsible for defining keyboard actions

###jello.h/cpp
- Main program execution starts here.
- Sets up OpenGL and Window
