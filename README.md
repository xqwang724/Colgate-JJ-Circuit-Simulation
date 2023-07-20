# ColgateJJSim

Here we present the source codes for simulating classical superconducting circuits. Contrary to existing methods that use differential algebraic equations, this program simulates superconducting circuits using only differential equations. These differential equations can be solved using ```SciPy.integrate.solve_ivp```. More theoretical details and demonstrations can be found in the following paper:

> Simulation of JJ circuits with only ODEs.



## User guide

To use the program, put the SPICE circuit file in the same folder as the Python scripts. Change the name of the circuit accordingly in ```cir_to_ivp.py```.

```python
ColgateJJSim(circuit_name, simulation_length)
```



## Structure of the code

The program is separated into four python scripts: ```read_cir.py``` , ```findtree.py```, ```construct_matrices.py```, and ```cir_to_ivp.py```.

```read_cir.py``` converts the SPICE circuit files into a directed multigraph using NetworkX. The graph contains information about the component type and model information of the Josephson junction. Currently, the only quasiparticle resistance model we include is the zero shunt conductance (rtype=0). ```findtree.py``` will separate the graph into two subgraphs: tree and chords. The partitioned graph is then converted into fundamental loop matrices and matrix equations in ```construct_matrtices.py```. Finally, the simulation is run in the ```cir_to_ivp.py```. The outputs included in the source code are junction phases and voltages, capacitor phases and voltages, and resistor phases. Other outputs can be acquired from these five sets of outputs post-simulation through algebra. 
