January 27, 2019

This is version 1.0 of the "small test case" for the experiments to be performed in Q6.

There are two top-level directories: GML (for JHU) and FOL (for PNNL and UVM).


GML
===

This directory contains 6 minutely demand and solar csv files corresponding to the transformer banks of Allendale (4: 2 for February 2016 and 2 for August 2016) and Blooming Grove (2: 1 for February 2016 and 1 for August 2016). The overall format of each file is more or less the same as from Q5:

column 1 => timestamp.
column 2 => overall kW consumption of the bank, summed across phases 1, 2, and 3.
column 3 => overall "beta" of the bank (based on the power factor), such that S = P + j*beta*P, i.e., Q = beta*P.
column 4 => ckt1_kW => kW consumption of circuit 1, summed across phases 1, 2, and 3.
column 5 => ckt1_asum => Amperes drawn by circuit 1, summed across phases 1, 2, and 3.
column 6 => ckt1_fr1 => the fraction of kW or Amperes consumed by phase 1 of circuit 1.
column 7 => ckt1_fr2 => the fraction of kW or Amperes consumed by phase 2 of circuit 1.
column 8 => ckt1_fr3 => the fraction of kW or Amperes consumed by phase 3 of circuit 1.
columns 9-13 => columns 4-8 for circuit 2.
columns 14-18 => columns 4-8 for circuit 3.
columns 19-23 => columns 4-8 for circuit 4.
column 24 => panel_output_kW => the output in kW of a 5 kW solar array, based on the irradiance data provided by Clean Power Research.

The seventh file in this directory is "solar_installed_capacity.txt". It contains the installed capacity (column: capacity_MW) and number of 5 kW solar arrays (column: number_of_5kW_arrays) needed to meet 50% of the annual energy demand for each circuit. Note that we use a capacity factor of 15% (which is a reasonable assumption for the NY/NJ area) to arrive at these number.

The net demand for each circuit (N) and timestamp, consistent with that used by the FOL, is:

net_demand_kW = ckt{N}_kW*0.95 - panel_output_kW*number_of_5kW_arrays

Where 0.95 = (1-0.05), corresponding to 5% system-aggregate losses.

net_demand_kVAr = bnk_beta*net_demand_kW


FOL/allendale/circuit_1
=======================

It has 2 subdirectories: "full" and "kron".

The latter directory, i.e., kron, which contains 101 nodes (51 "super-nodes" and 50 "super-loads"), is what users should use for Q6. Under this directory, there is a file named "graph.pdf", which contains a layout of the Kron-reduced Allendale circuit 1. Note that not shown are the 50 super-loads connected to super-nodes 1-50 via nominal impedances, i.e., directly at the medium voltage. Also, this is a radialized version of the actual Kron-reduced graph, which is weakly-meshed.

MODEL.glm is the GridLab-D markup for the circuit. Under "players/feb_2016" and "players/aug_2016", there are GLD-ready net demand csv files corresponding to the 50%-annual-energy-demand-met-by-solar scenario.

I've successfully run the above scenarios in GLD. A week at a minutely timestep takes about 45 seconds to complete on my somewhat dated ThinkPad laptop. I can separtely provide the outputs for these simulations as a point of reference.


That's all I can think of for now. Should you have any questions, please do not hesitate to email me at pracherl@uvm.edu.

Thanks,
Pavan
