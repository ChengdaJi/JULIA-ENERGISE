February 6, 2019

the Johns Hopkins University Energise Team

The Julia codes that come along with this README file were developed by the Johns Hopkins University Energise Team for the "small test case" deliverables in Q6.

=======================
# Software Versions
-Julia 1.0.2
-JuMP 0.18.5
-Gurobi 0.5.6

=======================
# Required Julia Packages

- CSV 0.4.3
- LinearAlgebra (Julia 1.0.2 built-in package) # We use functions dot() and norm() from that package in our code.
- DataFrames 0.15.0

=======================
# Code structure
-GML.jl
The GML.jl is the main file which takes the demand, solar, real-time real-time price, solar (curtailment) penalty and ancillary market price data as input data and generates the optimal trajectory of solar and battery dispatch for each feeder as output. It contains all the function listed below.

-GML_struct.jl
The GML_struct.jl is the main file contains all the required Julia constructor.

-example.jl
This file is a example of utilizing the GML.jl and GML_struct.jl file.

=======================
# GML.jl
# Function GML(Pd, Pg_trend, Lambda, Beta, Alpha) is the main function of the GML layer.
Input:
-Pd: real-time demand data. Each row represents a feeder whereas each column represents a time slot, e.g., in the small test case, there are 12 feeders with 288 slots (every 5 minutes for 24 hours). Therefore Pd is a 12 by 288 matrix.
-Sg_max: solar data. Each row represents a feeder whereas each column represents a time slot, e.g., in the small test case, there are 12 feeders with 288 slots (every 5 minutes for 24 hours). Therefore Pd is a 12 by 288 matrix.
-Lambda: Real-time market price data. A vector of size 288.
-Beta: Solar (curtailment) penalty. A positive scalar.
-alpha: Ancillary market price data. A column vector, same form as Lambda.

Output:
- Pf power at each feeders
 - size 12 x 288 (12 feeders and 288 time slots)
- Power.csv a csv file contains Pf

=======================
# example.jl
it is an example of applying the GML() function. All the required data are in the director ~/data.










=====================================================================
=====================================================================
The detailed sub-functions' descriptions are listed below.
=====================================================================
=====================================================================

=======================
# Function GML_boundaries(Pd, Pg_trend, Lambda, Beta, Alpha) is a sub-function of GML() which generate boundaries for Pg, Pd, R, C, B, and generate data for W.
Inputs:
Same as GML();
Outputs:
-bd: composed by constructor "bd_struct" and it contains the following information:
	#number of time slots T;
        #number of feeders F;
        #number of banks BN;
        # demand Pd;
        # reactive demand Qd;
        # maximum power output Pg_max;
        # minimum power output Pg_min;
        # maximum reactive power output Qg_max;
        # minimum reactive power output Qg_min;
        # the maximum power of storage B_max;
        # the minimum power of storage B_min;
        # the maximum power discharge rate R_max;
        # the minimum power discharge rate R_min;
        # power storage C_max;
        # the minimum power of storage C_min;
        # power storage b_0;
        # the exogenous change  W;
        # time interval delta_t;
        # The maximum ancillary power P_rsrv_min;
        # the apparent power S;
        # the maximum voltage V_max;
        # the minimum voltage V_min;
        # price of ancillary power alpha;
        # price of cost beta;
        # the grid spot price lambda;
        # ancillary respond time tau;
        # impedance r;
        # impedance  x;
        # time asso/w ancellary k;
Assumption on bounds:
 	Pg_max = 0.5 max(Pd)*Pg_trends/sum(Pg_trends)*sum(Pd)
	Qg_max = 0.05*Pg_max
	Pg_min = 0;
	Qg_min = -Qg_max;
	B_max = 0.1*Pd
	B_min=0;
	R_max = B_max;
	R_min = - B_max;
	C_max = 0.05*R_max;
	C_min = -0.05*R_max;
	B_f(0)=0.5*B_max_f(0);
	B_f(t_f)=0.5*B_max_f(t_f);	(t_f is the last slot)
	if B_max(t)-B_max(t-1)>0
		W_f_t=0.1*(B_max(t)-B_max(t-1))
	otherwise
		W_f_t=0.9*(B_max(t)-B_max(t-1))
	P_rsrv_min=0;

=======================
# Function GML_optimization(bd) is a sub-function of GML() and it takes the upper and lower bounds on the power and energy as inputs and returns optimization variables.
Input:
-bd: same as the output of GML_boundaries
Output:
-optimized_value: composed by constructor "Optimization_output_struct" and it contains the following information:
    # solar generation Pg; size: 12* 288
    # reactive solar generation Qg; size: 12* 288
    # battery SOC  B; size: 12* 288
    # battery power discharge  R; size: 12* 288
    # battery reactive power discharge C; size: 12* 288
    # branch power flow P_hat; size: 3 * 288
    # branch reactive power flow Q_hat; size: 3* 288
    # current square l; size: 3* 288
    # voltage square v; size: 3* 288
    # schedule ancillary power P_rsrv; 288
    # schedule ancillary energy B_rsrv; 288
    # power consuming at head node P_0; 288
    # operating cost cost; 1

=======================
# Function GML_output(Opt_value) is a sub-function of GML() and it takes optimization variable as input and writes the Power.CSV files that are outputs of the function GML().

Input:
-Opt_value: same as the output of GML_optimization
