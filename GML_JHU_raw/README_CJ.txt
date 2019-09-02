February 6, 2019

the Johns Hopkins University Energise Team

The Julia codes that come along with this README file were developed by the Johns Hopkins University Energise Team for the "small test case" deliverables in Q6. 

=======================
# Software Versions
-Julia 1.0.2
-JuMP 0.18.5
-Gurobi 0.5.6

=======================
# Code structure
-GML.jl  
The GML.jl is the main file which takes the demand, solar, real-time real-time price, solar (curtailment) penalty and ancillary market price data as input data and generates the optimal trajectory of solar and battery usages at feeder level as output.
-GML_struct.jl
The GML_struct.jl is the main file contains all the required Julia constructor. 

=======================
# Required Julia Package

- CSV 0.4.3
- LinearAlgebra (Julia 1.0.2 built-in package) # We use function from that package dot() and norm() in our code.
- DataFrames 0.15.0

=======================
# Function GML(Pd, Pg_trend, Lambda, Beta, Alpha) is the main function of the GML layer. 
Input:
-Pd: real-time demand data. Each column represents a feeder, whereas and each row represents a time slot, e.g., in our test case, it should be 12 feeders with 288 slots (every 5 minutes for 24 hours). Therefore Pd is a 288 by 12 matrix. 
-Pg_trend: solar data.  A vector. E.g., in our test case we have 288 slots, therefore, Pg_trend a vector of length 288. 
-Lambda: Real-time market price data. A vector, same form as Pg_trend.
-Beta: Solar (curtailment) penalty. A positive scalar.
-alpha: Ancillary market price data. A column vector, same form as Pg_trend.

Output:
A solar CSV file "solar.csv" comprising the optimal solar trajectory for the FOL.
 	-column 1: timestamp.
	-column 2 to column 13: optimal solar generation P_g at each feeder.
	-column 14 to column 25: optimal reactive solar generation Q_g

A battery CSV file "battery.csv" comprising the optimal battery trajectory for the FOL.
    -column 1: timestamp.
    -column 2 to column 13: power discharge R_f at each feeder
    -column 14 to column 25: reactive power discharge C_f at each feeder
    -column 25 to column 36: battery SOC B_f at each feeder

=======================
# Function GML_boundaries(Pd, Pg_trend, Lambda, Beta, Alpha) is the function that returns the required boundaries.
Input:
Same as GML();
Output:
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
Assumption on boundaries:
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
# Function GML_optimization(bd) takes the boundaries as input and returns optimization variables. 
Input: 
-bd: same as the output of GML_boundaries
Output:
-optimized_value: composed by constructor "Optimization_output_struct" and it contains the following information:
    # solar generation Pg;
    # reactive solar generation Qg;
    # battery SOC  B;
    # battery power discharge  R;
    # battery reactive power discharge C;
    # branch power flow P_hat;
    # branch reactive power flow Q_hat;
    # current square l;
    # voltage square v;
    # schedule ancillary power P_rsrv;
    # schedule ancillary energy B_rsrv;
    # power consuming at head node P_0;
    # operating cost cost;

# Function GML_output(Opt_value) takes optimization variable as input and write the two CSV files.

Input:
-Opt_value: same as the output of GML_optimization





