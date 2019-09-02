February 6, 2019

the Johns Hopkins University Energise Team

The Julia codes come along with this README file are developed by the Johns Hopkins University Energise Team which is specific for the "small test case" in Q6. 

=======================
# Julia Package Version
-Julia 1.0.2
-JuMP 0.18.5
-Gurobi 0.5.6

=======================
# Code structure
-GML_main.jl  
	-GML_boundaries.jl
	-GML_optimization.jl
	-GML_output.jl
	-GML_strcut.jl
The GML_main.jl is the main code which take the demand, solar and price data as input and generate the optimal trajectory as output, and other codes are generated to assist the GML_main.jl achieve its function.

=======================
# Function GML_main(Pd, Pg_trend, Lambda, Beta, Alpha) is the main function of the GML layer. 
Input:
-Pd: real-time demand data. Each row represents a feeder, and each column represents a time slot. E.g., in our test case, it should be 12 feeders with 288 slots (every 5 minutes for 24 hours). Therefore Pd is a 12 by 288 matrix. 
-Pg_trend: solar data.  A row vector. E.g. in our test case we have 288 slots, therefore, Pg_trend is a 1 by 288 vector.
-Lambda: Real-time market price data. A row vector, same form as Pg_trend.
-Beta: solar penalty. A positive scalar.
-alpha: Ancillary market price data. A row vector, same form as Pg_trend.

Output:
A solar Excel file "solar.xls" contains information of optimal solar trajectory for FOL.
 	-column 1: timestamp.
	-column 2 to column 13: solar generation P_g at each feeder.
	-column 14 to column 25: reactive solar generation Q_g

A battery Excel file "battery.xls" contains information of optimal battery trajectory for FOL.
     -column 1: timestamp.
    -column 2 to column 13: power discharge R_f at each feeder
    -column 14 to column 25: reactive power discharge C_f at each feeder
    -column 25 to column 36: battery SOC B_f at each feeder

=======================
# Function GML_boundaries(Pd, Pg_trend, Lambda, Beta, Alpha) is the function generate the required boundaries.
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
