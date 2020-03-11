# This runs all of the jobs necessary for our trade-off analysis

include("GML.jl")

# # of timeslots
T=288;
# # of banks
BN=3;
# # of feeders
F=BN*4;
# # of secnarios
SN=1;

################################################################################

# Ancillar Markets Considered: "without", "10min", "30min"
ancillary_type = "10min"
# penetation levels
# p_rate = [0.25 0.5 0.75, 1]; #[0.25, 0.5, 0.75, 1];
p_rate = 0.25;

# [90.95,99]
# icdf = [-1.2816, -1.6449, -2.3263]
icdf = -1.2816;
# prediction length [2hr, 1hr, .5hr]
pred_length = 24;

# max solar error
solar_error_max = 0.25;

# max solar error
B_cap = 3;
################################################################################

# load price data
price_raw = read_price_data()
delta_rt_raw=matread("../data/price_prediction.mat");

# read demand data
pd_raw = read_demand_data()
pd_noise = matread("../data/demand_noise.mat")["demand_noise"];

# read pg noise
pg_noise = matread("../data/solar_noise_0025.mat")["solar_noise"];

# read pg_raw
pg_raw = read_solar_data()



folder="results_single";
filename="testrun.csv";

GML(ancillary_type, T, BN, F, SN, p_rate, icdf,
     pred_length, solar_error_max, B_cap,
     price_raw, delta_rt_raw, pd_raw, pd_noise, pg_noise, pg_raw,
     folder, filename);
