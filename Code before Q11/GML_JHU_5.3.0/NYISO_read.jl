    function shunt_data(data_trace_bus,data_trace_shunt)
    # data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
    baseKV_bus = collect(data_trace_bus[:,Symbol("baseKV")]);
    # data_trace_shunt=CSV.File("../data/NYISO-data/shunt.csv") |> DataFrame
    P = collect(data_trace_shunt[:,Symbol("P")]).+0.01;
    find_bus = collect(data_trace_shunt[:,Symbol("bus")]);
    Q_max = collect(data_trace_shunt[:,Symbol("Qmax")]);
    kind = collect(data_trace_shunt[:,Symbol("kind")]);
    baseKV = zeros(209,1)
    baseMVA = zeros(209,1)
    type = zeros(209,1)

    for shunt = 1:length(find_bus)
        baseKV[shunt]=baseKV_bus[find_bus[shunt]];
        baseMVA[shunt]=(baseKV[shunt]/20)^2;
        if kind[shunt]=="LOAD"
            type[shunt]=1;
        elseif kind[shunt]=="CTRL_V"
            type[shunt]=2;
        else
            type[shunt]=3;
        end
    end

    return shunt = (find_bus = (find_bus), baseKV=(baseKV),baseMVA=(baseMVA),
    Q_max=(Q_max), P=(P), type=(type))
end

function bus_data()
    # # data_trace_bus=CSV.File("./data/NYISO-data/bus.csv") |> DataFrame
    # baseKV = collect(data_trace_bus[:,Symbol("baseKV")]);
    # kind = collect(data_trace_bus[:,Symbol("kind")]);
    # # data_trace_shunt=CSV.File("../data/NYISO-data/shunt.csv") |> DataFrame
    # type = zeros(1, length(baseKV))
    # baseMVA = zeros(1, length(baseKV))
    # for bus = 1:length(baseKV)
    #     if kind[bus]=="PQ"
    #         type[bus]=1;
    #     elseif kind[bus]=="PV"
    #         type[bus]=2;
    #     elseif kind[bus]=="RE"
    #         type[bus]=3;
    #     end
    # end
    # #	bus_id(1)	type(2)	  Pd(3)	  Qd(4)	  Gs(5)	Bs(6)	area(7)	Vm(8)
    # # Va(9)	baseKV(10)	zone(11)	Vmax(12)	Vmin(13)

    bus = [
    1  3  0.0000  0.0000  0.0000  0.0000  1  1  0  12.5  1  1.1  0.9
    2  1  0.0000  0.0000  0.0000  1.2000  1  1  0  12.5  1  1.1  0.9
    3  1  0.0000  0.0000  0.0000  0.0000  1  1  0  12.5  1  1.1  0.9
    4  1  0.2000  0.1200  0.0000  1.0500  1  1  0  12.5  1  1.1  0.9
    5  1  0.4000  0.2500  0.0000  0.6000  1  1  0  12.5  1  1.1  0.9
    6  1  1.5000  0.9300  0.0000  0.6000  1  1  0  12.5  1  1.1  0.9
    7  1  3.0000  2.2600  0.0000  1.8000  1  1  0  12.5  1  1.1  0.9
    8  1  0.8000  0.5000  0.0000  0.0000  1  1  0  12.5  1  1.1  0.9
    9  1  0.2000  0.1200  0.0000  0.6000  1  1  0  12.5  1  1.1  0.9
   10  1  1.0000  0.6200  0.0000  0.0000  1  1  0  12.5  1  1.1  0.9
   11  1  0.5000  0.3100  0.0000  0.0000  1  1  0  12.5  1  1.1  0.9
   12  1  1.0000  0.6200  0.0000  0.6000  1  1  0  12.5  1  1.1  0.9
   13  1  0.3000  0.1900  0.0000  1.2000  1  1  0  12.5  1  1.1  0.9
   14  1  0.2000  0.1200  0.0000  0.0000  1  1  0  12.5  1  1.1  0.9
   15  1  0.8000  0.5000  0.0000  0.0000  1  1  0  12.5  1  1.1  0.9
   16  1  0.5000  0.3100  0.0000  1.5000  1  1  0  12.5  1  1.1  0.9
   17  1  1.0000  0.6200  0.0000  0.9000  1  1  0  12.5  1  1.1  0.9
   18  1  0.2000  0.1200  0.0000  0.0000  1  1  0  12.5  1  1.1  0.9
   ];

   baseKV = (bus[:,10]);
   Pd_max = (bus[:,3])

   bus_sum_pd = (bus[:,3]);
   frac = bus_sum_pd./sum(bus_sum_pd);

   type = (bus[:,2]);
   # println(frac)
   # println(sum(frac))
   # println(size(frac))
   return bus = (baseKV = (baseKV), Pd_max=(Pd_max), type=(type), frac=(frac))
end

function branch_data()
    # data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
    # id = collect(data_trace_branch[:,Symbol("id")]);
    # fbus = collect(data_trace_branch[:,Symbol("f")]);
    # tbus = collect(data_trace_branch[:,Symbol("t")]);
    # status = collect(data_trace_branch[:,Symbol("status")]);
    # r = collect(data_trace_branch[:,Symbol("R")]);
    # x = collect(data_trace_branch[:,Symbol("X")]);


#	fbus(1)	tbus(2)	r(3)	x(4)	b(5)	rateA(6)	rateB(7)
#   rateC(8)	ratio(9)	angle(10)	status(11)	angmin(12)	angmax(13)
    branch = [
    1   2  0.00004998  0.00035398  0.00000000  999  999  999  0  0  1  -360  360
    2   3  0.00031200  0.00675302  0.00000000  999  999  999  0  0  1  -360  360
    3   4  0.00043098  0.00120403  0.00035000  999  999  999  0  0  1  -360  360
    4   5  0.00060102  0.00167699  0.00049000  999  999  999  0  0  1  -360  360
    5   6  0.00031603  0.00088198  0.00026000  999  999  999  0  0  1  -360  360
    6   7  0.00089600  0.00250202  0.00073000  999  999  999  0  0  1  -360  360
    7   8  0.00029498  0.00082400  0.00024000  999  999  999  0  0  1  -360  360
    8   9  0.00172000  0.00212000  0.00046000  999  999  999  0  0  1  -360  360
    9  10  0.00407002  0.00305299  0.00051000  999  999  999  0  0  1  -360  360
    4  11  0.00170598  0.00220902  0.00043000  999  999  999  0  0  1  -360  360
    3  12  0.00291002  0.00376800  0.00074000  999  999  999  0  0  1  -360  360
   12  13  0.00222202  0.00287699  0.00056000  999  999  999  0  0  1  -360  360
   13  14  0.00480301  0.00621798  0.00122000  999  999  999  0  0  1  -360  360
   13  15  0.00398502  0.00516000  0.00101000  999  999  999  0  0  1  -360  360
   15  16  0.00291002  0.00376800  0.00074000  999  999  999  0  0  1  -360  360
   15  17  0.00372698  0.00459302  0.00100000  999  999  999  0  0  1  -360  360
   17  18  0.00110400  0.00136000  0.00118000  999  999  999  0  0  1  -360  360
    ];


    fbus = (branch[:,1]);
    tbus = (branch[:,2]);
    status = (branch[:,11]);
    r = (branch[:,3]);
    x = (branch[:,4]);

    return branch = (fbus= (fbus), tbus= (tbus), status= (status),
    r=(r), x=(x))
end

function gen_data()
    # data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
## bus(1)	Pg(2)	Qg(3)	Qmax(4)	Qmin(5)	Vg(6)	mBase(7)	status(8)
## Pmax(9)	Pmin(10)	Pc1(11)	Pc2(12)	Qc1min(13)	Qc1max(14)	Qc2min(15)
## Qc2max(16)	ramp_agc(17)	ramp_10(18)	ramp_30(19)	ramp_q(20)	apf(21)
    gen = [
    1  0.0000  0.0000  999  -999  1.0500  100  1   999  0  0  0  0  0  0  0  0  0  0  0  0
    ];


    bus = (gen[:,1]);
    Pmin = (gen[:,10]);
    Pmax = (gen[:,9]);
    Qmin = (gen[:,5]);
    Qmax = (gen[:,4]);

    # id = collect(data_trace_gen[:,Symbol("id")]);
    # bus = collect(data_trace_gen[:,Symbol("bus")]);
    # baseMVA = zeros(length(id),1);
    # Pmin = collect(data_trace_gen[:,Symbol("Pmin")]);
    # Pmax = collect(data_trace_gen[:,Symbol("Pmax")]);
    # Qmin = collect(data_trace_gen[:,Symbol("Qmin")]);
    # Qmax = collect(data_trace_gen[:,Symbol("Qmax")]);
    return gen = (bus= (bus), Pmin= (Pmin), Pmax= (Pmax),Qmin= (Qmin),
    Qmax= (Qmax))
end

function demand_multiplier(bus_struct, pd_raw)
    # p_max = sum(bus_struct.Pd_max)
    # println(p_max)
    # raw_data_max = maximum(sum(pd_raw.pd_rt[:,1:288], dims=1));
    # total_demand_mult = p_max/(raw_data_max);
    # println(total_demand_mult)
    # feeder_mult = bus_struct.Pd_max/p_max*total_demand_mult;
    # println(sum(feeder_mult))
    # println(size(feeder_mult))
    # # println(feeder_mult)
    # println(total_demand_mult*maximum(sum(pd_raw.pd_rt, dims=1)))
    NoBus=length(bus_struct.Pd_max)
    bus_mult = zeros(1,NoBus)
    for bus=1:NoBus
        if bus_struct.Pd_max[bus]==0
            bus_mult[1,bus] = 0;
        else
            bus_mult[1,bus] = bus_struct.Pd_max[bus]/maximum(pd_raw.pd_rt[bus,:]);
        end
    end
    return multiplier = (bus_mult=(bus_mult), total_demand_mult=(sum(bus_mult)))
end
