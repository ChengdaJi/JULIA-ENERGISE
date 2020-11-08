function NYISO_shunt_data(data_trace_bus,data_trace_shunt)
    # number of bus in the system
    NoShunt = length(collect(data_trace_shunt[:,Symbol("id")]));
    NoBus = length(collect(data_trace_bus[:,Symbol("id")]));
    # Real demand at that instance
    P = reshape(collect(data_trace_shunt[:,Symbol("P")]), NoShunt, 1);
    # Reactive demand at that instance
    Q = reshape(collect(data_trace_shunt[:,Symbol("Q")]), NoShunt, 1);
    find_bus = reshape(collect(data_trace_shunt[:,Symbol("bus")]), NoShunt, 1);
    Q_max = reshape(collect(data_trace_shunt[:,Symbol("Qmax")]), NoShunt, 1);
    kind = reshape(collect(data_trace_shunt[:,Symbol("kind")]), NoShunt, 1);
    Vsp = reshape(collect(data_trace_shunt[:,Symbol("Vsp")]), NoShunt, 1);
    type = zeros(NoShunt,1);
    P_agg = sum(P);
    # Q_agg = sum(Q);
    frac_P = P/P_agg;
    # frac_Q = Q/Q_agg;
    # println(size(frac_P))
    bus_switch_shunt = zeros(NoBus, 1)
    for shunt = 1:NoShunt
        if kind[shunt]=="LOAD"
            type[shunt]=1;
        elseif kind[shunt]=="CTRL_V"
            type[shunt]=2;
            # println(shunt)
            # println(find_bus[shunt,1])
            bus_switch_shunt[find_bus[shunt,1],1]=shunt;
            # println(bus_switch_shunt[find_bus[shunt,1],1])
        else
            type[shunt]=3;
        end
    end

    return shunt = (Q_max=(Q_max), bus_switch_shunt=(bus_switch_shunt),
    P=(P), Q=(Q), type=(type), Vsp=(Vsp), frac_P=(frac_P), find_bus = (find_bus))
end

# function bus_data(system)
#     # # data_trace_bus=CSV.File("./data/NYISO-data/bus.csv") |> DataFrame
#     # baseKV = collect(data_trace_bus[:,Symbol("baseKV")]);
#     # kind = collect(data_trace_bus[:,Symbol("kind")]);
#     # # data_trace_shunt=CSV.File("../data/NYISO-data/shunt.csv") |> DataFrame
#     # type = zeros(1, length(baseKV))
#     # baseMVA = zeros(1, length(baseKV))
#     # for bus = 1:length(baseKV)
#     #     if kind[bus]=="PQ"
#     #         type[bus]=1;
#     #     elseif kind[bus]=="PV"
#     #         type[bus]=2;
#     #     elseif kind[bus]=="RE"
#     #         type[bus]=3;
#     #     end
#     # end
#     # #	bus_id(1)	type(2)	  Pd(3)	  Qd(4)	  Gs(5)	Bs(6)	area(7)	Vm(8)
#     # # Va(9)	baseKV(10)	zone(11)	Vmax(12)	Vmin(13)
#
#
#
#    baseKV = (system["bus"][:,10]);
#    Pd_max = (system["bus"][:,3])
#    # println(sum(Pd_max))
#
#    bus_sum_pd = (system["bus"][:,3]);
#    frac = bus_sum_pd./sum(bus_sum_pd);
#
#    type = (system["bus"][:,2]);
#    # println(frac)
#    # println(sum(frac))
#    # println(size(frac))
#    return bus = (baseKV = (baseKV), Pd_max=(Pd_max), type=(type), frac=(frac))
# end

function NYISO_bus_data(data_trace_bus)
    # data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
    NoBus=length(collect(data_trace_bus[:,Symbol("id")]));
    baseKV = reshape(collect(data_trace_bus[:,Symbol("baseKV")]),NoBus,1);
    kind = reshape(collect(data_trace_bus[:,Symbol("kind")]),NoBus,1);
    # data_trace_shunt=CSV.File("../data/NYISO-data/shunt.csv") |> DataFrame
    type = zeros(length(baseKV),1)
    for bus = 1:length(baseKV)
        if kind[bus]=="PQ"
            type[bus]=1;
        elseif kind[bus]=="PV"
            type[bus]=2;
        elseif kind[bus]=="RE"
            type[bus]=3;
        end
    end
    id=reshape(1:length(type),length(type),1)
   return bus = (id=(id),baseKV = (baseKV), type=(type))
end


# function branch_data(system)
#
#     fbus = (system["branch"][:,1]);
#     tbus = (system["branch"][:,2]);
#     status = (system["branch"][:,11]);
#     r = (system["branch"][:,3]);
#     x = (system["branch"][:,4]);
#
#     return branch = (fbus= (fbus), tbus= (tbus), status= (status),
#     r=(r), x=(x))
# end

function NYISO_branch_data(data_trace_branch)
    id = collect(data_trace_branch[:,Symbol("id")]);
    fbus = collect(data_trace_branch[:,Symbol("f")]);
    tbus = collect(data_trace_branch[:,Symbol("t")]);
    status = collect(data_trace_branch[:,Symbol("status")]);
    r = collect(data_trace_branch[:,Symbol("R")]);
    x = collect(data_trace_branch[:,Symbol("X")]);
    NoBr=length(id);
    id=reshape(id, NoBr, 1)
    fbus=reshape(fbus, NoBr, 1)
    tbus=reshape(tbus, NoBr, 1)
    status = reshape(status, NoBr, 1)
    r=reshape(r, NoBr, 1)
    x=reshape(x, NoBr, 1)
    return branch = (fbus= (fbus), tbus= (tbus), status= (status),
    r=(r), x=(x))
end

# function gen_data(system)
#     # data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
# ## bus(1)	Pg(2)	Qg(3)	Qmax(4)	Qmin(5)	Vg(6)	mBase(7)	status(8)
# ## Pmax(9)	Pmin(10)	Pc1(11)	Pc2(12)	Qc1min(13)	Qc1max(14)	Qc2min(15)
# ## Qc2max(16)	ramp_agc(17)	ramp_10(18)	ramp_30(19)	ramp_q(20)	apf(21)
#     # gen = [
#     # 1  0.0000  0.0000  999  -999  1.0500  100  1   999  0  0  0  0  0  0  0  0  0  0  0  0
#     # ];
#
#
#     bus = (system["gen"][:,1]);
#     Pmin = (system["gen"][:,10]);
#     Pmax = (system["gen"][:,9]);
#     # println(sum(Pmax))
#
#     Qmax = (system["gen"][:,4]);
#     # println(sum(Qmax))
#
#     Qmin = (system["gen"][:,5]);
#     # Qmin=-Qmax
#     # println(sum(Qmin))
#     # id = collect(data_trace_gen[:,Symbol("id")]);
#     # bus = collect(data_trace_gen[:,Symbol("bus")]);
#     # baseMVA = zeros(length(id),1);
#     # Pmin = collect(data_trace_gen[:,Symbol("Pmin")]);
#     # Pmax = collect(data_trace_gen[:,Symbol("Pmax")]);
#     # Qmin = collect(data_trace_gen[:,Symbol("Qmin")]);
#     # Qmax = collect(data_trace_gen[:,Symbol("Qmax")]);
#     return gen = (bus= (bus), Pmin= (Pmin), Pmax= (Pmax),Qmin= (Qmin),
#     Qmax= (Qmax))
# end

function NYISO_gen_data(data_trace_gen)
    # println(sum(Qmin))
    NoGen = length(collect(data_trace_gen[:,Symbol("id")]));
    bus = reshape(collect(data_trace_gen[:,Symbol("bus")]), NoGen,1);
    Pmin = reshape(collect(data_trace_gen[:,Symbol("Pmin")]), NoGen,1);
    Pmax = reshape(collect(data_trace_gen[:,Symbol("Pmax")]), NoGen,1);
    Qmin = reshape(collect(data_trace_gen[:,Symbol("Qmin")]), NoGen,1);
    Qmax = reshape(collect(data_trace_gen[:,Symbol("Qmax")]), NoGen,1);
    return gen = (bus= (bus), Pmin= (Pmin), Pmax= (Pmax), Qmin= (Qmin),
    Qmax= (Qmax))
end
