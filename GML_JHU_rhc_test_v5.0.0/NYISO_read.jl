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

function bus_data(data_trace_bus)
    # data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
    baseKV = collect(data_trace_bus[:,Symbol("baseKV")]);
    kind = collect(data_trace_bus[:,Symbol("kind")]);
    # data_trace_shunt=CSV.File("../data/NYISO-data/shunt.csv") |> DataFrame
    type = zeros(1, length(baseKV))
    baseMVA = zeros(1, length(baseKV))
    for bus = 1:length(baseKV)
        if kind[bus]=="PQ"
            type[bus]=1;
        elseif kind[bus]=="PV"
            type[bus]=2;
        elseif kind[bus]=="RE"
            type[bus]=3;
        end
    end
    bus_sum_pd = [0.0; 569.38; 508.487; 0.0;
    41.242; 50.008; 74.2912; 50.7684; 58.2395;
    0.0; 0.0; 36.4817; 78.1665; 117.133;
    0.0; 0.0; 0.0; 0.0; 67.2924; 0.0; 1.72044;
    24.04; 0.415092; 0.0; 4.25288; 0.0209048;
    52.3195; 48.9072; 18.3068; 9.25575; 1.47497;
    1.69484; 54.9096; 44.8021; 52.2773; 0.0; 16.7049;
    0.0; 40.1738; 33.7294; 13.7582; 40.843; 6.31713;
    0.0; 0.306036; 70.0847; 52.7064; 40.814;
    44.5069; 6.57194; 48.2553; 0.0; 0.0; 18.9136;
    43.4704; 7.10993; 0.0; 0.0835084; 0.0408293;
    63.9603; 132.388; 26.3296; 0.0; 0.0209232;
    0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0;
    0.0; 37.4859; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

    frac = bus_sum_pd./sum(bus_sum_pd);
    # println(frac)
    # println(sum(frac))
    # println(size(frac))
    return bus = (baseKV = (baseKV), type=(type), frac=(frac))
end

function branch_data(data_trace_branch)
    # data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
    id = collect(data_trace_branch[:,Symbol("id")]);
    fbus = collect(data_trace_branch[:,Symbol("f")]);
    tbus = collect(data_trace_branch[:,Symbol("t")]);
    status = collect(data_trace_branch[:,Symbol("status")]);
    r = collect(data_trace_branch[:,Symbol("R")]);
    x = collect(data_trace_branch[:,Symbol("X")]);
    return branch = (id=(id), fbus= (fbus), tbus= (tbus), status= (status),
    r=(r), x=(x))
end

function gen_data(data_trace_gen, bus_struct)
    # data_trace_bus=CSV.File("../data/NYISO-data/bus.csv") |> DataFrame
    id = collect(data_trace_gen[:,Symbol("id")]);
    bus = collect(data_trace_gen[:,Symbol("bus")]);
    baseMVA = zeros(length(id),1);
    Pmin = collect(data_trace_gen[:,Symbol("Pmin")]);
    Pmax = collect(data_trace_gen[:,Symbol("Pmax")]);
    Qmin = collect(data_trace_gen[:,Symbol("Qmin")]);
    Qmax = collect(data_trace_gen[:,Symbol("Qmax")]);
    return gen = (id=(id), bus= (bus), Pmin= (Pmin), Pmax= (Pmax),Qmin= (Qmin),
    Qmax= (Qmax))
end
