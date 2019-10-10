function pg_traj_raw(t, pg_raw_one, pd_traj, p_rate, BN,T, Pred_length);
    pg_traj1_raw=hcat(pg_raw_one["Pg1"]);
    pd_current_sum_1=sum(pd_traj[1:8,1:T]);
    pd_portion1=zeros(8,1)
    for feeder=1:8
        pd_portion1[feeder] = sum(pd_traj[feeder,1:T])/pd_current_sum_1;
    end
    println(sum(pd_portion1))
    mu1=pd_portion1*reshape(pg_traj1_raw[t:t+T-1],1,T)/sum(pg_raw_one["Pg1"])*pd_current_sum_1*p_rate;
    println(sum(mu1)/pd_current_sum_1)

    pg_traj2_raw=hcat(pg_raw_one["Pg2"], pg_raw_one["Pg2"]);
    pd_current_sum_2=sum(pd_traj[9:12,1:T]);
    pd_portion2=zeros(4,1)
    for feeder=9:12
        pd_portion2[feeder-8] = sum(pd_traj[feeder,1:288])/pd_current_sum_2;
    end
    println(sum(pd_portion2))

    mu2=pd_portion2*reshape(pg_traj2_raw[t:t+T-1],1,T)/sum(pg_raw_one["Pg2"])*pd_current_sum_2*p_rate;
    println(sum(mu2)/pd_current_sum_2)
    mu=vcat(mu1,mu2);
    println(sum(mu)/(sum(pd_traj[1:12,1:T])))
    pg=pg_struct(mu,0,0,0);
    return pg
end

function solar_output(pg)
    # write the solar file
    println("writing solar file");
    time_stamp=zeros(576,1)
    for i=1:576
        time_stamp[i]=i
    end
    println(size(pg))
    Solar=hcat(time_stamp, vcat(transpose(pg), 0.8*transpose(pg[:, 1:144]), 1.2*transpose(pg[:, 145:288])));
    println(size(Solar))
    CSV.write("Solar.csv", DataFrame(Solar, [:timestamp, :bnk139_Pf_1,
    :bnk139_Pf_2, :bnk139_Pf_3, :bnk139_Pf_4, :bnk239_Pf_5, :bnk239_Pf_6,
    :bnk239_Pf_7, :bnk239_Pf_8, :bnk276_Pf_9, :bnk276_Pf_10, :bnk276_Pf_11,
    :bnk276_Pf_12]));
    # Reactive_Power=[time_stamp; opt_value.Qf]';
    # CSV.write("Reactive_Power.csv", DataFrame(Power, [:timestamp, :bnk139_Pf_1,
    # :bnk139_Pf_2, :bnk139_Pf_3, :bnk139_Pf_4, :bnk239_Pf_5, :bnk239_Pf_6,
    # :bnk239_Pf_7, :bnk239_Pf_8, :bnk276_Pf_9, :bnk276_Pf_10, :bnk276_Pf_11,
    # :bnk276_Pf_12]));
    # V1=(hcat(opt_value.v[1,:],opt_value.v[1,:],opt_value.v[1,:],opt_value.v[1,:]))
    # V2=(hcat(opt_value.v[2,:],opt_value.v[2,:],opt_value.v[2,:],opt_value.v[2,:]))
    # V3=(hcat(opt_value.v[3,:],opt_value.v[3,:],opt_value.v[3,:],opt_value.v[3,:]))
    # V_sqrt=sqrt.(hcat(V1,V2,V3));
    # # println(size(v_final))
    # Voltage=hcat(time_stamp',V_sqrt);
    # CSV.write("Voltage.csv", DataFrame(Voltage, [:timestamp, :bnk139_Pf_1,
    # :bnk139_Pf_2, :bnk139_Pf_3, :bnk139_Pf_4, :bnk239_Pf_5, :bnk239_Pf_6,
    # :bnk239_Pf_7, :bnk239_Pf_8, :bnk276_Pf_9, :bnk276_Pf_10, :bnk276_Pf_11,
    # :bnk276_Pf_12]));
    # println("GML - Finish writting files! ")
end
