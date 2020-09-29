# using CSV
# using DataFrames
#
# read_old_data(2)

function read_old_data(t, multiplier)
    filename1 = "P_rsrv_feedback.csv"
    data_trace_P_rsrv = CSV.File(filename1) |> DataFrame
    P_rsrv_feedback = data_trace_P_rsrv[1:t-1, 1];
    println(P_rsrv_feedback)
    # println("here")
    # filename2 = string("M15P_rate1Time",t,".csv")
    filename2 = string("result/M1P_rate0.5Time",t,".csv")
    data_trace_all = CSV.File(filename2) |> DataFrame
    B_sum = collect(data_trace_all[:,Symbol("B")])
    # println(B_sum)
    B_feedback_one = vcat(zeros(3,1),
        (B_sum[1]/(12*multiplier))*ones(12*multiplier,1))
    # B_feedback_one = ones(12*multiplier,1)
    B_feedback = vcat(0, repeat(B_feedback_one, multiplier))
    # println(B_feedback)
    return old = (B_feedback=(B_feedback), P_rsrv_feedback=(P_rsrv_feedback))
end

    # filename2 = string("M15P_rate1Time",t,".csv")
    # data_trace_all = CSV.File(filename2) |> DataFrame
