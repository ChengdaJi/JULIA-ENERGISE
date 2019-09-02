A=[1 -1 2; -2 3 4]
sqrt.(A)
# include("GML_RHC.jl")
# A_new = positive_array(A)
# println(A_new)

# function positive_array(Array_)
#     for row=1:length(Array_[:,1])
#         for column=:length(Array_[1,:])
#             if Array_[row, column]<=0
#                 Array_[row, column]=0;
#             end
#         end
#     end
#     return Array_
# end
