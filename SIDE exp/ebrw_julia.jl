using CSV
using Statistics
using DataFrames, DataFramesMeta
import Dates
# import Lazy # enable @>

# --- read files
curdic = "/home/shulai/Google Drive/IUB/Shiffrin Lab/mix_experiment/SIDE exp"#dirname(@__FILE__)

path = string(curdic,"/sidedata")
df = CSV.File(joinpath(path,"df.csv"))|> DataFrame!
df_org = CSV.File(joinpath(path,"df_org.csv"))|> DataFrame!
df_err = CSV.File(joinpath(path,"df_err.csv"))|> DataFrame!
df_crt = CSV.File(joinpath(path,"df_crt.csv"))|> DataFrame!;


# global alpha, beta, s, c ,Old_crit ,New_crit ,
#     t0 ,kao , F, L, c2

alpha=Dict(); beta =Dict(); s = Dict(); F = Dict(); L = Dict(); c=Dict(); t0=Dict(); Old_crit=Dict(); New_crit=Dict()


# boost = 1.05
alpha["all"] = 0.19
# alpha['CM'] = 0.1969
# alpha['VM'] = 0.1969
# alpha["AN"]=0.84

beta["all"] =1.2
# beta['CM'] = 1.228
# beta['VM'] = 1.228
# beta["AN"]= 1.708

s["all"]=0.01
s["ss24"]=0.01
s["ss8"]=0.01
# s["AN"] = 0.01
# s["CM"] = 0.02
# s["VM"] = 0.02
c["ANpure"] = 0.3938
c["CMpure"] = 0.3938
c["VMpure"] = 0.3938
c["MIX"] = 0.3938
c["MIX2"] = 0.3938
# c["CMat"] = 0.3938
c2=0.3
Old_crit["ANpure"] = 1.9197
New_crit["ANpure"] = -2.33
Old_crit["CMpure"] = 1.9197
New_crit["CMpure"] = -2.33
Old_crit["VMpure"] = 1.9197
New_crit["VMpure"] = -2.33
Old_crit["MIX"] = 1.9197
New_crit["MIX"] = -2.33
Old_crit["MIX2"] = 1.9197
New_crit["MIX2"] = -2.33

t0["all"] = 699.98
t0["ann"] =  699.98
kao = 37

#---CM
F["ANpure_AN_oldiold_oldinew"] = 0.2

F["CMpure_CM_oldiold_oldinew"] = 0.2
L["CMpure_CM_oldiold_newinew"] = 0.2
L["CMpure_CM_oldinew_newiold"] = 0

F["VMpure_VM_oldiold_oldinew"] = 0.2
L["VMpure"] = 0.2

F["MIX_CM_oldiold_oldinew"] = 0.2
L["MIX_CM_oldiold_newinew"] = 0.2
L["MIX_CM_oldinew_newiold"] = 0.2
F["MIX_AN_oldiold_oldinew"] = 0.2

F["MIX2_CM_oldiold_oldinew"] = 0.2
L["MIX2_CM_oldiold_newinew"] = 0.2
L["MIX2_CM_oldinew_newiold"] = 0.2
F["MIX2_AN_oldiold_oldinew"] = 0.2

# global vary_ss, give_c2
vary_ss=0
give_c2=0

# global nameset,is_search_allcondi_besideCMat,search_MIX1
is_search_allcondi_besideCMat=0
search_MIX1 = 1
search_MIX2 = 1
nameset = CategoricalArray(df_org.FileCondi) |> levels;

#=====
LTM global
======#

function assign_LTM_global(Filecondi,item_condi, walk, item)

    global F,L,Fnow,Lnow
    if Filecondi != "VMpure"
        if item_condi == "CM"

            if walk * "i" * item in ["oldiold","oldinew"]

                Fnow = F[Filecondi*"_CM_oldiold_oldinew"]
            else Fnow = 0
            end

            if walk * "i" * item in ["oldiold","newinew"]

                Lnow = L[Filecondi*"_CM_oldiold_newinew"]
            elseif walk * "i" * item in ["oldinew","newiold"]

                Lnow = L[Filecondi*"_CM_oldinew_newiold"]
            else Lnow=0
            end

        elseif item_condi == "AN"

            if walk * "i" * item in ["oldiold","oldinew"]

                Fnow = F[Filecondi*"_AN_oldiold_oldinew"]
            else Fnow=0
            end

            Lnow = 0

        else
            println("wrong condi 1",Filecondi,item_condi, walk*"i"*item)
        end
    elseif Filecondi == "VMpure"

        if item_condi == "VM"
            if walk*"i"*item in ["oldiold","oldinew"]
                Fnow = F["VMpure_VM_oldiold_oldinew"]
            else
                Fnow = 0
            end
        else print("wrong condi 2")
        end

        Lnow = L["VMpure"]
    else
        print("wrong filecondi")
    end
    return(Fnow + Lnow)
end
assign_LTM_global("MIX","AN","old","old")


function calcA(df)

    dnow = df
    # a = Array{Union{Nothing,Float64}}(nothing,size(dnow)[1], 8)  #activation
    a =   zeros(Float64, size(dnow)[1], 8)
    betanow=beta["all"]
    alphanow=alpha["all"]

    m = [(j^(-betanow) + alphanow) for j in 1:8]
    for j in 1:8

        if vary_ss==0
            a[(dnow.Lag .== j )|> Array,j] .= m[j]
            a[(dnow.Lag .!= j )|> Array,j] .= m[j] .* s["all"]
        else
            a[(dnow.Lag .!= j .& dnow.Setsize .!= 8)  |> Array,j] .= m[j] .* s["ss24"]
            a[(dnow.Lag .!= j .& dnow.Setsize .== 8)  |> Array,j] .= m[j] .* s["ss8"]
            a[(dnow.Lag .== j) |> Array,j] .= m[j]

        end
    end

    for i in 1:size(a)[1]

        a[i,df.Setsize[i]+1:end] .= 0 #a_ij suit for the correct amount of setsize
    end

    debug = 0
    if debug==1
        for i in range(a.shape[0])
            println(i, "beginn","a is",a[i].round(5),
                  "n m is", np.array(m).round(3),"n Probtype is",
                  "n Lag is ",dnow.Lag[i],
                  "n Setsize is", dnow.Setsize[i],
                  "n Probtype:",
                  dnow.Probtype[i],
                  "n Oldnew: ",dnow.Oldnew[i],
                  "n Ai is", a[i,:].sum().round(2),
                  "n snow", s["all"],
                 "n------------------------------------------------" )
        end
    end

    A = sum(a, dims = 2)

    return(A)
end

A0=calcA(df)
present_now = [df.FileCondi[i] *
        "-- " * (A0[i] |> string) *
      "--" *  (df.Setsize[i] |> string)*
      "--" * df.Probtype[i] *
      "--" * df.Oldnew[i] for i in 1:size(A0)[1]];
display(sort(present_now))


# F["AN_oldiold_oldinew"] = 0.099
function calcp(df)

    A = calcA(df)
    dnow = df
    p = Iterators.repeated(3.0, size(A)[1]) |> collect
    filecondis = CategoricalArray(dnow.FileCondi) |> levels; #get categories name

    for icondi in filecondis

        #only ieterate through probs in that condition
        probs = filter( t -> t[:FileCondi] == icondi, dnow).Probtype |> CategoricalArray |> levels
        for iprob in probs

            for ion in ["old","new"]

                tf_ion = (dnow.Oldnew .== ion)
                tf_iprob = (dnow.Probtype .== iprob)
                tf_icondi = (dnow.FileCondi .== icondi)

                tf_combi1 = (tf_ion .& tf_iprob) |> Array
                tf_all = (tf_combi1 .& tf_icondi) |> Array

                IR_old_current = assign_LTM_global(icondi,iprob,"old",ion)
                IR_new_current = assign_LTM_global(icondi, iprob,"new",ion)

                p[tf_all] = (A[tf_all] .+ IR_old_current)./
                (A[tf_all] .+ IR_old_current .+ c[icondi] .+ IR_new_current)
            end
        end
    end

    return(p)

end

p=calcp(df)
present_now = [df.FileCondi[i] *
        "-- " * (p[i] |> string) *
      "--" *  (df.Setsize[i] |> string)*
      "--" * df.Probtype[i] *
      "--" * df.Oldnew[i] for i in 1:size(p)[1]];
display(sort(present_now))



function calc_theoretical_RW(df, which)

    which_return = "crt"
    A = calcA(df)
    p = calcp(df)
    q = 1 .- p

    dnow = df
    filecondis = CategoricalArray(dnow.FileCondi) |> levels; #get categories name

    p_resp_old = zeros(size(p)[1])
    pred_correct = zeros(size(p)[1])
    pred_rt = zeros(size(p)[1])

    theta1 = zeros(size(p)[1])
    theta2 = zeros(size(q)[1])
    theta11 = zeros(size(p)[1])
    theta22 = zeros(size(q)[1])
    exp_nstep = zeros(size(p)[1])


    for ifile in filecondis

        gen_tf = (dnow.FileCondi .== ifile)|> Array
        AA = Old_crit[ifile]
        BB = -New_crit[ifile] # bb is a postive number

        if which_return=="crt"

            pq = (p[gen_tf]./q[gen_tf])
            theta1[gen_tf] = (pq .^(AA+BB) .+1)./(pq .^(AA+BB) .-1)
            theta2[gen_tf] = (pq .^BB .+1)./(pq .^BB .-1)

            tf = ((p .!= q) .& ((dnow.Oldnew .== "old").& gen_tf )) |> Array
            exp_nstep[tf] = (1 ./ (p[tf] .-q[tf])).*(theta1[tf].*(AA+BB) .- theta2[tf].*BB)

            tf = ((p .== q) .& ((dnow.Oldnew .== "old").& gen_tf)) |> Array
            exp_nstep[tf] .= (AA/3)*(2*BB+AA)

            theta11[gen_tf] = (pq .^ (-(AA+BB)) .+1) ./ (pq.^(-(AA+BB)) .-1)
            theta22[gen_tf] = (pq .^ -AA .+1)./(pq .^ -AA .-1)

            tf = ((p .!= q) .& (dnow.Oldnew .== "new") .& gen_tf) |> Array
            exp_nstep[tf] = (1 ./(q[tf] .- p[tf])).*(theta11[tf].*(AA+BB) .- theta22[tf].*AA)

            tf = ((p .== q) .& (dnow.Oldnew .== "new") .& gen_tf) |> Array
            exp_nstep[tf] .= (BB/3)*(2*AA+BB);

        #------------------- correct response
        elseif which_return == "err"

            qp = (q[gen_tf]./p[gen_tf])
            qptfn = ((p .!= q) .& gen_tf) |> Array
            qptf = ((p .== q) .& gen_tf) |> Array

            p_resp_old[qptfn] = ((1 .- qp.^BB) ./ (1 .-qp.^(AA+BB)))[Array(p[gen_tf] .!= q[gen_tf])];
            p_resp_old[qptf] .= BB/(AA +BB);
        end
    end
    #------------------


    if which_return=="crt"
        suprise = (((dnow.FileCondi.=="MIX") .| (dnow.FileCondi.=="MIX2")) .&
                (dnow.Probtype.=="AN") .& (dnow.Oldnew.=="new") )|> Array

        notsuprise = .!suprise

        pred_rt[suprise] = t0["ann"] .+ kao .* exp_nstep[suprise]
        pred_rt[notsuprise] = t0["all"] .+ kao .* exp_nstep[notsuprise]
    end

    #------------------

    if which_return == "err"
        pred_correct[Array(dnow.Oldnew.=="old")] = p_resp_old[Array(dnow.Oldnew.=="old")]
        pred_correct[Array(dnow.Oldnew.=="new")] = 1-p_resp_old[Array(dnow.Oldnew.=="new")]

    end


    if which == "crt"
        return(pred_rt)
    elseif which == "err"
        return(pred_correct)
    end

end

rw = calc_theoretical_RW(df,"crt")
present_now = [df.FileCondi[i] *
        "-- " * (rw[i] |> string) *
      "--" *  (df.Setsize[i] |> string)*
      "--" * df.Probtype[i] *
      "--" * df.Oldnew[i] for i in 1:size(rw)[1]];
display(sort(present_now))


w=Dict()

w["new_rt"] = 4*2
w["old_rt"] = 1*2
w["new_err"] = 4
w["old_err"] = 1

function calc_wssd(df_err,df_crt)


    df_sub_err = df_err |> copy
    df_sub_crt = df_crt |> copy

    if search_MIX1 == 1
        df_sub_err = filter(t -> t[:FileCondi]=="MIX", df_sub_err)
        df_sub_crt = filter(t -> t[:FileCondi]=="MIX", df_sub_crt)
    end
    pred_correct = calc_theoretical_RW(df_sub_err,"err")
    pred_crt = calc_theoretical_RW(df_sub_crt,"crt")

    df_sub_crt.pred_crt = pred_crt ./ 1000 #translate crt to seconds.
    df_sub_err.pred_error = 1 .- pred_correct

    df_sub_crt.RT = df_sub_crt.RT ./1000

    df_sub_crt.SSD_RT = ((df_sub_crt.RT .- df_sub_crt.pred_crt) .^2)
    df_sub_err.SSD_err = ((df_sub_err.Error .- df_sub_err.pred_error).^2)

    df_sub_crt.wSSD_RT = [df_sub_crt.SSD_RT[i] * w[df_sub_crt.Oldnew[i] * "_rt"]
        for i in 1:size(df_sub_crt)[1]]

    df_sub_err.wSSD_err = [df_sub_err.SSD_err[i] * w[df_sub_err.Oldnew[i] * "_err"]
        for i in 1:size(df_sub_err)[1]]

    WSSD = (df_sub_err.wSSD_err |> sum)  .+ (df_sub_crt.wSSD_RT |> sum)

    return(WSSD)
end

time1=Dates.now()
calc_wssd(df_err,df_crt)
println(Dates.now()-time1)
