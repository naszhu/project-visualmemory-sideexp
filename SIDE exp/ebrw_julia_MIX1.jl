using CSV
using Statistics
using DataFrames, DataFramesMeta
import Dates
using NLopt
using JuMP
import Distributions: Uniform
include("core_grant_func.jl")
# import Lazy # enable @>

# --- read files
curdic = "/home/shulai/Google Drive/IUB/Shiffrin Lab/mix_experiment/SIDE exp"#dirname(@__FILE__)

path = string(curdic,"/sidedata")
df = CSV.File(joinpath(path,"df.csv"))|> DataFrame!
df_org = CSV.File(joinpath(path,"df_org.csv"))|> DataFrame!
df_err = CSV.File(joinpath(path,"df_err.csv"))|> DataFrame!
df_crt = CSV.File(joinpath(path,"df_crt.csv"))|> DataFrame!;

df_err = filter(t -> t[:FileCondi]=="MIX", df_err)
df_crt = filter(t -> t[:FileCondi]=="MIX", df_crt)

alpha=Dict(); beta =Dict(); s = Dict(); F = Dict(); L = Dict(); c=Dict(); t0=Dict(); Old_crit=Dict(); New_crit=Dict()

alpha["all"] = 0.19
beta["all"] =1.2
s["all"] =0.01
c["MIX"] = 0.3938
Old_crit["MIX"] = 1.9197
New_crit["MIX"] = -2.33

t0["all"] = 699.98
t0["ann"] =  699.98
kao = 37

#---CM
F["MIX_CM_oldiold_oldinew"] = 0.2
L["MIX_CM_oldiold_newinew"] = 0.2
L["MIX_CM_oldinew_newiold"] = 0.2
F["MIX_AN_oldiold_oldinew"] = 0.2


all_vars = [alpha, beta, s, c ,Old_crit ,New_crit , t0 ,kao , F, L]

# global vary_ss, give_c2
vary_ss=0
give_c2=0

# global nameset,is_search_allcondi_besideCMat,search_MIX1
is_search_allcondi_besideCMat=1
search_MIX1 = 0
search_MIX2 = 0
nameset = CategoricalArray(df_org.FileCondi) |> levels;

#===temp===#
wssdn = calc_wssd(df_err,df_crt, all_vars)
println(wssdn)
#===temp===#


bdd = (
    (0.1, 2.9),  # alpha,
    (0.01, 4),  # beta,
    (0.01, 0.999),  # sall

    (0.01, 3),  # c_MIX,
    (1, 10),  # BDDo_MIX,
    (-10, 1),  # BDDn_MIX,
    (100, 900),  # t0,
    (100, 900),  # t0an,
    (10, 100),  # kappa,
    (0.0001, 4),  # F_MIX_CM_oldiold_oldinew,
    (0.0001, 4),  # L_MIX_CM_oldiold_newinew,
    (0.0001, 4),  # L_MIX_CM_oldinew_newiold,
    (0.0001, 4),  # F_MIX_AN_oldiold_oldinew,
      )
bdd_lower = [i[1] for i in bdd];
bdd_upper = [i[2] for i in bdd];

function random_start0()
    global alpha, beta, s, c ,Old_crit ,New_crit , t0 ,kao , F, L

    alpha=Dict(); beta =Dict(); s = Dict(); F = Dict(); L = Dict()

    alpha["all"] = rand(Uniform(0.1,3))
    beta["all"] = rand(Uniform(0.01,4))
    s["all"] = rand(Uniform(0.01,0.99))

    c["MIX"] = rand(Uniform(0.01,0.99))
    Old_crit["MIX"] = rand(Uniform(1,10))
    New_crit["MIX"] = rand(Uniform(-10,1))

    t0["all"] =  rand(Uniform(1,900))
    t0["ann"] =  rand(Uniform(1,900))
    kao = rand(Uniform(1,100))

    F["MIX_CM_oldiold_oldinew"] = rand(Uniform(0.0001,1))
    L["MIX_CM_oldiold_newinew"] = rand(Uniform(0.0001,1))
    L["MIX_CM_oldinew_newiold"] = rand(Uniform(0.0001,1))
    F["MIX_AN_oldiold_oldinew"] = rand(Uniform(0.0001,1))


    param_dic=Array([alpha["all"], beta["all"], s["all"],
            c["MIX"], Old_crit["MIX"], New_crit["MIX"],
            t0["all"],t0["ann"], kao,
            F["MIX_CM_oldiold_oldinew"], L["MIX_CM_oldiold_newinew"],
            L["MIX_CM_oldinew_newiold"],F["MIX_AN_oldiold_oldinew"]])

    return(param_dic)
end

param_dic = random_start0()
print(param_dic)

random_start() = [rand(Uniform(bdd_lower[i],bdd_upper[i])) for i in 1:length(random_start0())]
random_start();




function nlopt_opt(x::Vector, grad::Vector)

    alpha=Dict(); beta =Dict(); s = Dict(); F = Dict(); L = Dict(); c=Dict(); t0=Dict(); Old_crit=Dict(); New_crit=Dict()

    alpha["all"], beta["all"], s["all"],
    c["MIX"], Old_crit["MIX"], New_crit["MIX"],
    t0["all"],t0["ann"], kao,
    F["MIX_CM_oldiold_oldinew"], L["MIX_CM_oldiold_newinew"],
    L["MIX_CM_oldinew_newiold"],F["MIX_AN_oldiold_oldinew"] = x

    all_vars = [alpha, beta, s, c ,Old_crit ,New_crit , t0 ,kao , F, L ]

    return(calc_wssd(df_err,df_crt, all_vars))
end

function myconstraint1(x::Vector, grad::Vector)
    return(x[13]-x[10])
end

function myconstraint2(x::Vector, grad::Vector)
    return(x[12]-x[11])
end


function optim_all()
    println("\n ************start optimization now****************")
    opt = Opt(:GN_ISRES, length(random_start()))#Opt(:LD_SLSQP, 3)
    # opt = Opt(:LD_MMA, length(random_start()))
    # opt = Opt(:LN_COBYLA, length(random_start()))
    # opt = Opt(:LD_SLSQP, length(random_start()))
    # opt = Opt(:LN_SBPLX, length(random_start()))

    opt.lower_bounds = bdd_lower
    opt.upper_bounds = bdd_upper

    opt.xtol_rel = 1e-4
    opt.maxtime = 60*3
    # opt.maxeval = 10000
    opt.min_objective = nlopt_opt

    inequality_constraint!(opt,(x,g) -> myconstraint1(x,g), 1e-8)
    inequality_constraint!(opt,(x,g) -> myconstraint2(x,g), 1e-8)

    time1=Dates.now()
    (minf,minx,ret) = optimize(opt, random_start())
    println(Dates.now()-time1)


    numevals = opt.numevals
    println("got $minf at $minx after $numevals iterations (returned $ret)")

    return(minf, minx)
end
# optim_all()

num_iter = 5
results=zeros(num_iter)
resultsx=zeros(num_iter,length(bdd))
for i in 1:num_iter
    print("**",i)
    (results[i],resultsx[i,:]) = optim_all()
end
print("\n******* results **********\n",
findmin(results)[1],"**",resultsx[findmin(results)[2],:])


all_results = (cat(dims=2, resultsx, reshape(results,length(results),1))) |> DataFrame
names!(all_results,[:alpha_all, :beta_all, :s_ss2,
:c_MIX, :Old_crit_MIX, :New_crit_MIX,
:t0_all, :t0_ann, :kao,
:F_MIX_CM_oldiold_oldinew, :L_MIX_CM_oldiold_newinew,
:L_MIX_CM_oldinew_newiold, :F_MIX_AN_oldiold_oldinew, :wssd])

filename = "ebrw_julia_MIX1_Nss_" * string(time) * ".csv"
CSV.write(filename,all_results)
