using CSV
using Statistics
using DataFrames, DataFramesMeta
# import Lazy # enable @>

# --- read files
curdic = dirname(@__FILE__)
path = string(curdic,"/sidedata")
df = CSV.File(joinpath(path,"df.csv"))|> DataFrame!
df_org = CSV.File(joinpath(path,"df_org.csv"))|> DataFrame!
df_err = CSV.File(joinpath(path,"df_err.csv"))|> DataFrame!
df_crt = CSV.File(joinpath(path,"df_crt.csv"))|> DataFrame!

#==============

===============#
global alpha, beta, s, c ,Old_crit ,New_crit ,
    t0 ,kao , F, L, c2

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

#=====
LTM global
====#

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
            elseif walk+"i"+item in ["oldinew","newiold"]

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
