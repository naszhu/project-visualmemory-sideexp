
df_all = @> begin
    df_org
    @where(:FileCondi .!= "CMat")
    @transform(
    Probtype = map(:Stimkind) do t
        t == 1 ? "CM" :
        t == 0 ? "AN" :
        t == 3 ? "VM" :
        missing
    end,
    Oldnew = map(:Old) do t
        t == 1 ? "old" :
        t == 2 ? "new" :
        missing
    end,
    Error = 1 .- :Correctness
    )
    end
g = DataFrames.groupby(df_all, [:Oldnew,:Setsize,:Probtype,:Lag,:Error,:FileCondi,:RT])
df = aggregate(g,sum)
aggregate(df_all,[:Oldnew,:Setsize,:Probtype,:Lag,:Error,:FileCondi],RT=mean(:RT))
