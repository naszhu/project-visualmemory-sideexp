using CSV
using Statistics
using DataFrames, DataFramesMeta
using Lazy # enable @>
df = DataFrame(a = repeat([1, 2, 3, 4], outer=[2]),
                             b = repeat([2, 1], outer=[4]),
                             c = 1:8);

 new = by(DataFrames.groupby(df, :a), [sum, x->mean(skipmissing(x))])
 new

 by(df, :a, :c => sum)

 by(df, [:a,:b], :c=>sum)

module now

a = 10
function change()
    global a
    a=101
end
change()
println(a)


end
 # module now
