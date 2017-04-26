# Includes etc. Imports and blah

sec = [1,2]

data = rand

gem = model(IpoptSolver())

@variables gem begin
    pd[sec] >= 0,


end

@constraints gem beging
    C[sec] .== muH[sec] + alphaHLES[sec]/(PD[sec]) * (Y - sum(PD[s] * muH[s] for s in sec)),

end
