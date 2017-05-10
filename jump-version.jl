using Ipopt, JuMP, DataFrames

M = Model(solver = IpoptSolver(print_level = 0))

data = readtable("data.csv")

### PARAMETERS ###

## Production block

KZ = data[:KZ]              # inital capital demand
KSZ = sum(KZ)               # capital endowment
LZ  = data[:LZ]             # initial labor demand
LSZ = sum(LZ)               # labor endowment
PKZ = 1.                    # initial capital price
PLZ = 1.                    # initial labor price (wages)
PDZ  = data[:PDZ]           # inital commodity price (price of domestically-produced commodities)
IOZ = [5. 40. ; 15. 20.]    # intermediate inputs
XDZ = (PKZ.*KZ)./PDZ + (PLZ.*LZ)./PDZ + [sum(IOZ[1:2,]),sum(IOZ[3:4,])].*PDZ
                            # initial commodity production level
io = IOZ./XDZ               # Technical coefficients
sigmaF  = data[:sigmaF]     # elasticity of substition between factors in production function
gammaF = 1./(1.+((PLZ/PKZ)*((KZ./LZ).^(-1./sigmaF))))
                            # distribution parameter of capital
aF = XDZ ./ (gammaF.*KZ.^((sigmaF-1)./sigmaF)+(1-gammaF).*LZ.^((sigmaF-1)./sigmaF)).^(sigmaF./(sigmaF-1))
                            # efficiency parameter in production function

## Consumption block

CZ  = data[:CZ]                         # household demand for commodities
YZ = PKZ*KSZ + PLZ*LSZ			            # household's total income
frisch = -1.1			                      # Frisch parameter
elasY  = data[:elasY]                   # income elasticities of demand for commodities
alphaHLES = elasY.*((PDZ.*CZ)./YZ)
alphaHLES = alphaHLES./sum(alphaHLES)   # marginal budget shares of utility function
muH = CZ + (alphaHLES.*YZ)./(PDZ.*frisch) 	            # Subsistence level
UZ = prod((CZ - muH).^alphaHLES)	      # utility level of household


### VARIABLES ###

sector = [1,2]
com = [1,2]

## Production block

@variables M begin
  K[i = sector], (start = KZ[i])
  PK, (start = PKZ)
  L[i = sector], (start = LZ[i])
  PL, (start = PLZ)
  PD[i = sector], (start = PDZ[i])
  KS, (start = KSZ)
  LS, (start = LSZ)
  XD[i = sector], (start = XDZ[i])
end

## Consumption block

@variables M begin
  C[i = sector], (start = CZ[i])
  Y, (start = YZ)
  U, (start = UZ)
end


### EQUATIONS (constraints) ###

## Consumption

@NLconstraints M begin

  # HOUSEHOLDS
  EQC[i = sector], C[i]  == muH[i]  + alphaHLES[i] /(PD[i]  * (Y - sum(PD[j]  * muH[j] for j in sector )))    	 #consumer consumption
  EQY, Y == PK*KS + PL*LS	                                 #income balance
  EQU[i = sector], U == prod((C[j]  - muH[j])^alphaHLES[j] for j in sector)	                     #household utility
  
  # MARKET CLEARING
  EQXD[i = sector, j = com], XD[i]  == sum(io[i,j]  * XD[j] for i in sector for j in com) + C[i]	                            #market clearing consumption

end

## Production
@NLconstraints M begin
  # FIRMS
  EQK[i = sector], K[i] 	== gammaF[i] ^sigmaF[i]  * PK[i]^(-sigmaF[i] ) * (gammaF[i] ^sigmaF[i]  * PK[i]^(1-sigmaF[i]) +
        (1-gammaF[i] )^sigmaF[i]  * PL[i]^(1-sigmaF[i] ))^(sigmaF[i] /(1-sigmaF[i] )) * (XD[i] /aF[i] )                #firm demand for capital

  EQL[i = sector], L[i]  == (1-gammaF[i] )^sigmaF[i]  * PL[i]^(-sigmaF[i] ) * (gammaF[i] ^sigmaF[i]  * PK[i]^(1-sigmaF[i] ) +
        (1-gammaF[i] )^sigmaF[i]  * PL[i]^(1-sigmaF[i] )^(sigmaF[i] /(1-sigmaF[i] )) * (XD[i] /aF[i] )                 #firm demand for labor

  EQZPC[i = sector], PD[i]  * XD[i]  == PK[i]*K[i]  + PL[i]*L[i]  + sum(io[j] *PD[j] for j in sector)*XD[i]                                          #zero-profit condition

  # MARKET CLEARING
  EQKS[i = sector], sum(K[j] for j in sector) == KS[i]           #capital market clearing  == sum(K)
  EQLS[i = sector], sum(L[j] for j in sector) == LS[i]           #labor market clearing == sum(L)
end

### SOLVER ###

#TRICK == 0  #objective function -> unnecessary because using solver
solve(M)
