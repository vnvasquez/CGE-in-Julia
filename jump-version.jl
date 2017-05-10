using Ipopt, JuMP, DataFrames

M = Model(solver = IpoptSolver(print_level = 0))
data = readtable("data.csv")

### PARAMETERS ###

## Production block

KZ = [60.,50.]              # inital capital demand
KSZ = sum(KZ)               # capital endowment
LZ  = [20.,90.]              # initial labor demand
LSZ = sum(LZ)               # labor endowment
PKZ = 1.                    # initial capital price
PLZ = 1.                    # initial labor price (wages)
PDZ  = [1.,1.]               # inital commodity price (price of domestically-produced commodities)
IOZ = [5. 40. ; 15. 20.]    # intermediate inputs
XDZ = (PKZ.*KZ)./PDZ + (PLZ.*LZ)./PDZ + [sum(IOZ[1:2,]),sum(IOZ[3:4,])].*PDZ
                            # initial commodity production level
io = IOZ./XDZ               # Technical coefficients
sigmaF  = [0.8, 1.2]         # elasticity of substition between factors in production function
gammaF = 1./(1.+((PLZ/PKZ)*((KZ./LZ).^(-1./sigmaF))))
                            # distribution parameter of capital
aF = XDZ ./ (gammaF.*KZ.^((sigmaF-1)./sigmaF)+(1-gammaF).*LZ.^((sigmaF-1)./sigmaF)).^(sigmaF./(sigmaF-1))
                            # efficiency parameter in production function


## Consumption block (household = consumer and provider of labor)
CZ  = [55.,165.]                          # household demand for commodities
YZ = PKZ*KSZ + PLZ*LSZ			             # household's total income
frisch = -1.1			                       # Frisch parameter
elasY  = [0.9,1.1]                        # income elasticities of demand for commodities
alphaHLES = elasY.*((PDZ.*CZ)./YZ)
alphaHLES = alphaHLES./sum(alphaHLES) # marginal budget shares of utility function
muH = CZ + (alphaHLES.*YZ)./(PDZ.*frisch) 	            # Subsistence level
UZ = prod((CZ - muH).^alphaHLES)	       # utility level of household


# Variables
sector = [1,2]
## Production
@variables M begin
  K[i = sector], (start = KZ[i])
  PK[i = sector], (start = PKZ[i])
  L[i = sector], (start = LZ[i])
  PL[i = sector], (start = PLZ[i])
  PD[i = sector], (start = PDZ[i])
  KS[i = sector], (start = KSZ[i])
  LS[i = sector], (start = LSZ[i])
  XD[i = sector], (start = XDZ[i])
end

## Consumption
@variables M begin
  C[i = sector], (start = CZ[i])
  Y[i = sector], (start = YZ[i])
  U[i = sector], (start = UZ[i])
end


### EQUATIONS (constraints) ###

## Consumption
@NLconstraints M begin
  # HOUSEHOLDS
  EQC[i = sector], C[i]  == muH[i]  + alphaHLES[i] /(PD[i]  * (Y[i] - sum(PD[j]  * muH[j] for j in sector ))    	 #consumer consumption
  EQY[i = sector], Y[i] == PK[i]*KS[i] + PL[i]*LS[i]	                                 #income balance
  EQU[i = sector], U[i] == prod((C[i]  - muH[i] )^alphaHLES[i])	                     #household utility

  # MARKET CLEARING
  EQXD[i = sector], XD[i]  == sum(io[j]  * XD[j] for j in sector) + C[i]	                            #market clearing consumption

end

## Production
@NLconstraints M begin
  # FIRMS
  EQK, K 	== gammaF ^sigmaF  * PK^(-sigmaF ) * (gammaF ^sigmaF  * PK^(1-sigmaF ) +
        (1-gammaF )^sigmaF  * PL^(1-sigmaF ))^(sigmaF /(1-sigmaF )) * (XD /aF )                #firm demand for capital

  EQL, L  == (1-gammaF )^sigmaF  * PL^(-sigmaF ) * (gammaF ^sigmaF  * PK^(1-sigmaF ) +
        (1-gammaF )^sigmaF  * PL^(1-sigmaF )^(sigmaF /(1-sigmaF )) * (XD /aF )                 #firm demand for labor

  EQZPC, PD  * XD  == PK*K  + PL*L  + sum(io *PD)*XD                                          #zero-profit condition

  # MARKET CLEARING
  EQKS[i = sector], sum(K[j] for j in sector) == KS[i]           #capital market clearing  == sum(K)
  EQLS[i = sector], sum(L[j] for j in sector) == LS[i]           #labor market clearing == sum(L)
end

### SOLVER ###

#TRICK == 0  #objective function -> unnecessary because using solver
solve(M)
