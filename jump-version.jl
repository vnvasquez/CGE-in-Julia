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
  EQC[i = sector], C[i]  == muH[i]  + alphaHLES[i] /(PD[i]  * (Y[i] - sum(PD[j]  * muH[j] for j in sector ))    	 #consumer consumption
  EQY, Y == PK*KS + PL*LS	                                 #income balance
  EQU, U == prod((C  - muH )^alphaHLES )	                     #household utility

  # MARKET CLEARING
  EQXD, XD  == sum(io  * XD) + C	                            #market clearing consumption

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
  EQKS, sum(K ) == KS           #capital market clearing  == sum(K)
  EQLS, sum(L ) == LS           #labor market clearing == sum(L)
end

### SOLVER ###

#TRICK == 0  #objective function -> unnecessary because using solver
solve(M)
