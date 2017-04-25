using Mimi

@defcomp Test begin

# Create index for the two model sectors, sec1 and sec2

  sector = Index(sector)

# Declare parameters

## Production block

  KSZ = Parameter() # capital endowment
  KZ = Parameter(index=[sector]) # inital capital demand
  LSZ = Parameter() # labor endowment
  LZ = Parameter(index=[sector]) # initial labor demand
  PKZ = Parameter() # initial capital price
  PLZ = Parameter() # initial labor price (wages)
  PDZ = Parameter(index=[sector]) # inital commodity price
  XDZ = Parameter(index=[sector]) # initial commodity production level
  gammaF = Parameter(index=[sector]) # distribution parameter of capital
  sigmaF = Parameter(index=[sector]) # elasticity of substition between factors in production function
  aF = Parameter(index=[sector]) # efficiency parameter in production function

## Consumption block

  CZ
  UZ
  YZ
  muH
  alphaHLES
  frisch
  elasY
