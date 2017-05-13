# CGE-in-Julia
*Completed in partial fulfillment of ER290A: Computational Methods in Economics*

*Summary:* This repository contains a simple CGE model in Julia based on an EcoMod CGE instruction module. The final product represents a two-sector economy with two firms and one household that maximize profits and utility according to a CES production function and a LES utility function, respectively.  The model uses the JuMP package with Ipopt as a non-linear optimization solver to demonstrate an open-source alternative to GAMS. The GitHub for this project can be accessed at: https://github.com/vnvasquez/CGE-in-Julia

*Longform:* The project herein replicates a computable general equilibrium (CGE) model that is drawn from a Global Economic Modeling Network (EcoMod) training example. It translates this simple non-linear model from the General Algebraic Modeling System (GAMS) into julia code.

The model has two sectors (with one firm each) and one household. The firms each produce a single commodity for their sector. The household is the sole consumer of these two commodities. Firms maximize profit according to a constant elasticity of substitution (CES) production function and the household maximizes utility according to a Linear Expenditure System (LES) utility function. The LES utility function is subject to a household budget constraint.  

Importantly, the model uses the JuMP package and the Ipopt solver in place of the GAMS solvers. Therefore when converting the GAMS configuration into Julia, we used the latter language’s Ipopt and DataFrames packages in addition to JuMP. The DataFrames package enabled us to create and view our results in a manageable format.

For ease of use, we divided the model economy into two main activities – production and consumption - and annotated the code in comment blocks. We selected this block format approach with the intention of generating code that is both easily read and simple to expand.

Parameters were constructed to match the sectoral index used to organize input data. Parameters are organized into comment blocks by production and consumption. Parameter values include both initial values such as capital and labor demand as well as endowment figures for capital and labor. Technical coefficients such as the elasticity of substitution between factors in the production function, and the Frisch parameter – the elasticity of the marginal utility of expenditure – are also added to the model in this subsection of the project, among other elements.

The model specifies a sector index for variables that matches the index in the parameter data. The model features block format to enter variables, taking advantage of the JuMP macro for declaring multiple variables. Each variable within each block is indexed appropriately and connected to its starting parameter values, as needed. Some variables are fixed, including the numeraire of a CGE model (in this case, labor price).

Similarly, constraints – known as equations in the GAMS version of the model – are entered in block format. This furnishes an accessible means of categorizing the relationships between variables according to whether they are in the consumption or production block.

To calibrate the model to other data and/or expand the number of sectors, users can update data in the data files. The index specified before the variable blocks can be updated to represent an expanded number of sectors. Individual parameters, such as the Frisch parameter, can be updated by hand.

Future expansions of this model could separate economic activities into different files and functionalize components of a larger model linked together in a main file (see the Mimi package in julia). This is an alternative form of model organization to the GAMS set-parameter-variable-equation format. Larger CGE models would likely need a more firm organizational structure than the one presented here.

The model successfully replicates baseline data. Users can specify shocks to the economy and check results using the getvalue(variable_name) function.
