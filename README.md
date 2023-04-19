## LBM Everywhere 
This space is for someone who wants to learn LBM Hands-on.

This covers my journey , solvers developed by me.

I am average , so I beleive anyone can do LBM if provided with proper material for learning

LBM can be applied to lot of areas. 
It is just collision , streaming , equilibrium , moments.

In STR/BGK , collision is between probability distribution functions

In MRT , collision is between raw moments.

In Cascade , collison is between central moments.

In Cumulants , collision is between cumulants.

I would suggest begineer to start with solving diffusion equation.(There is hand calculation available to give more clarity)

Later do isothermal flows like couette flow , flow in channel with various boundary conditions , lid driven cavity ,
flow past single , multiple obstalces of various shapes

Move to non-isothermal flows starting with natural convection/Rayleigh bernard convection , forced convection , mixed convection

Than I moved to aeroacoustics , slip flows , MHD natural convection , Nano fluids , Ferro fluids , electro kinetic flows , EHD , ETHD , conjugate heat transfer, Double diffusive natural convection.

I have also done isothermal , non-isothermal VP-LBM where bounce back on obstacle(no slip BC or in general dirichlet BC) can be avoided using external force.

## Can LBM do turbulent flows or turbulence modelling ?

Yes , It is not as hard as traditional turbulence modelling.

It can do LES , k-epsilon 

## Can LBM do flow of power law fluids ?

Yes , similar to turbulence modelling.
Viscosity/relaxation is computed locally using strain rate tensor (which happens to be second order moment)
advantage is avoiding computing velocity gradients


## Area left to explored ?
Multiphase ,multicomponent flows , flow in porous media , image processing , Battery simulation , combustion , Schrodinger equation ,Fire dynamics.

Writing report for all my work is also on my to do list.
