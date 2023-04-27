## LBM Everywhere 
This space is for someone who wants to learn LBM Hands-on.

This covers my journey , solvers developed by me.

I am average , so I beleive anyone can do LBM if provided with proper material for learning

LBM can be applied to lot of areas. 
It is just collision , streaming , equilibrium , moments.

In STR/BGK , collision is between probability distribution functions/ single particle velocity distribution function

In MRT , collision is between raw moments.

In Cascade , collison is between central moments.

In Cumulants , collision is between cumulants.

In Regularized LBM , collision operator is modified.

To get rid of collision operator researchers have also come up with Lattice kinetic schemes(LKS) , MACLAB.

I would suggest begineer to start with solving diffusion equation.(There is hand calculation available to give more clarity)

Later do isothermal flows like couette flow , flow in channel with various boundary conditions , lid driven cavity ,
flow past single , multiple obstalces of various shapes

Move to non-isothermal flows starting with natural convection/Rayleigh bernard convection , forced convection , mixed convection

Than I moved to aeroacoustics , slip flows(TMAC ,combination of bounce back , specular bc) ,natural convection with (MHD, Nano fluids , Ferro fluids) , electro kinetic flows , EHD , ETHD , conjugate heat transfer, Double diffusive natural convection.

I have also done isothermal , non-isothermal VP-LBM where bounce back on obstacle(no slip BC or in general dirichlet BC) can be avoided using external force.

## Can LBM do turbulent flows or turbulence modelling ?

Yes , It is not as hard as traditional turbulence modelling.

It can do LES , k-epsilon (RANS) or even hybrid.

## Can LBM do flow of power law fluids ?

Yes , similar to turbulence modelling.
Viscosity/relaxation is computed locally using strain rate tensor (which happens to be second order moment)
advantage is avoiding computing velocity gradients


## Area left to explored ?
Multiphase ,multicomponent flows , flow in porous media , image processing , Battery simulation , combustion , Schrodinger equation ,Fire dynamics.

Writing report for all my work is also on my to do list.

## Credits 
Who taught me LBM ?
I am nowhere near them.

All credit goes to Prof Dr. Krafczyk (my first mentor in life, never found anyone in my country )
For free materials

I do miss him , remember all his advice like explaining things to layman , setting time for problems & moving on , focus.

Had hard time in Germany(bad luck due to weather , cooking which I was not used to earlier)

Prof Dr Geier (for top class FREE lecture series)
https://www.youtube.com/@tubs-irmb6980

