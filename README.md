## LBM Everywhere 
This space is for someone who wants to learn LBM Hands-on.

This covers my journey , solvers developed by me.

I am average , so I beleive anyone can do LBM if provided with proper material for learning

LBM can be applied to lot of areas. 
It is just collision , streaming , equilibrium , moments.

In STR/BGK , collision is between probability distribution functions/ single particle velocity distribution function.
All digital particles are relaxed at same rate hence the name Single.

In MRT , collision is between raw moments.
Relaxation is different hence the name multi

In Entropic ,  which is based on H-theorem there are two terms to control collision/relaxation , equilbrium df is modified.

In Cascade , collison is between central moments.

In Cumulants , collision is between cumulants.

In Regularized LBM , collision operator is modified.

To get rid of collision operator (which is reason for instability) researchers have also come up with Lattice kinetic schemes(LKS) , MACLAB where relaxation is fixed to 1 , equilibrium distribution function is altered

I would suggest beginner to start with solving diffusion equation.(There is hand calculation available to give more clarity , algorithm is also explained)

Later do  1d shock wave , isothermal flows like couette flow , flow in channel with various boundary conditions , lid driven cavity ,
flow past single , multiple obstalces of various shapes

Move to non-isothermal flows starting with natural convection/Rayleigh bernard convection , forced convection , mixed convection

Than I moved to aeroacoustics , slip flows(TMAC ,combination of bounce back , specular bc) ,natural convection with (MHD, Nano fluids , Ferro fluids) , electro kinetic flows , EHD , ETHD , conjugate heat transfer, Double diffusive natural convection , phase change/Melting natural convection.

I have also done isothermal , non-isothermal VP-LBM where bounce back on obstacle(no slip BC or in general dirichlet BC) can be avoided using external force.
There is alternate to bounce back BC which is counter slip , it can even handle non-linear robin BC (radiation BC). I have also written codes for this as well

## Currently I have been applying LB models for climate studies during my doctoral research
Available is videos of cylinder turning to cone(just one of my recent works)

## Can LBM do turbulent flows or turbulence modelling ?

Yes , It is not as hard as traditional turbulence modelling.

It can do LES , k-epsilon (RANS) or even hybrid.

## Can LBM do flow of power law fluids ?

Yes , similar to turbulence modelling.
Viscosity/relaxation is computed locally using strain rate tensor (which happens to be second order moment)
advantage is avoiding computing velocity gradients


## Area left to explored ?
Image processing , Battery simulation , combustion , Schrodinger equation ,Fire dynamics.

Writing report , adding all my solvers is also on my to do list.

## Credits 
Who taught me LBM ?
I am nowhere near them.

All credit goes to Prof Dr. Krafczyk (my first mentor in life, never found anyone in my country )

For free materials , replying to emails inspite of holding such top positions

I do miss him , remember all his advice like explaining things to layman , setting time for problems & moving on , focus.

Had hard time in Germany(bad luck due to weather , cooking which I was not used to earlier)

Prof Dr Geier who is one of few researcher with expertise in fire simulation , developer of cumulant , cascade operator

According to him only cumulant collison operator is right way to do LBM (reference lectures)

(His top class FREE lecture series with hands on expereince)
https://www.youtube.com/@tubs-irmb6980

## Things to do (May be by 2024 it will be ready)
Added is 3D diffusion solver which can be used as base to do any 3D simulation using LBM , it uses D3Q15 lattice model.
Attached is link of lectures (did it 1 year before) will continue 
https://www.youtube.com/watch?v=BUL-crp1DHc


