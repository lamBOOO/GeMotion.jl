using GFluxx

out = GFluxx.simulate(Pr=0.7, Ra=1E3, n=1.0, 50, linesearch=BackTracking(), levels=(;T=[0.1*i for i=1:10],psi=([0.01,0.05,0.1,0.15] |> x->vcat(x,-x))); turan...)
