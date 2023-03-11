# SimulationBTL0

## Single trial (with plots)

```
$> root
root [0] .L SimuBT.C+
Info in <TMacOSXSystem::ACLiC>: creating shared library /Users/bozzo/Library/CloudStorage/OneDrive-IstitutoNazionalediFisicaNucleare/Work/AMS-L0/SlideEVarie/./SimuBT_C.so
ld: warning: directory not found for option '-L/opt/homebrew/gfortran/lib/'
ld: warning: -undefined dynamic_lookup may not work with chained fixups
ld: warning: directory not found for option '-L/opt/homebrew/gfortran/lib/'
ld: warning: -undefined dynamic_lookup may not work with chained fixups
root [1] Single()
DUT X) sigma measured: 10.969250, smearing measured: 4.511834
DUT X) resolution: 9.998389
```

## Many trials (only final plot)
$> root
root [0] .L SimuBT.C+
Info in <TMacOSXSystem::ACLiC>: creating shared library /Users/bozzo/Library/CloudStorage/OneDrive-IstitutoNazionalediFisicaNucleare/Work/AMS-L0/SlideEVarie/./SimuBT_C.so
ld: warning: directory not found for option '-L/opt/homebrew/gfortran/lib/'
ld: warning: -undefined dynamic_lookup may not work with chained fixups
ld: warning: directory not found for option '-L/opt/homebrew/gfortran/lib/'
ld: warning: -undefined dynamic_lookup may not work with chained fixups
root [1] Full()
0) (0.009998 -0.010000)/0.010000 = -0.000161
1) (0.009982 -0.010000)/0.010000 = -0.001769
2) (0.009979 -0.010000)/0.010000 = -0.002115
3) (0.010007 -0.010000)/0.010000 = 0.000700
...
```