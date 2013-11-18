from main.gearbox import Tool, Gear
from main.pair import Pair

from standards.agma import Pitting as agma_pitting
from standards.agma import Bending as agma_bending

perfil1 = Tool(
    ha_p=1.,
    hf_p=1.25,
    rho_fp=0.38,
    hao=1.4760,
    x=0.,
    rho_ao=0.4092,
    delta_ao=0.0061
)

perfil2 = Tool(
    ha_p=1,
    hf_p=1.25,
    rho_fp=0.38,
    hao=1.4760,
    x=0.,
    rho_ao=0.4092,
    delta_ao=0.0240
)

gear1 = Gear(
    profile=perfil1,
    z=21.,
    beta=15,
    alpha=20,
    m=0.166667,
    x=0.5343,
    b=3.75,
    bs=3.75,
    sr=0.0,
    e=206000.0,
    poisson=0.3,
    sigmaHLimit=1500.0,
    sigmaFLimit=460.0,
    material='NV(nitrocar)',
    rz=3.67,
    hb=286.6667,
    precision_grade=6.0,
    shaft_diameter=35.0,
    schema=3.0,
    l=60.0,
    s=15.0,
    rho=7.83e-6,
    backlash = 0.0240
)

gear2 = Gear(
    profile=perfil2,
    z=86.,
    beta=15,
    alpha=20,
    m=0.166667,
    x=0.,
    b=3.75,
    bs=3.75,
    sr=0.0,
    e=206000.0,
    poisson=0.3,
    sigmaHLimit=1500.0,
    sigmaFLimit=460.0,
    material='NV(nitrocar)',
    rz=3.67,
    hb=286.6667,
    precision_grade=6.0,
    shaft_diameter=35.0,
    schema=3.0,
    l=60.0,
    s=15.0,
    rho=7.83e-6,
    backlash = 0.0240

)

transmision = Pair(
    gear1,
    gear2,
    rpmGearOne=1450.0,
    rpmGearTwo=243.5,
    n=40.0,
    l=10000.0,
    v40=160.0,
    fav=1.0,
    helixModiffication=1.0,
    gearBoxType=2,
    gearCrown=1,
    gearCondition=1
)

agma_picadura = agma_pitting(transmition=transmision, ka=1.3, sHMin=1)
agma_flexion = agma_bending(transmition=transmision, ka=1.3, sFMin=1)

print agma_picadura.calculate()['sigmaH']
print agma_flexion.calculate()['sigmaF']
