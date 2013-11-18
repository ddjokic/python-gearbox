from main.gearbox import *

from standards.iso import Pitting as iso_pitting
from standards.iso import Bending as iso_bending
from standards.agma import Pitting as agma_pitting
from standards.agma import Bending as agma_bending


Lubricante1 = Lubricant(
    name='Kiruna',
    v40=160
)

Material1 = Material(
    name='AISI 2010',
    classification='NV(nitrocar)',
    sh_limit=1500.,
    sf_limit=460.,
    e=206000.,
    poisson=0.3,
    density=7.83e-6,
    brinell=286.6667
)

herramienta = Tool(
    ha_p=1,
    hf_p=1.25,
    rho_fp=0.38,
    x=0,
    rho_ao=0,
    delta_ao=0,
    nc=10.
)

gear1 = Gear(
    profile=herramienta,
    material=Material1,
    z=22.,
    beta=16.0,
    alpha=20.0,
    m=2.5,
    x=0.0,
    b=34.0,
    bs=34.0,
    sr=0.0,
    rz=3.67,
    precision_grade=6.0,
    shaft_diameter=35.0,
    schema=3.0,
    l=60.0,
    s=15.0,
    backlash=0.017
)

gear2 = Gear(
    profile=herramienta,
    material=Material1,
    z=131.,
    m=2.5,
    beta=16.0,
    alpha=20.0,
    x=0.0,
    b=34.0,
    bs=34.0,
    sr=0.0,
    rz=3.67,
    precision_grade=6.0,
    shaft_diameter=20.0,
    schema=3.0,
    l=60.0,
    s=35.0,
    backlash=-0.017
)

print gear1.da
print gear2.da

pares = [[gear1,gear2]]

transmision = Transmition(
    gears=pares,
    lubricant=Lubricante1,
    rpm_gear_one=1450.0,
    rpm_gear_two=243.5,
    n=40.0,
    l=10000.0,
    gear_box_type=2,
    ka=1.3,
    sh_min=1,
    sf_min=1
)

a = iso_pitting(transmition=transmision)
b = a.calculate()
