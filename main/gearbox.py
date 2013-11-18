from math import pi, atan, tan, radians, cos, sin, degrees, sqrt, floor, ceil
from main import involute, arcinvolute


class Tool:
    def __init__(self, ha_p, hf_p, rho_fp, x, rho_ao, delta_ao, nc):
        self.ha_p = ha_p
        self.hf_p = hf_p
        self.rho_fp = rho_fp
        self.c = 0.25
        self.nc = nc
        self.x = x
        self.rho_ao = rho_ao
        self.delta_ao = delta_ao


class Material:
    def __init__(self, sh_limit, sf_limit, brinell, classification, name='', e=206000., poisson=0.3, density=7.83e-6):
        self.sh_limit = sh_limit
        self.sf_limit = sf_limit
        self.type = classification
        self.name = name
        self.e = e
        self.poisson = poisson
        self.density = density
        self.brinell = brinell


class Lubricant:
    def __init__(self, v40, name=''):
        self.name = name
        self.v40 = v40


class Gear:
    def __init__(self, profile, material, z, beta, b, bs, alpha=20, m=1, x=0.0, sr=0, rz=0, precision_grade=6,
                 shaft_diameter=0, schema=0, l=0, s=0, backlash=0, gear_crown=1, gear_condition=1,
                 helix_modiffication=1, favorable_contact=True):

        self.profile = profile
        self.material = material

        self.z = z
        self.beta = beta
        self.alpha = alpha
        self.m = m
        self.x = x
        self.b = b
        self.bs = bs
        self.sr = sr
        self.rz = rz
        self.precision_grade = precision_grade
        self.shaft_diameter = shaft_diameter
        self.schema = schema
        self.l = l
        self.s = s
        self.backlash = backlash

        self.gear_crown = gear_crown
        self.helix_modiffication = helix_modiffication

        if favorable_contact:
            self.favorable_contact = 1
        else:
            self.favorable_contact = 0

        self.gear_condition = gear_condition

        self.alpha_t = degrees(atan(tan(radians(self.alpha)) / cos(radians(self.beta))))
        self.d = self.m * self.z / cos(radians(self.beta))
        self.da = self.m * (self.z / cos(radians(self.beta)) + 2 * (self.profile.ha_p + self.x))
        self.df = self.d - 2 * self.m * (self.profile.ha_p + self.profile.c - self.x)
        self.db = self.d * cos(radians(self.alpha_t))
        self.addendum = self.m * (self.profile.ha_p + self.x)
        self.dedendum = self.m * (self.profile.hf_p - self.x)
        self.h = self.dedendum + self.addendum
        self.rho_f = self.profile.rho_fp * self.m
        self.mt = self.m / cos(radians(self.beta))
        self.p_b = self.m * cos(radians(self.alpha)) * pi
        self.p_n = pi * cos(radians(self.alpha))

        #FIXME self.sn = ((pi/2) + 2 * 0 *self.m*self.x*tan(self.Alpha))
        self.beta_b = degrees(atan(self.db * tan(radians(self.beta)) / self.d))

        self.zn = self.z / (cos(radians(self.beta)) * cos(radians(self.beta_b)) ** 2)

        mc = self.__interval_calc([0, 0.5, 2.0, 3.5, 6.0, 25, 40, 70], self.m)
        dc = self.__interval_calc([0, 5, 20, 50, 125, 280, 560, 1000, 1600, 2500, 4000, 6000, 8000, 10000], self.d)
        bc = self.__interval_calc([0, 4, 10, 20, 40, 80, 160, 250, 400, 650, 1000], self.b)

        f_pt = 0.3 * (mc + 0.4 * sqrt(dc)) + 4.
        f_p = 0.3 * mc + 1.25 * sqrt(dc) + 7.
        f_a = 3.2 * sqrt(mc) + 0.22 * sqrt(dc) + 0.27
        f_beta = 0.1 * sqrt(dc) + 0.63 * sqrt(bc) + 4.2
        f_f_alpha = 2.5 * sqrt(mc) + 0.17 * sqrt(dc) + 0.5
        f_h_alpha = 2 * sqrt(mc) + 0.14 * sqrt(dc) + 0.5
        f_h_beta = 0.07 * sqrt(dc) + 0.45 * sqrt(bc) + 3.

        self.f_pt = self.__q(f_pt, self.precision_grade)
        self.f_p = self.__q(f_p, self.precision_grade)
        self.f_alpha = self.__q(f_a, self.precision_grade)
        self.f_beta = self.__q(f_beta, self.precision_grade)
        self.f_f_alpha = self.__q(f_f_alpha, self.precision_grade)
        self.f_h_alpha = self.__q(f_h_alpha, self.precision_grade)
        self.f_h_beta = self.__q(f_h_beta, self.precision_grade)
        self.f_f_beta = self.__q(f_h_beta, self.precision_grade)
        self.f_h_beta5 = self.__q(f_h_beta, 5)

    def __interval_calc(self, intervals, attr):
        for i in range(1, intervals.__len__()):
            if intervals[i - 1] <= attr < intervals[i]:
                return sqrt(intervals[i] * intervals[i - 1])

    def __q(self, x, Qiso):
        x *= 2 ** (0.5 * (Qiso - 5))
        if x >= 10:
            x = round(x)
        elif 5 <= x < 10:
            if x % 1 <= 0.25 or (0.5 <= x % 1 % 1 <= 0.75):
                x = floor(x * 2) * 0.5
            else:
                x = ceil(x * 2) * 0.5
        else:
            x = round(x, 1)
        return x


class Transmition:
    def __init__(self, lubricant, rpm_gear_one, rpm_gear_two, gear_box_type, n, l, gears, ka, sf_min, sh_min):
        self.rpm_gear_one = rpm_gear_one
        self.rpm_gear_two = rpm_gear_two
        self.ka = ka
        self.sh_min = sh_min
        self.sf_min = sf_min

        self.v40 = lubricant.v40
        self.gear_box_type = gear_box_type

        self.u = rpm_gear_one / rpm_gear_two
        self.n = n
        self.l = l
        self.pairs = []
        self.__calculate_single_pairs(gears)

    def __calculate_single_pairs(self, gears):
        for i in gears:
            self.pairs.append(self.__single_pair(i[0], i[1], self.rpm_gear_one, self.rpm_gear_two))

    def __single_pair(self, gear_one, gear_two, rpm_in, rpm_out):
        if gear_one.m is not gear_two.m:
            raise Exception("the modulus of the two gears most be equals")
        if gear_one.alpha is not gear_two.alpha:
            raise Exception("the pressure angle of the two gears most be equals")
        if gear_one.beta is not gear_two.beta:
            raise Warning("the helix angle of the two gears are different")

        u_real = gear_one.z / gear_two.z
        u = rpm_out / rpm_in
        u_error = abs(1 - (u_real / u)) * 100
        inv = involute(gear_one.alpha_t) + 2 * (gear_one.x + gear_two.x) / (gear_one.z + gear_two.z) * tan(
            radians(gear_one.alpha))
        alpha_wt = arcinvolute(inv)
        a = ((gear_one.z + gear_two.z) * gear_one.m) / (2 * cos(radians(gear_one.beta)))
        aw = a * cos(radians(gear_one.alpha)) / cos(radians(alpha_wt))
        epsilon_alpha = (0.5 * (
            sqrt(gear_one.da ** 2 - gear_one.db ** 2) + sqrt(gear_two.da ** 2 - gear_two.db ** 2)) - a * sin(
            radians(alpha_wt))) / (pi * gear_one.m * cos(radians(gear_one.alpha_t)) / (cos(radians(gear_one.beta))))
        epsilon_beta = gear_one.b * sin(radians(gear_one.beta)) / (gear_one.m * pi)
        epsilon_gama = epsilon_alpha + epsilon_beta
        v = rpm_in * gear_one.d * pi / 60000
        ft = 1000. * self.n * 60000 / (pi * gear_one.d * rpm_in)
        if self.ka * ft / gear_one.b < 100:
            fmt = 100
        else:
            fmt = self.ka * ft / gear_one.b

        return [gear_one, gear_two, u_real, u, u_error, v, ft, fmt, rpm_in, rpm_out, a, aw, alpha_wt, epsilon_alpha,
                epsilon_beta, epsilon_gama]
