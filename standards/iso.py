from math import acos, log, exp

from numpy import interp
from main.gearbox import *


class Pitting:
    def __init__(self, transmition):
        self.transmition = transmition

    def calculate(self):
        results = []
        for i in self.transmition.pairs:
            results.append(self.calculate_by_pair(i))
        return 0

    def calculate_by_pair(self, pair):
        u = pair[2]

        zh = self.__zh(pair)
        zb, zd = self.__zb(pair)
        ze = self.__ze(pair)
        z_epsilon = self.__z_epsilon(pair)
        z_beta = 1 / sqrt(cos(radians(pair[0].beta)))
        znt_one = self.__znt(pair[0].material.type, pair[8])
        znt_two = self.__znt(pair[1].material.type, pair[9])
        zl = self.__zl(pair)
        zv = self.__zv(pair)
        zr = self.__zr(pair)
        zw_one, zw_two = self.__zw(pair)
        zx = 1

        kv = __kv__(pair)
        khb = khb(pair)
        kha = kha(pair)

        sigmaH0 = zh * ze * z_epsilon * z_beta * sqrt((self.Ft * (u + 1)) / (pair.gear_one.d * pair.gear_one.b * u))
        sigmaHOne = zb * sigmaH0 * sqrt(self.ka * kv * khb * kha)
        sigmaHTwo = zd * sigmaH0 * sqrt(self.ka * kv * khb * kha)

        sigmaHPOne = pair.gear_one.material.sh_limit * znt_one * zl * zv * zr * zw_one * zx / self.sh_min
        sigmaHPTwo = pair.gear_two.material.sh_limit * znt_two * zl * zv * zr * zw_two * zx / self.sh_min

        sHOne = znt_one * zl * zv * zr * zw_one * zx / sigmaHOne
        sHTwo = znt_two * zl * zv * zr * zw_two * zx / sigmaHTwo

        return {
            'sigmaHOne': sigmaHOne,
            'sigmaHTwo': sigmaHTwo,
            'sigmaHPOne': sigmaHPOne,
            'sigmaHPTwo': sigmaHPTwo,
            'zh': zh,
            'zb': zb,
            'zd': zd,
            'ze': ze,
            'z_epsilon': z_epsilon,
            'z_beta': z_beta,
            'znt_one': znt_one,
            'znt_two': znt_two,
            'zl': zl,
            'zv': zv,
            'zr': zr,
            'zw_one': zw_one,
            'zw_two': zw_two,
            'zx': zx,
            'kv': kv,
            'khb': khb,
            'kha': kha
        }

    def __zh(self,pair):
        beta_b = radians(pair[0].beta_b)
        alpha_wt = radians(pair[12])
        alpha_t = radians(pair[0].alpha_t)
        return sqrt((2. * cos(beta_b) * cos(alpha_wt)) / (cos(alpha_t) ** 2. * sin(alpha_wt)))

    def __zb(self,pair):
        da_one = pair[0].da
        db_one = pair[0].db
        z_one = pair[0].z
        beta = pair[0].beta
        da_two = pair[1].da
        db_two = pair[1].db
        z_two = pair[1].z
        alpha_wt = radians(pair[12])
        epsilon_alpha = pair[13]
        epsilon_beta = pair[14]
        Zb = 0
        Zd = 0

        M1 = tan(alpha_wt) / sqrt((sqrt((da_one ** 2 / db_one ** 2) - 1) - (2 * pi) / z_one) * (
            sqrt((da_two ** 2 / db_two ** 2) - 1) - (epsilon_alpha - 1) * (2 * pi) / z_two))
        M2 = tan(alpha_wt) / sqrt((sqrt((da_two ** 2 / db_two ** 2) - 1) - (2 * pi) / z_two) * (
            sqrt((da_one ** 2 / db_one ** 2) - 1) - (epsilon_alpha - 1) * (2 * pi) / z_one))

        if beta is 0:
            if epsilon_alpha > 1:
                if M1 <= 1:
                    Zb = 1
                else:
                    Zb = M1
                if M2 <= 1:
                    Zd = 1
                else:
                    Zd = M2
            if (z_one / z_two) > 1.5:
                Zd = 1
        else:
            if epsilon_alpha > 1 and epsilon_beta >= 1:
                Zb = 1
                Zd = Zb
            elif epsilon_alpha > 1 > epsilon_beta:
                if M1 - epsilon_beta * (M1 - 1) < 1:
                    Zb = 1
                else:
                    Zb = M1 - epsilon_beta * (M1 - 1)
                if M2 - epsilon_beta * (M2 - 1) < 1:
                    Zd = 1
                else:
                    Zd = M2 - epsilon_beta * (M2 - 1)

        return Zb, Zd

    def __ze(self,pair):
        e_one = pair[0].material.e
        e_two = pair[1].material.e
        poisson_one = pair[0].material.poisson
        poisson_two = pair[1].material.poisson

        if e_one == e_two:
            return sqrt(e_one / (2 * pi * (1 - poisson_one ** 2)))
        else:
            return sqrt(1 / (pi * (((1 - poisson_one) / e_one) + (1 - poisson_two) / e_two)))

    def __z_epsilon(self,pair):
        epsilon_alpha = pair[13]
        epsilon_beta = pair[14]

        if epsilon_beta is 0:
            return sqrt((4 - epsilon_alpha) / 3)
        elif 1 > epsilon_beta > 0:
            return sqrt(((4 - epsilon_alpha) / 3) * (1 - epsilon_beta) + (epsilon_beta / epsilon_alpha))
        else:
            return sqrt(1 / epsilon_alpha)

    def __znt(self, material, rpm):

        nl = self.transmition.l * 60 * rpm

        if material == 'NV(nitrocar)':
            Y = [1.1, 1.1, 1.02, 1, 0.97, 0.93, 0.89, 0.85]
        elif material == 'GG' or material == 'GGG(ferr)' or material == 'NT' or material == 'NV(nitr)':
            Y = [1.3, 1.3, 1.07, 1, 0.97, 0.93, 0.89, 0.85]
        else:
            Y = [1.6, 1.6, 1.36, 1.14, 1, 0.98, 0.915, 0.85]
        X = [1e4, 1e5, 1e6, 2e6, 1e7, 1e8, 1e9, 1e10]

        return interp(nl, X, Y)

    def __r_red(self,pair):
        db_one = pair[0].db
        db_two = pair[1].db
        alpha_wt = radians(pair[12])
        r1 = 0.5 * db_one * tan(alpha_wt)
        r2 = 0.5 * db_two * tan(alpha_wt)

        return (r1 * r2) / (r1 + r2)

    def __czl(self,pair):
        sh_limit_1 = pair[0].material.sh_limit
        sh_limit_2 = pair[1].material.sh_limit
        if sh_limit_1 < sh_limit_2:
            sh_min = sh_limit_1
        else:
            sh_min = sh_limit_2

        if 850 <= sh_min <= 1200:
            czl = (sh_min / 437.5) + 0.6357
        elif 850 > sh_min:
            czl = 0.83
        else:
            czl = 0.91
        return czl

    def __zl(self,pair):
        czl = self.__czl(pair)
        v40 = self.transmition.v40
        return czl + 4.0 * (1.0 - czl) / (1.2 + 134.0 / v40) ** 2

    def __zv(self,pair):
        v = pair[5]
        czv = self.__czl(pair) + 0.02
        return czv + (2 * (1 - czv) / sqrt(0.8 + 32 / v))

    def __zr(self,pair):
        rz_one = pair[0].rz
        rz_two = pair[1].rz
        sh_limit_1 = pair[0].material.sh_limit
        sh_limit_2 = pair[0].material.sh_limit
        czr = 0

        if sh_limit_1 < sh_limit_2:
            sh_min = sh_limit_1
        else:
            sh_min = sh_limit_2

        rz = (rz_one + rz_two) / 2.

        Rz10 = rz * ((10.0 / self.__r_red(pair)) ** (1. / 3.))

        if 850 <= sh_min <= 1200:
            czr = 0.32 - 0.0002 * sh_min
        if 850 > sh_min:
            czr = 0.15
        if 1200 < sh_min:
            czr = 0.08

        return (3.0 / Rz10) ** czr

    def __zw(self, pair):
        rz_one = pair[0].rz
        rz_two = pair[1].rz
        hb_1 = pair[0].material.brinell
        hb_2 = pair[1].material.brinell
        v40 = self.transmition.v40
        v = pair[5]

        rzH = ((rz_one * (10 / self.__r_red(pair)) ** 0.33) * (rz_one / rz_two) ** 0.66) / ((v40 * v / 1500) ** 0.33)

        if rzH > 16:
            rzH = 16
        if rzH < 3:
            rzH = 3

        if hb_1 < 130:
            return 1.2 * (3 / rzH) ** 0.15
        elif hb_1 > 470:
            return (3 / rzH) ** 0.15

        if hb_2 < 130:
            return 1.2 * (3 / rzH) ** 0.15
        elif hb_2 > 470:
            return (3 / rzH) ** 0.15

        return (1.2 - (hb_1 - 130) / 1700) * (3 / rzH) ** 0.15,(1.2 - (hb_2 - 130) / 1700) * (3 / rzH) ** 0.15


#iso 6336-3
class Bending:
    def __init__(self, ka, sFMin, pair):
        self.ka = ka
        self.sFMin = sFMin
        self.pair = pair
        self.Ft = 1000. * pair.n * 60000 / (pi * pair.gear_one.d * pair.rpmgear_one)
        if self.ka * self.Ft / pair.gear_one.b < 100:
            self.fmt = 100
        else:
            self.fmt = ka * self.Ft / pair.gear_one.b

    def calculate(self):
        gear_one = self.pair.gear_one
        gear_two = self.pair.gear_two
        pair.gear_one.b = gear_one.b
        bTwo = gear_two.b
        m = gear_one.m
        sigmaFLimitOne = gear_one.sigmaFLimit
        sigmaFLimitTwo = gear_two.sigmaFLimit

        Yst = self.__yst()
        YxOne = self.__yx(gear_one)
        YxTwo = self.__yx(gear_two)
        YfOne, YfTwo = self.__yf()
        YsOne, YsTwo = self.__ys()
        Ybeta = self.__ybeta()
        YbOne = self.__yb(gear_one)
        YbTwo = self.__yb(gear_two)

        YdeltaOne, YdeltaTwo = self.__ydelta()
        Ydt = self.__ydt()
        YntOne, YntTwo = self.__ynt()
        YRelOne = self.__yrel(gear_one)
        YRelTwo = self.__yrel(gear_two)

        kv = __kv(self.pair, self.fmt)
        kfa = __kfa(self.pair, self.fmt)
        kfb = __kfb(self.pair, self.fmt)

        sigmaF0One = self.Ft / (pair.gear_one.b * m) * YfOne * YsOne * Ybeta * YbOne * Ydt
        sigmaF0Two = self.Ft / (bTwo * m) * YfTwo * YsTwo * Ybeta * YbTwo * Ydt

        sigmaFOne = sigmaF0One * self.ka * kv * kfb * kfa
        sigmaFTwo = sigmaF0Two * self.ka * kv * kfb * kfa

        sigmaFPOne = sigmaFLimitOne * Yst * YntOne * YdeltaOne * YRelOne * YxOne / self.sFMin
        sigmaFPTwo = sigmaFLimitTwo * Yst * YntTwo * YdeltaTwo * YRelTwo * YxTwo / self.sFMin

        sFOne = sigmaFLimitOne * YsOne * YntOne * YdeltaOne * YRelOne / sigmaFOne
        sFTwo = sigmaFLimitTwo * YsTwo * YntTwo * YdeltaTwo * YRelTwo / sigmaFTwo

        return {
            'sigmaFOne': sigmaFOne,
            'sigmaFTwo': sigmaFTwo,
            'sigmaFPOne': sigmaFPOne,
            'sigmaFPTwo': sigmaFPTwo,
            'Yst': Yst,
            'YxOne': YxOne,
            'YxTwo': YxTwo,
            'YfOne': YfOne,
            'YfTwo': YfTwo,
            'YsOne': YsOne,
            'YsTwo': YsTwo,
            'Ybeta': Ybeta,
            'YbOne': YbOne,
            'YbTwo': YbTwo,
            'YdeltaOne': YdeltaOne,
            'YdeltaTwo': YdeltaTwo,
            'Ydt': Ydt,
            'YntOne': YntOne,
            'YntTwo': YntTwo,
            'YRelOne': YRelOne,
            'YRelTwo': YRelTwo,
            'kv': kv,
            'kfa': kfa,
            'kfb': kfb
        }

    def __yrel(self, gear):
        rz = gear.rz
        material = gear.material
        if rz < 1:
            if material == 'V' or material == 'GGG(perl)' or material == 'Eh' or material == 'IF' or material == 'GTS':
                return 1.12
            elif material == 'St':
                return 1.07
            elif material == 'GG' or material == 'GGG(ferr)' or material == 'NT' or material == 'NV(nitrocar)' or material == 'NV(nitr)':
                return 1.025
        else:
            if material == 'V' or material == 'GGG(perl)' or material == 'Eh' or material == 'IF' or material == 'GTS':
                return 1.674 - 0.529 * (rz + 1) ** 0.1
            elif material == 'St':
                return 5.306 - 4.203 * (rz + 1) ** 0.01
            elif material == 'GG' or material == 'GGG(ferr)' or material == 'NT' or material == 'NV(nitrocar)' or material == 'NV(nitr)':
                return 4.299 - 3.259 * (rz + 1) ** 0.0058

    def __ynt(self):
        materialOne = self.pair.gear_one.material
        materialTwo = self.pair.gear_two.material
        rpmOne = self.pair.rpmgear_one
        rpmTwo = self.pair.rpmgear_two
        L = self.pair.l
        Y = 0

        nlOne = L * 60 * rpmOne
        nlTwo = L * 60 * rpmTwo

        result = []

        for material, nl in [materialOne, nlOne], [materialTwo, nlTwo]:

            if material == 'V' or material == 'GGG(perl)' or material == 'GTS' or material == 'St':
                Y = [2.5, 2.5, 2.5, 1.79, 1.2, 1, 0.98, 0.94, 0.895, 0.85]

            elif material == 'Eh' or material == 'IF':
                Y = [2.5, 2.5, 1.97, 1.53, 1.15, 1, 0.98, 0.94, 0.895, 0.85]

            elif material == 'GG' or material == 'GGG(ferr)' or material == 'NT' or material == 'NV(nitr)':
                Y = [1.6, 1.6, 1.4, 1.21, 1.07, 1, 0.98, 0.94, 0.895, 0.85]

            elif material == 'NV(nitrocar)':
                Y = [1.1, 1.1, 1.07, 1.05, 1.01, 1, 0.98, 0.94, 0.895, 0.85]

            X = [1e2, 1e3, 1e4, 1e5, 1e6, 3e6, 1e7, 1e8, 1e9, 1e10]
            result.append(interp(nl, X, Y))

        return result

    def __ydelta(self):
        gear_one = self.pair.gear_one
        gear_two = self.pair.gear_two

        alphaFenOne, alphaFenTwo, sFnOne, sFnTwo, rhoFOne, rhoFTwo, hFeOne, hFeTwo = self.__aux()
        rhoOne = self.__rho(gear_one)
        rhoTwo = self.__rho(gear_two)

        qsOne = sFnOne / (2 * rhoFOne)
        qsTwo = sFnTwo / (2 * rhoFTwo)
        xpOne = (1 / 5.) * (1 + 2 * qsOne)
        xpTwo = (1 / 5.) * (1 + 2 * qsTwo)
        xt = (1 / 5.) * (1 + 2 * 2.5)

        ydeltaOne = (1 + sqrt(rhoOne * xpOne)) / (1 + sqrt(rhoOne * xt))
        ydeltaTwo = (1 + sqrt(rhoTwo * xpTwo)) / (1 + sqrt(rhoTwo * xt))

        return ydeltaOne, ydeltaTwo

    def __rho(self, gear):
        sigmaFLimit = gear.sigmaFLimit
        material = gear.material
        rho = 0

        if material == 'GG' or material == 'GGG(ferr)':
            if sigmaFLimit <= 150:
                rho = 0.3124
            elif sigmaFLimit >= 300:
                rho = 0.3095
            else:
                rho = interp(sigmaFLimit, [150, 300], [0.3124, 0.3095])

        elif material == 'GGG(perl)' or material == 'V' or material == 'GTS':
            if sigmaFLimit <= 500:
                rho = 0.0281
            elif sigmaFLimit >= 1000:
                rho = 0.0014
            else:
                rho = interp(sigmaFLimit, [500, 600, 800, 1000], [0.0281, 0.0194, 0.0064, 0.0014])

        elif material == 'NT' or material == 'NV(nitrocar)' or material == 'NV(nitr)':
            rho = 0.1005

        elif material == 'St':
            if sigmaFLimit <= 300:
                rho = 0.0833
            elif sigmaFLimit >= 400:
                rho = 0.0445
            else:
                rho = interp(sigmaFLimit, [300, 400], [0.0833, 0.445])

        elif material == 'Eh' or material == 'IF':
            rho = 0.003

        return rho

    def __ydt(self):
        epsilonAlpha = self.pair.epsilonAlpha
        betaB = radians(self.pair.gear_one.betaB)
        gPrecision = self.pair.gear_one.gPrecision
        epsilonAlphaN = epsilonAlpha / (cos(betaB) ** 2)

        if epsilonAlphaN <= 2.05 or epsilonAlphaN > 2.05 and gPrecision > 4:
            return 1
        elif 2.05 < epsilonAlphaN <= 2.5 and gPrecision <= 4:
            return -0.666 * epsilonAlphaN + 2.366
        elif epsilonAlphaN > 2.5 and gPrecision <= 4:
            return 0.9

    def __yb(self, gear):
        sr = gear.sr
        h = gear.h
        srH = sr / h
        m = self.pair.gear_one.m

        if srH >= 1.2 or sr / m <= 1.75:
            return 1
        elif 0.5 < srH < 1.2:
            return 1.6 * log(2.242 * 1 / srH)

    def __ybeta(self):
        epsilonBeta = self.pair.epsilonBeta
        beta = self.pair.gear_one.beta

        if epsilonBeta > 1:
            eb = 1
        else:
            eb = epsilonBeta

        if beta > 30:
            be = 30
        else:
            be = beta

        return 1 - eb * (be / 120)

    def __ys(self):
        alphaFenOne, alphaFenTwo, sFnOne, sFnTwo, rhoFOne, rhoFTwo, hFeOne, hFeTwo = self.__aux()
        LOne = sFnOne / hFeOne
        LTwo = sFnTwo / hFeTwo
        qsOne = sFnOne / (2 * rhoFOne)
        qsTwo = sFnTwo / (2 * rhoFTwo)
        ysOne = (1.2 + 0.13 * LOne) * (qsOne ** (1 / (1.21 + (2.3 / LOne))))
        ysTwo = (1.2 + 0.13 * LTwo) * (qsTwo ** (1 / (1.21 + (2.3 / LTwo))))

        return ysOne, ysTwo

    def __aux(self):
        m = self.pair.gear_one.m
        epsilonAlpha = self.pair.epsilonAlpha
        znOne = self.pair.gear_one.zn
        znTwo = self.pair.gear_two.zn
        zTwo = self.pair.gear_two.z
        xOne = self.pair.gear_one.x
        xTwo = self.pair.gear_two.x
        zOne = self.pair.gear_one.z
        beta = radians(self.pair.gear_one.beta)
        betaB = radians(self.pair.gear_one.betaB)
        alpha = radians(self.pair.gear_one.alpha)
        hfP = self.pair.gear_one.profile.hfP * m
        haP = self.pair.gear_one.profile.haP
        rhoFP = self.pair.gear_one.profile.RhoFP * m
        thetaOne = 0
        thetaTwo = 0

        e = ((pi / 4) * m) - hfP * tan(alpha) - (1 - sin(alpha)) * (rhoFP / cos(alpha))
        gOne = (rhoFP / m) - (hfP / m) + xOne
        gTwo = (rhoFP / m) - (hfP / m) + xTwo
        hOne = (2 / znOne) * ((pi / 2) - (e / m)) - (pi / 3)
        hTwo = (2 / znTwo) * ((pi / 2) - (e / m)) - (pi / 3)

        thetaTOne = pi / 6
        thetaTTwo = pi / 6
        for i in range(1, 10):
            thetaOne = ((2 * gOne) / znOne) * tan(thetaTOne) - hOne
            thetaTwo = ((2 * gTwo) / znTwo) * tan(thetaTTwo) - hTwo
            thetaTOne = thetaOne
            thetaTTwo = thetaTwo

        pair.gear_one.d = m * zOne / cos(beta)
        dTwo = m * zTwo / cos(beta)
        daOne = m * (zOne / cos(beta) + 2. * (haP + xOne))
        daTwo = m * (zTwo / cos(beta) + 2. * (haP + xTwo))

        epsilonAlphaN = epsilonAlpha / (cos(betaB) ** 2)
        dnOne = pair.gear_one.d / (cos(betaB) ** 2)
        dnTwo = dTwo / (cos(betaB) ** 2)
        dbnOne = dnOne * cos(alpha)
        dbnTwo = dnTwo * cos(alpha)
        danOne = dnOne + daOne - pair.gear_one.d
        danTwo = dnTwo + daTwo - dTwo
        denOne = 2. * sqrt((sqrt((danOne / 2.) ** 2 - (dbnOne / 2.) ** 2) - (
            (pi * pair.gear_one.d * cos(beta) * cos(alpha)) / zOne) * (epsilonAlphaN - 1)) ** 2 + (dbnOne / 2) ** 2)
        denTwo = 2. * sqrt((sqrt((danTwo / 2.) ** 2 - (dbnTwo / 2.) ** 2) - (
            (pi * dTwo * cos(beta) * cos(alpha)) / zTwo) * (epsilonAlphaN - 1)) ** 2 + (dbnTwo / 2) ** 2)
        alphaEnOne = acos(dbnOne / denOne)
        alphaEnTwo = acos(dbnTwo / denTwo)
        gamaEOne = ((0.5 * pi + 2. * tan(alpha) * xOne) / znOne) + degrees(involute(alpha)) - degrees(involute(
            alphaEnOne))
        gamaETwo = ((0.5 * pi + 2. * tan(alpha) * xTwo) / znTwo) + degrees(involute(alpha)) - degrees(involute(
            alphaEnTwo))

        alphaFenOne = alphaEnOne - gamaEOne
        alphaFenTwo = alphaEnTwo - gamaETwo

        sFnOne = m * (znOne * sin((pi / 3) - thetaOne) + sqrt(3) * ((gOne / cos(thetaOne)) - (rhoFP / m)))
        sFnTwo = m * (znTwo * sin((pi / 3) - thetaTwo) + sqrt(3) * ((gTwo / cos(thetaTwo)) - (rhoFP / m)))

        rhoFOne = m * (
            rhoFP / m + ((2 * gOne ** 2) / (cos(thetaOne) * (znOne * cos(thetaOne) ** 2 - 2 * gOne))))
        rhoFTwo = m * (
            rhoFP / m + ((2 * gTwo ** 2) / (cos(thetaTwo) * (znTwo * cos(thetaTwo) ** 2 - 2 * gTwo))))

        hFeOne = 0.5 * m * (
            (cos(gamaEOne) - sin(gamaEOne) * tan(alphaFenOne)) * denOne / m - znOne * cos(pi / 3 - thetaOne) - (
                gOne / cos(thetaOne) - rhoFP / m))
        hFeTwo = 0.5 * m * (
            (cos(gamaETwo) - sin(gamaETwo) * tan(alphaFenTwo)) * denTwo / m - znTwo * cos(pi / 3 - thetaTwo) - (
                gTwo / cos(thetaTwo) - rhoFP / m))

        return alphaFenOne, alphaFenTwo, sFnOne, sFnTwo, rhoFOne, rhoFTwo, hFeOne, hFeTwo

    def __yf(self):
        alpha = radians(self.pair.gear_one.alpha)
        alphaFenOne, alphaFenTwo, sFnOne, sFnTwo, rhoFOne, rhoFTwo, hFeOne, hFeTwo = self.__aux()
        m = self.pair.gear_one.m

        yfOne = ((6 * hFeOne / m) * cos(alphaFenOne)) / (((sFnOne / m) ** 2) * cos(alpha))
        yfTwo = ((6 * hFeTwo / m) * cos(alphaFenTwo)) / (((sFnTwo / m) ** 2) * cos(alpha))

        return yfOne, yfTwo

    def __yx(self, gear):
        material = gear.material
        m = gear.m
        X = 0
        Y = 0

        if material == 'St' or material == 'V' or material == 'GGG(perl)' or material == 'GTS':
            Y = [1, 1, 0.85, 0.85]
            X = [0, 5, 30, 60]
        elif material == 'GG' or material == 'GGG(ferr)':
            Y = [1, 1, 0.85, 0.7, 0.7]
            X = [0, 5, 15, 25, 60]
        elif material == 'NV(nitrocar)' or material == 'NT' or material == 'NV(nitr)' or material == 'Eh' or material == 'IF':
            Y = [1, 1, 0.95, 0.9, 0.85, 0.8, 0.8]
            X = [0, 5, 10, 15, 20, 25, 60]

        return interp(m, X, Y)

    def __yst(self):
        return 2

def __c( pair ):
    x = [0, 0.04723, 0.15551, 0.25791, -0.00635, -0.11654, -0.00193, -0.24188, 0.00529, 0.00182]
    qp = x[1] + (x[2] / pair.gear_one.zn) + ( x[3] / pair.gear_two.zn ) + (
        x[4] * pair.gear_one.x) + ( x[5] * pair.gear_one.x / pair.gear_one.zn ) + (
             x[6] * pair.gear_two.x) + (x[7] * pair.gear_two.x / pair.gear_two.zn) + (
             x[8] * pair.gear_one.x ** 2) + (x[9] * pair.gear_two.x ** 2)

    beq = pair.gear_one.bs / pair.gear_one.b
    if beq < 0.2:
        beq = 0.2
    elif beq > 1.2:
        beq = 1.2

    SrM = pair.gear_one.sr / pair.gear_one.m
    if SrM < 1:
        SrM = 1

    cth = 1 / qp

    CM = 0.8
    CR = 1 + ( log(beq) / ( 5 * exp(SrM / ( 5 * pair.gear_one.m ))))
    CB = 0.975
    cp = cth * CM * CR * CB * cos(radians(pair.gear_one.beta))
    cGammaAlpha = cp * (0.75 * pair.epsilonAlpha + 0.25)
    if pair.epsilonAlpha < 1.2:
        cGammaAlpha *= 0.9

    cGammaBeta = 0.85 * cGammaAlpha

    return cGammaAlpha, cGammaBeta, cp


def __yb( gear, pair, fbX ):
    material = gear.material
    sigmaHLimit = gear.sigmaHLimit
    v = pair.v
    yb = 0

    if material == 'V' or material == 'St' or material == 'GGG(perl)' or material == 'GTS':
        yb = ( 320. / sigmaHLimit) * fbX
        if 5 < v <= 10 and yb > 25600. / sigmaHLimit:
            yb = 25600. / sigmaHLimit
        elif v > 10 and yb > 12800. / sigmaHLimit:
            yb = 12800. / sigmaHLimit

    elif material == 'GG' or material == 'GGG(ferr)':
        yb = 0.55 * fbX
        if 5 < v <= 10 and yb > 45:
            yb = 45
        elif v > 10 and yb > 22:
            yb = 22

    elif material == 'Eh' or material == 'IF' or material == 'NT' or material == 'NV(nitr)' or material == 'NV(nitrocar)':
        yb = 0.15 * fbX

    return yb


def __ya( material, sigmaHLimit, v, fpbOne, fpbTwo ):
    ya = 0
    if fpbOne > fpbTwo:
        fpb = fpbOne
    else:
        fpb = fpbTwo

    if material == 'V' or material == 'St' or material == 'GGG(perl)' or material == 'GTS':
        ya = 160 / (sigmaHLimit * fpb)
        if 5 < v <= 10 and ya > 12800 / sigmaHLimit:
            ya = 12800 / sigmaHLimit
        elif 10 < v and ya > 6400 / sigmaHLimit:
            ya = 6400 / sigmaHLimit

    elif material == 'GG' or material == 'GGG(ferr)':
        ya = 0.275 * fpb
        if 5 < v <= 10 and ya > 12800 / sigmaHLimit:
            ya = 22
        elif 10 < v and ya > 6400 / sigmaHLimit:
            ya = 11

    elif material == 'Eh' or material == 'IF' or material == 'NT' or material == 'NV(nitr)' or material == 'NV(nitrocar)':
        ya = 0.075 * fpb
        if ya > 3:
            ya = 3

    return ya


def __kv__(pair):
    fmt = pair[7]
    sr_one = pair.gear_one.sr
    sr_two = pair.gear_one.sr
    epsilon_gama = pair.epsilonGama
    sh_limit_1 = pair.gear_one.sigmaHLimit
    sh_limit_2 = pair.gear_two.sigmaHLimit
    rho_one = pair.gear_one.rho
    rho_two = pair.gear_two.rho
    da_one = pair.gear_one.da
    da_two = pair.gear_two.da
    df_one = pair.gear_one.df
    df_two = pair.gear_two.df
    db_one = pair.gear_one.db
    rpm_one = pair.rpmgear_one
    f_f_alpha_one = pair.gear_one.FfAlpha
    f_f_alpha_two = pair.gear_two.FfAlpha
    material_one = pair.gear_one.material
    material_two = pair.gear_two.material
    v = pair.v
    precision_grade = pair.gear_one.gPrecision
    u = pair.uReal
    z_one = pair.gear_one.z

    c_gamma_alpha, c_gamma_beta, cp = __c(pair)
    di_one = df_one - 2 * sr_one
    di_two = df_two - 2 * sr_two
    dm_one = (da_one + df_one) / 2
    dm_two = (da_two + df_two) / 2
    q_one = di_one / dm_one
    q_two = di_two / dm_two

    cv1 = 0
    cv2 = 0
    cv3 = 0
    cv4 = 0
    cv5 = 0
    cv6 = 0
    cv7 = 0
    kv = 0

    if sr_one == 0:
        a_one = 1
    else:
        a_one = 1 / ( 1 - q_one ** 4 )

    if sr_two == 0:
        a_two = 1
    else:
        a_two = 1 / ( 1 - q_two ** 4 )

    m_red = ( pi / 8 ) * (( dm_one / db_one ) ** 2 ) * (
        dm_one ** 2 / (( 1 / ( a_one * rho_one )) + ( 1 / (a_two * rho_two * u ** 2))))

    ne_one = (30000 / (pi * z_one)) * sqrt(c_gamma_alpha / m_red)

    n = rpm_one / ne_one

    if fmt == 100:
        ns = 0.5 + 0.35 * sqrt(fmt / 100)
    else:
        ns = 0.85

    if 1 < epsilon_gama <= 2:
        cv1 = 0.32
        cv2 = 0.34
        cv3 = 0.23
        cv4 = 0.9
        cv5 = 0.47
        cv6 = 0.47
    elif epsilon_gama > 2:
        cv1 = 0.32
        cv2 = 0.57 / (epsilon_gama - 0.3)
        cv3 = 0.096 / (epsilon_gama - 1.56)
        cv4 = (0.57 - 0.05 * epsilon_gama) / (epsilon_gama - 1.44)
        cv5 = 0.47
        cv6 = 0.12 / (epsilon_gama - 1.74)

    if 1 < epsilon_gama <= 1.5:
        cv7 = 0.75
    elif 1.5 < epsilon_gama <= 2.5:
        cv7 = 0.125 * sin(pi * (epsilon_gama - 2)) + 0.875
    elif epsilon_gama > 2.5:
        cv7 = 1

    cay_one = 1. / 8 * (((pair.gear_one.material.sh_limit / 97.) - 18.45) ** 2) + 1.5
    cay_two = 1. / 8 * (((pair.gear_two.material.sh_limit / 97.) - 18.45) ** 2) + 1.5
    cay = 0.5 * (cay_one + cay_two)
    ya_one = __ya(material_one, pair.gear_one.material.sh_limit, v, f_f_alpha_one, f_f_alpha_two)
    ya_two = __ya(material_two, pair.gear_two.material.sh_limit, v, f_f_alpha_one, f_f_alpha_two)
    ya = 0.5 * (ya_two + ya_one)

    #FIXME
    if f_f_alpha_one > f_f_alpha_two:
        f_f_alpha = f_f_alpha_one
        fpb = f_f_alpha_one
    else:
        f_f_alpha = f_f_alpha_two
        fpb = f_f_alpha_two

    fpbeff = fpb - ya
    ffaeff = f_f_alpha - ya
    bp = cp * fpbeff / fmt
    bf = cp * ffaeff / fmt
    #FIXME

    if precision_grade >= 6:
        bk = 1
    else:
        bk = abs(1 - (cp * cay / fmt))

    if n <= ns:
        k = cv1 * bp + cv2 * bf + cv3 * bk
        kv = n * k + 1
    elif ns < n <= 1.15:
        kv = cv1 * bp + cv2 * bf + cv4 * bk + 1
    elif n >= 1.5:
        kv = cv5 * bp + cv6 * bf + cv7
    elif 1.15 < n < 1.5:
        kv = ( cv5 * bp + cv6 * bf + cv7 ) + (
            (( cv1 * bp + cv2 * bf + cv4 * bk + 1) - ( cv5 * bp + cv6 * bf + cv7)) / 0.35) * (1.5 - n)

    return kv


def khb(pair):
    fmt = pair[7]
    b = pair.gear_one.b
    helixModiffication = pair.helixModiffication
    d = pair.gear_one.d
    shaftDiameter = pair.gear_one.shaftDiameter
    schema = pair.gear_one.schema
    l = pair.gear_one.l
    s = pair.gear_one.s
    FhbOne = pair.gear_one.FhBeta
    FhbTwo = pair.gear_two.FhBeta
    FhBeta5One = pair.gear_one.FhBeta5
    FhBeta5Two = pair.gear_two.FhBeta5
    fav = pair.fav
    kp = 0
    B1 = 0
    B2 = 0

    kv = __kv(pair, fmt)
    cGammaBeta = __c(pair)[1]

    if kv * fmt < 100:
        fm_b = 100
    else:
        fm_b = kv * fmt

    if helixModiffication == 1:
        B1 = 1.
        B2 = 1.
    if helixModiffication == 2:
        B1 = 1.
        B2 = 0.5
    if helixModiffication == 3:
        B1 = 0.1
        B2 = 1.
    if helixModiffication == 4:
        B1 = 0.1
        B2 = 0.5
    if helixModiffication == 5:
        B1 = 0.7
        B2 = 0.7

    stiff = d / shaftDiameter

    if stiff < 1.15:
        if schema == 1:
            kp = 0.8
        if schema == 2:
            kp = -0.8
        if schema == 3:
            kp = 1.33
        if schema == 4:
            kp = -0.6
        if schema == 5:
            kp = -1
    else:
        if schema == 1:
            kp = 0.48
        if schema == 2:
            kp = -0.48
        if schema == 3:
            kp = 1.33
        if schema == 4:
            kp = -0.36
        if schema == 5:
            kp = -0.6

    facTemp = (stiff ** 4 ) * (( l * s ) / ( d ** 2))
    fsh = fm_b * 0.023 * ( abs(1 + kp * facTemp - 0.3) + 0.3 ) * (( b / d ) ** 2 )
    fma = sqrt(FhbOne ** 2 + FhbTwo ** 2)

    if FhBeta5One > FhBeta5Two:
        fhb5 = FhBeta5One
    else:
        fhb5 = FhBeta5Two

    if fav == 1:
        fbX = abs(1.33 * B1 * fsh - fhb5)
    else:
        fbX = 1.33 * B1 * fsh + B2 * fma

    ybOne = __yb(pair.gear_one, pair, fbX)
    ybTwo = __yb(pair.gear_two, pair, fbX)
    yb = round(0.5 * ( ybOne + ybTwo), 1)
    fbY = fbX - yb

    if fbY * cGammaBeta / (2 * fm_b) >= 1:
        khb = sqrt(( 2. * fbY * cGammaBeta ) / fm_b)
    else:
        khb = 1 + ( fbY * cGammaBeta ) / ( 2 * fm_b )

    return khb


def kfb( pair, fmt ):
    khb = __khb(pair, fmt)
    bhOne = pair.gear_one.b / pair.gear_one.h
    bhTwo = pair.gear_two.b / pair.gear_two.h

    if bhOne < bhTwo:
        bh = bhOne
    else:
        bh = bhTwo
    if bh < 3:
        bh = 3

    nf = (bh ** 2) / ( 1 + bh + ( bh ** 2))

    return khb ** nf


def var( pair, fmt ):
    pair.gear_one.material.sh_limit = pair.gear_one.sigmaHLimit
    pair.gear_two.material.sh_limit = pair.gear_two.sigmaHLimit
    FfAlphaOne = pair.gear_one.FfAlpha
    FfAlphaTwo = pair.gear_two.FfAlpha
    materialOne = pair.gear_one.material
    materialTwo = pair.gear_two.material
    v = pair.v

    if FfAlphaOne > FfAlphaTwo:
        fpb = FfAlphaOne
    else:
        fpb = FfAlphaTwo

    khb = __khb(pair, fmt)
    kv = __kv(pair, fmt)
    cGammaAlpha = __c(pair)[0]

    yaOne = __ya(materialOne, pair.gear_one.material.sh_limit, v, FfAlphaOne, FfAlphaTwo)
    yaTwo = __ya(materialTwo, pair.gear_two.material.sh_limit, v, FfAlphaOne, FfAlphaTwo)
    ya = 0.5 * ( yaTwo + yaOne )

    if kv * fmt < 100:
        fm_b = 100
    else:
        fm_b = kv * fmt

    fthb = fm_b * khb

    return ( cGammaAlpha * ( fpb - ya ) ) / fthb


def kha( pair, fmt):
    epsilonAlpha = pair.epsilonAlpha
    epsilonGama = pair.epsilonGama
    epsilonBeta = pair.epsilonBeta

    var = __var(pair, fmt)

    if epsilonBeta is 0:
        zEpsilon = sqrt(( 4 - epsilonAlpha ) / 3)
    elif 1 > epsilonBeta > 0:
        zEpsilon = sqrt((( 4 - epsilonAlpha ) / 3 ) * ( 1 - epsilonBeta ) + ( epsilonBeta / epsilonAlpha ))
    else:
        zEpsilon = sqrt(1 / epsilonAlpha)

    if epsilonGama <= 2:
        kha = ( epsilonGama / 2 ) * ( 0.9 + 0.4 * var )
        if kha > epsilonGama / ( epsilonAlpha * zEpsilon ** 2 ):
            kha = epsilonGama / ( epsilonAlpha * zEpsilon ** 2 )
        elif kha < 1:
            kha = 1
    else:
        kha = 0.9 + 0.4 * var * sqrt(2 * ( epsilonGama - 1 ) / epsilonGama)
        if kha > epsilonGama / ( epsilonAlpha * zEpsilon ** 2 ):
            kha = epsilonGama / ( epsilonAlpha * zEpsilon ** 2 )
        elif kha < 1:
            kha = 1

    return kha


def kfa( pair, fmt ):
    epsilonAlpha = pair.epsilonAlpha
    epsilonGama = pair.epsilonGama

    var = __var(pair, fmt)

    if epsilonGama <= 2:
        kfa = ( epsilonGama / 2 ) * ( 0.9 + 0.4 * var )
        if kfa > epsilonGama / ( 0.25 * epsilonAlpha + 0.75 ):
            kfa = epsilonGama / (0.25 * epsilonAlpha + 0.75 )
        elif kfa < 1:
            kfa = 1
    else:
        kfa = 0.9 + 0.4 * var * sqrt(2 * ( epsilonGama - 1 ) / epsilonGama)
        if kfa > epsilonGama / ( 0.25 * epsilonAlpha + 0.75 ):
            kfa = epsilonGama / ( 0.25 * epsilonAlpha + 0.75 )
        elif kfa < 1:
            kfa = 1

    return kfa