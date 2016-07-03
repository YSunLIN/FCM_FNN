# encoding: utf-8
import math
import random
import copy


sqrt2 = math.sqrt(2)
sqrtPi = math.sqrt(math.pi)

def erf(x):
    return math.erf(x) * 0.5


def erfs2(ld, C, sigma):
    return erf(sqrt2 * (ld - C) / sigma)


def subsethood(inC, inSigma, outC, outSigma):
    if inC == outC:
        if inSigma > outSigma:
            return outSigma / inSigma
        else:
            return inSigma / outSigma
    elif inC > outC:
        tmp = inC; inC = outC; outC = tmp
        tmp = inSigma; inSigma = outSigma; outSigma = tmp

    # avoid outSigma - inSigma is 0
    if outSigma != inSigma:
        ld1 = (outSigma * inC - inSigma * outC) / (outSigma - inSigma)
        erfs2_ld1_i = erfs2(ld1, inC, inSigma)
        erfs2_ld1_o = erfs2(ld1, outC, outSigma)

    ld2 = (outSigma * inC + inSigma * outC) / (outSigma + inSigma)
    erfs2_ld2_i = erfs2(ld2, inC, inSigma)
    erfs2_ld2_o = erfs2(ld2, outC, outSigma)

    if inSigma == outSigma:
        numerator = inSigma * (0.5 - erfs2_ld2_i) + outSigma * (0.5 + erfs2_ld2_o)
        denominator = inSigma * (0.5 + erfs2_ld2_i) + outSigma * (0.5 - erfs2_ld2_o)
    elif inSigma < outSigma:
        numerator = outSigma * (erfs2_ld2_o - erfs2_ld1_o) + inSigma * (1 + erfs2_ld1_i - erfs2_ld2_i)
        denominator = outSigma * (1 + erfs2_ld1_o - erfs2_ld2_o) + inSigma * (erfs2_ld2_i - erfs2_ld1_i)
    else:
        numerator = outSigma * (1 + erfs2_ld2_o - erfs2_ld1_o) + inSigma * (erfs2_ld1_i - erfs2_ld2_i)
        denominator = outSigma * (erfs2_ld1_o - erfs2_ld2_o) + inSigma * (1 + erfs2_ld2_i - erfs2_ld1_i)

    return numerator / denominator


def card(inC, inSigma, outC, outSigma):
    if inC == outC:
        if inSigma < outSigma:
            return inSigma * sqrtPi
        else:
            return outSigma * sqrtPi
    elif inC > outC:
        tmp = inC; inC = outC; outC = tmp
        tmp = inSigma; inSigma = outSigma; outSigma = tmp
    
    # avoid outSigma - inSigma is 0
    if outSigma != inSigma:
        ld1 = (outSigma * inC - inSigma * outC) / (outSigma - inSigma)
        erfs2_ld1_i = erfs2(ld1, inC, inSigma)
        erfs2_ld1_o = erfs2(ld1, outC, outSigma)

    ld2 = (outSigma * inC + inSigma * outC) / (outSigma + inSigma)
    erfs2_ld2_i = erfs2(ld2, inC, inSigma)
    erfs2_ld2_o = erfs2(ld2, outC, outSigma)

    if inSigma == outSigma:
        return outSigma * sqrtPi * (erfs2_ld2_o + 0.5) + inSigma * sqrtPi * (0.5 - erfs2_ld2_i)
    elif inSigma < outSigma:
        return inSigma * sqrtPi * (1 + erfs2_ld1_i - erfs2_ld2_i) + outSigma * sqrtPi * (erfs2_ld2_o - erfs2_ld1_o)
    else:
        return outSigma * sqrtPi * (1 + erfs2_ld2_o - erfs2_ld1_o) + inSigma * sqrtPi * (erfs2_ld1_i - erfs2_ld2_i)


def deltaCardC(inC, inSigma, outC, outSigma):
    if inC == outC:
        return 0

    # avoid outSigma - inSigma is 0
    if outSigma != inSigma:
        ld1 = (outSigma * inC - inSigma * outC) / (outSigma - inSigma)
        ld1_i = (ld1 - inC) / inSigma
        ld1_o = (ld1 - outC) / outSigma
        exp_ld1_i = math.exp(-ld1_i**2)
        exp_ld1_o = math.exp(-ld1_o**2)

    ld2 = (outSigma * inC + inSigma * outC) / (outSigma + inSigma)
    ld2_i = (ld2 - inC) / inSigma
    ld2_o = (ld2 - outC) / outSigma
    exp_ld2_i = math.exp(-ld2_i**2)
    exp_ld2_o = math.exp(-ld2_o**2)

    if inC < outC:
        if inSigma == outSigma:
            return exp_ld2_i
        elif inSigma < outSigma:
            return exp_ld2_i - exp_ld1_i
        else:
            return exp_ld2_i - exp_ld1_i
    else:
        if inSigma == outSigma:
            return -exp_ld2_i
        elif inSigma < outSigma:
            return exp_ld1_i - exp_ld2_i
        else:
            return exp_ld1_i - exp_ld2_i


def deltaCardSigma(inC, inSigma, outC, outSigma):
    if inC == outC:
        if inSigma == outSigma:
            return 0
        elif inSigma < outSigma:
            return 1.0 / outSigma
        else:
            return -outSigma / inSigma ** 2

    # avoid outSigma - inSigma is 0
    if outSigma != inSigma:
        ld1 = (outSigma * inC - inSigma * outC) / (outSigma - inSigma)
        erfs2_ld1_i = erfs2(ld1, inC, inSigma)
        ld1_i = (ld1 - inC) / inSigma
        exp_ld1_i = math.exp(-ld1_i**2)

    ld2 = (outSigma * inC + inSigma * outC) / (outSigma + inSigma)
    erfs2_ld2_i = erfs2(ld2, inC, inSigma)
    ld2_i = (ld2 - inC) / inSigma
    exp_ld2_i = math.exp(-ld2_i**2)

    if inC < outC:
        if inSigma == outSigma:
            return ld2_i * exp_ld2_i + sqrtPi * (0.5 - erfs2_ld2_i)
        elif inSigma < outSigma:
            return ld2_i * exp_ld2_i - ld1_i * exp_ld1_i + sqrtPi * (1 + erfs2_ld1_i - erfs2_ld2_i)
        else:
            return ld2_i * exp_ld2_i - ld1_i * exp_ld1_i + sqrtPi * (erfs2_ld1_i - erfs2_ld2_i)
    else:
        if inSigma == outSigma:
            return sqrtPi * (0.5 + erfs2_ld2_i) - ld2_i * exp_ld2_i
        elif inSigma < outSigma:
            return ld1_i * exp_ld1_i - ld2_i * exp_ld2_i + sqrtPi * (1 + erfs2_ld2_i - erfs2_ld1_i)
        else:
            return ld1_i * exp_ld1_i - ld2_i * exp_ld2_i + sqrtPi * (erfs2_ld2_i - erfs2_ld1_i)


def deltaSubsetC(inC, inSigma, outC, outSigma):
    _deltaCardC = deltaCardC(inC, inSigma, outC, outSigma)
    numerator = card(inC, inSigma, outC, outSigma)
    denominator = sqrtPi * (outSigma + inSigma) - numerator
    return ( _deltaCardC * denominator - numerator * (-_deltaCardC) ) / denominator ** 2


def deltaSubsetSigma(inC, inSigma, outC, outSigma):
    _deltaCardSigma = deltaCardSigma(inC, inSigma, outC, outSigma)
    numerator = card(inC, inSigma, outC, outSigma)
    denominator = sqrtPi * (outSigma + inSigma) - numerator
    return ( _deltaCardSigma * denominator - numerator * (sqrtPi - _deltaCardSigma) ) / denominator ** 2


class Concept(object):
    def __init__(self, numOfTerms):
        self.numOfTerms = numOfTerms
        self.C = []
        self.sigma = []
        self.xi = []
        for i in range(numOfTerms):
            self.C.append(random.uniform(i * 2.0  / numOfTerms, (i + 1) * 2.0  / numOfTerms))
            self.sigma.append(random.uniform(0, 1))  # sigma > 0
            self.xi.append(random.uniform(0, 1))


class FCM_FNN(object):
    def __init__(self, *numOfTermsList):
        self.concepts = []
        for numOfTerms in numOfTermsList:
            self.concepts.append(Concept(numOfTerms))
        # 把下标0的层忽略掉
        self.layerf = [[] for i in range(5)]
        self.layerx = [[] for i in range(5)]
        # layer1
        self.layerf[1] = [0 for i in range(len(numOfTermsList))]
        self.layerx[1] = [0 for i in range(len(numOfTermsList))]
        # layer2
        self.layerf[2] = [[] for i in range(len(numOfTermsList))]
        self.layerx[2] = [[] for i in range(len(numOfTermsList))]
        for i in range(len(numOfTermsList)):
            for count in range(numOfTermsList[i]):
                self.layerf[2][i].append(0)
                self.layerx[2][i].append(0)
        # layer3
        self.layerf[3] = [[] for i in range(len(numOfTermsList))]
        self.layerx[3] = [[] for i in range(len(numOfTermsList))]
        for oj in range(len(numOfTermsList)):
            for omj in range(numOfTermsList[oj]):
                self.layerx[3][oj].append(0)
                self.layerf[3][oj].append([])
                for ii in range(len(numOfTermsList)):
                    self.layerf[3][oj][omj].append([])
                    for ini in range(numOfTermsList[ii]):
                        self.layerf[3][oj][omj][ii].append(0)
        # layer4
        self.layerf[4] = [0 for i in range(len(numOfTermsList))]
        self.layerx[4] = [0 for i in range(len(numOfTermsList))]


    def predict(self, X):
        if len(X) != len(self.concepts):
            raise Exception("The dimension of X is inavailable")
        # layer1
        for i, v in enumerate(X):
            self.layerf[1][i] = v
            self.layerx[1][i] = self.layerf[1][i]

        # layer2
        for ii, concept in enumerate(self.concepts):
            for ini in range(concept.numOfTerms):
                self.layerf[2][ii][ini] = -((self.layerx[1][ii] - concept.C[ini]) / concept.sigma[ini]) ** 2
                self.layerx[2][ii][ini] = math.exp(self.layerf[2][ii][ini])

        # layer3
        # f
        for oj, oconcept in enumerate(self.concepts):
            for omj in range(oconcept.numOfTerms):
                for ii, iconcept in enumerate(self.concepts):
                    for ini in range(iconcept.numOfTerms):
                        self.layerf[3][oj][omj][ii][ini] = self.layerx[2][ii][ini] * \
                            (1 - subsethood(iconcept.C[ini], iconcept.sigma[ini], oconcept.C[omj], oconcept.sigma[omj]))
        # x
        for oj, oconcept in enumerate(self.concepts):
            for omj in range(oconcept.numOfTerms):
                numerator = 0
                denominator = 0
                for ii, iconcept in enumerate(self.concepts):
                    if ii != oj:
                        for ini in range(iconcept.numOfTerms):
                            numerator += self.layerf[3][oj][omj][ii][ini] * iconcept.C[ini] * iconcept.sigma[ini]
                            denominator += self.layerf[3][oj][omj][ii][ini] * iconcept.sigma[ini]
                self.layerx[3][oj][omj] = numerator / denominator

        # layer4
        for oj, concept in enumerate(self.concepts):
            self.layerf[4][oj] = 0
            for omj in range(concept.numOfTerms):
                self.layerf[4][oj] += concept.xi[omj] * self.layerx[3][oj][omj]
            self.layerx[4][oj] = self.layerf[4][oj]

        return copy.copy(self.layerx[4])


    def test(self, dataSet):
        err = 0
        for X, D in dataSet:
            Y = self.predict(X)
            tmpErr = 0
            for j in range(len(Y)):
                tmpErr += 0.5 * (D[j] - Y[j]) ** 2
            err += tmpErr
        return math.sqrt(err * 2 / (len(dataSet) * len(dataSet[0][1])))


    def train(self, dataSet, learnRate=0.01):
        err = 0
        for X, D in dataSet:
            Y = self.predict(X)

            tmpErr = 0
            for j in range(len(Y)):
                tmpErr += 0.5 * (D[j] - Y[j]) ** 2
            err += tmpErr

            for ii, iconcept in enumerate(self.concepts):
                for ini in range(iconcept.numOfTerms):
                    # xi
                    deltaXi = - (D[ii] - Y[ii]) * self.layerx[3][ii][ini]
                    iconcept.xi[ini] -= learnRate * deltaXi

                    # C, sigma
                    deltaC = 0
                    deltaSigma = 0
                    for l, lconcept in enumerate(self.concepts):
                        if l == ii: continue
                        mlSumC = 0
                        mlSumSigma = 0
                        for ml in range(lconcept.numOfTerms):
                            numerator = denominator = 0
                            for ti, tconcept in enumerate(self.concepts):
                                if ti != l:
                                    for tni in range(tconcept.numOfTerms):
                                        numerator += self.layerf[3][l][ml][ti][tni] * tconcept.C[tni] * tconcept.sigma[tni]
                                        denominator += self.layerf[3][l][ml][ti][tni] * tconcept.sigma[tni]
                            delta_yl = lconcept.xi[ml]
                            delta_x = (iconcept.C[ini] * iconcept.sigma[ini] * denominator - iconcept.sigma[ini] * numerator) / denominator ** 2
                            delta_f_C = self.layerx[2][ii][ini] * (2 * (self.layerx[1][ii] - iconcept.C[ini]) / iconcept.sigma[ini] ** 2) \
                                * (1 - subsethood(iconcept.C[ini], iconcept.sigma[ini], lconcept.C[ml], lconcept.sigma[ml])) \
                                - self.layerx[2][ii][ini] * deltaSubsetC(iconcept.C[ini], iconcept.sigma[ini], lconcept.C[ml], lconcept.sigma[ml])
                            delta_f_Sigma = self.layerx[2][ii][ini] * (2 * (self.layerx[1][ii] - iconcept.C[ini]) ** 2 / iconcept.sigma[ini] ** 3) \
                                * (1 - subsethood(iconcept.C[ini], iconcept.sigma[ini], lconcept.C[ml], lconcept.sigma[ml])) \
                                - self.layerx[2][ii][ini] * deltaSubsetSigma(iconcept.C[ini], iconcept.sigma[ini], lconcept.C[ml], lconcept.sigma[ml])

                            mlSumC += delta_yl * delta_x * delta_f_C
                            mlSumSigma += delta_yl * delta_x * delta_f_Sigma

                        deltaC += mlSumC * (Y[l] - D[l])
                        deltaSigma += mlSumSigma * (Y[l] - D[l])

                    iconcept.C[ini] -= learnRate * deltaC
                    iconcept.sigma[ini] -= learnRate * deltaSigma
                    if iconcept.sigma[ini] < 0:
                        iconcept.sigma[ini] = 1e-100
                        print "Warn: maybe it's impossible to train with this DataSet and learning rate"
        return math.sqrt(err * 2 / (len(dataSet) * len(dataSet[0][1])))
