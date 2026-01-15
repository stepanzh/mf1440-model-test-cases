import numpy as np
import json
from typing import NamedTuple


class Parameters(object):
    pass


Const = Parameters()
Const.earthRadius = 6378135      # Экваториальный радиус Земли [m]
Const.earthGM = 3.986004415e+14  # Гравитационный параметр Земли [m3/s2]
Const.earthJ2 = 1.082626e-3      # Вторая зональная гармоника геопотенциала

group = Parameters()


class Walker(NamedTuple):
    inclination: float           # наклонение орбиты
    satsPerPlane: int            # число КА в каждой орбитальной плоскости группы
    planeCount: int              # число орбитальных плоскостей в группе
    f: int                       # фазовый сдвиг по аргументу широты между КА в соседних плоскостях
    altitude: float              # высота орбиты
    maxRaan: float               # максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей)
    startRaan: float             # прямое восхождение восходящего узла для первой плоскости


class WalkerGroup(Walker):

    def getTotalSatCount(self):
        return self.satsPerPlane * self.planeCount

    def getInitialElements(self):
        startRaan   = np.deg2rad(self.startRaan)
        maxRaan     = np.deg2rad(self.maxRaan)
        inclination = np.deg2rad(self.inclination)
        altitude    = self.altitude * 1000
        satCount    = self.getTotalSatCount()

        raans = np.linspace(startRaan, startRaan + maxRaan, self.planeCount + 1)
        raans = raans[:-1] % (2 * np.pi)

        elements = np.zeros((satCount, 6))
        idx = 0

        for raanIdx, raan in enumerate(raans):
            for satIdx in range(self.satsPerPlane):
                sma = Const.earthRadius + altitude
                aol = 2 * np.pi * (satIdx / self.satsPerPlane + self.f * raanIdx / satCount)

                elements[idx, :] = [sma, 0, 0, raan, inclination, aol]
                idx += 1

        return elements


class Constellation:

    # TODO. Deprecate name code -> pass to data to constructor.
    # like (self, count, groups, elements, state_eci)

    def __init__(self, nameCode):
        # TODO. Why below are lists? Why not numpy arrays?

        self.groups = None
        self.totalSatCount = 0  # NOTE. Should be property based on self.groups

        self.elements = None
        self.stateEci = None

        self.loadFromConfig(nameCode)

    # TODO. Not pure. Better to have loader of Constellation from Json.
    # like
    #   constellation = ConstellationJsonLoader(data_dir).create_from_code(code)

    def loadFromConfig(self, nameCode: str):
        assert self.groups is None, 'group of satellites was already loaded'

        with open('../ConstellationsTest.json') as io:
            jsonData = json.load(io)

        # Find first occurrence of constellation with nameCode
        data = filter(lambda x: x['name'].lower() == nameCode.lower(), jsonData)
        constellationData = next(data, None)
        if constellationData is None:
            raise Exception('Группировка не найдена в файле')

        print('Загружена группировка ' + nameCode)

        self.groups = list(map(lambda w: WalkerGroup(*w), constellationData['Walkers']))
        self.totalSatCount = sum(map(lambda g: g.getTotalSatCount(), self.groups))

        return None

    def getInitialState(self):
        assert self.elements is None, 'use loadFromConfig first'

        self.elements = np.zeros((self.totalSatCount, 6))
        shift = 0

        for singleGroup in self.groups:
            ending = shift + singleGroup.getTotalSatCount()
            self.elements[shift:ending, :] = singleGroup.getInitialElements()
            shift = ending

    def propagateJ2(self, epochs):
        self.stateEci = np.zeros((self.totalSatCount, 3, len(epochs)))

        inclination = self.elements[:, 4]
        sma = self.elements[:, 0]
        raan0 = self.elements[:, 3]
        aol0 = self.elements[:, 5]

        raanPrecessionRate = (
                -1.5 * (
                    Const.earthJ2 * np.sqrt(Const.earthGM) * Const.earthRadius**2
                ) / (sma**(7/2))
                * np.cos(inclination)
        )

        draconicOmega = (
            np.sqrt(Const.earthGM / sma**3)
            * (1 - 1.5 * Const.earthJ2 * (Const.earthRadius / sma)**2)
            * (1 - 4 * np.cos(inclination)**2)
        )

        for epoch in epochs:
            aol = aol0 + epoch * draconicOmega
            raanOmega = raan0 + epoch * raanPrecessionRate

            epochState = sma * [
                (np.cos(aol) * np.cos(raanOmega) - np.sin(aol) * np.cos(inclination) * np.sin(raanOmega)),
                (np.cos(aol) * np.sin(raanOmega) + np.sin(aol) * np.cos(inclination) * np.cos(raanOmega)),
                (np.sin(aol) * np.sin(inclination))]

            self.stateEci[:, :, epochs.index(epoch)] = np.array(epochState).T
