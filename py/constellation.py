import numpy as np
import json
from typing import NamedTuple


# TODO. Why not dataclass?
# FIX. Singleton.
class Parameters(object):
    pass


# TODO. It's better to keep unit in variable name or use unitful values (e.g. like in Unitful.jl).
Const = Parameters()
Const.earthRadius = 6378135      # Экваториальный радиус Земли [m]
Const.earthGM = 3.986004415e+14  # Гравитационный параметр Земли [m3/s2]
Const.earthJ2 = 1.082626e-3      # Вторая зональная гармоника геопотенциала


# TODO. Why not dataclass?
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


# TODO. Should be dataclass + builder from groups?
class Constellation:
    def __init__(self, groups):
        self.groups = groups
        self.totalSatCount = Constellation.countSatellitesTotal(groups)
        self.elements = Constellation._initializeOrbitElements(groups, self.totalSatCount)

    @classmethod
    def countSatellitesTotal(cls, groups):
        return sum(map(lambda g: g.getTotalSatCount(), groups))

    # NOTE. Consider to use loader class (DI) and/or constructor from dictionary.
    @classmethod
    def createFromJson(cls, path, nameCode: str):
        with open(path) as io:
            jsonData = json.load(io)

        # Find first occurrence of constellation with nameCode
        data = filter(lambda x: x['name'].lower() == nameCode.lower(), jsonData)
        constellationData = next(data, None)
        if constellationData is None:
            raise Exception('Группировка не найдена в файле')

        print('Загружена группировка ' + nameCode)

        groups = list(map(lambda w: WalkerGroup(*w), constellationData['Walkers']))

        return Constellation(groups)

    @classmethod
    def _initializeOrbitElements(cls, groups, totalSatCount = None):
        "вычисление элементов орбиты для всех КА в начальный момент"

        if totalSatCount is None:
            totalSatCount = cls.countSatellitesTotal(groups)

        elements = np.zeros((totalSatCount, 6))
        shift = 0

        for group in groups:
            ending = shift + group.getTotalSatCount()
            elements[shift:ending, :] = group.getInitialElements()
            shift = ending

        return elements

    def propagateJ2(self, epochs, constants: Parameters = Const):
        "расчёт положений всех КА в заданные моменты времени"

        stateEci = np.zeros((self.totalSatCount, 3, len(epochs)))

        inclination = self.elements[:, 4]
        sma = self.elements[:, 0]
        raan0 = self.elements[:, 3]
        aol0 = self.elements[:, 5]

        raanPrecessionRate = (
                -1.5 * (
                    constants.earthJ2 * np.sqrt(constants.earthGM) * constants.earthRadius**2
                ) / (sma**(7/2))
                * np.cos(inclination)
        )

        draconicOmega = (
            np.sqrt(constants.earthGM / sma**3)
            * (1 - 1.5 * constants.earthJ2 * (constants.earthRadius / sma)**2)
            * (1 - 4 * np.cos(inclination)**2)
        )

        for epoch in epochs:
            aol = aol0 + epoch * draconicOmega
            raanOmega = raan0 + epoch * raanPrecessionRate

            epochState = sma * [
                (np.cos(aol) * np.cos(raanOmega) - np.sin(aol) * np.cos(inclination) * np.sin(raanOmega)),
                (np.cos(aol) * np.sin(raanOmega) + np.sin(aol) * np.cos(inclination) * np.cos(raanOmega)),
                (np.sin(aol) * np.sin(inclination))
            ]

            stateEci[:, :, epochs.index(epoch)] = np.array(epochState).T

        return stateEci
