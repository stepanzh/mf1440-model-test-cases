#!/usr/bin/env python3

"""
(1) Реализуйте вычислительную процедуру, позволяющую назначать разные частоты
    соседним КА в группировке (для уменьшения интерференции сигнала). Соседними
    считаются КА, зоны обслуживания которых имеют общую границу.

- Зоны обслуживания на поверхности Земли получаются с помощью разбиения Вороного для подспутниковых точек КА в каждый момент времени.

(2) Предложите способ, позволяющий рационально минимизировать количество
    переключений частот на КА системы.

- Частоты на все КА назначаются из некоторого наперед заданного набора несущих частот.

(3) При каком минимальном размере этого набора (частот), переключения частоты на любом КА группировки будут происходить не чаще чем один раз за виток.

"""

"""
(1) 1. Choose some epoch
    2. Load coords in cosmos
    3. Project coords on Earth surface
    4. Voronoi on the surface => neighbors => graph
    5. Graph coloring with zones
    6. Map zones to radio frequencies
"""

import abc
import constellation
import gcol
import logging
import networkx
import numpy as np
import pathlib
import scipy


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

logger = logging.getLogger(__name__)


def projectOnSphere(P, r):
    scale = r / np.linalg.norm(P, axis=1)
    Q = P * scale[:, np.newaxis]
    return Q


class SatelliteShadowsGraph(networkx.Graph):
    def __init__(self, voronoi: scipy.spatial.SphericalVoronoi, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._voronoi = voronoi
        logger.debug('DONE: graph of neighbors')

    @property
    def voronoi(self):
        return self._voronoi

    @property
    def satellitePositions(self):
        return self._voronoi.points
        
    @classmethod
    def createFromCosmosPositions(cls, positionsAtEpoch, earthRadius = constellation.Const.earthRadius):
        P = projectOnSphere(positionsAtEpoch, earthRadius)
        origin = np.array([0, 0, 0])
        sv = scipy.spatial.SphericalVoronoi(P, earthRadius, origin)
        sv.sort_vertices_of_regions()
        edges = __class__.findNeighborsAsEdges(sv)
        return cls(sv, edges)

    @staticmethod
    def findNeighborsAsEdges(sv):
        # Two satellites are neighbors, when they share a pair of Voronoi points
    
        edges = []
        for i1 in range(len(sv.points)):
            iSatVoronoi = sv.regions[i1]
            targetNeighborsCount = len(iSatVoronoi)
            count = 0
            
            for i2 in range(i1 + 1, len(sv.points)):
                common = set(sv.regions[i1]).intersection(sv.regions[i2]) # TODO. Can be faster?
                assert len(common) <= 2
                assert len(common) != 1
    
                if len(common) == 0:
                    continue
    
                # i1 and i2 are neighbors!
                edges.append([i1, i2])
                count += 1
                # All neighbors of i1 point are found
                if count == targetNeighborsCount:
                    break
    
        return edges


class NodeKColorer(abc.ABC):
    @abc.abstractmethod
    def color(self, graph, k) -> dict:
        pass


class GcolKNodeColorer(NodeKColorer):
    def color(self, graph, k):
        coloring = gcol.node_k_coloring(graph, k)
        logger.debug(f'DONE: coloring')
        return coloring


class SatelliteFrequenceIndexMap(dict):
    @property
    def pool(self):
        # TODO. Should be memoized.
        return set(self.values())
    
    @classmethod
    def tryCreate(cls, G: SatelliteShadowsGraph, frequencePoolSize: int, colorer: NodeKColorer = GcolKNodeColorer()):
        coloring = colorer.color(G, frequencePoolSize)
        return cls(coloring)


class SatelliteFrequenceMap(dict):
    def __init__(self, indexMap: SatelliteFrequenceIndexMap, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._indexMap = indexMap
        logger.debug('DONE: frequency map')
        
    @property
    def pool(self):
        # TODO. Should be memoized.
        return set(self.values())

    @property
    def indexMap(self):
        return self._indexMap
    
    @classmethod
    def createFromFrequencyList(cls, idxMap: SatelliteFrequenceIndexMap, frequencyList: list):
        assert len(idxMap.pool) == len(frequencyList)
        trans = dict(enumerate(frequencyList))
        target = { k: trans[idxMap[k]] for k in idxMap.keys() }
        return cls(idxMap, target)


if __name__ == '__main__':
    constellation_ = constellation.Constellation.createFromJson(pathlib.Path('../constellationsTest.json'), 'Starlink')
    epochs = list(range(1002))
    stateEci = constellation_.propagateJ2(epochs)
    satIdx, epochIdx = 87, 912

    satellitePositions = stateEci[:, :, epochIdx]
    logger.debug(f'satellite position s/t/x/y/z = {satIdx} / {epochIdx} / {stateEci[satIdx, :, epochIdx]}')

    G = SatelliteShadowsGraph.createFromCosmosPositions(satellitePositions, earthRadius=constellation.Const.earthRadius)

    frequencyPool = [11, 12, 13, 14, 15, 16]
    freqMap = SatelliteFrequenceMap.createFromFrequencyList(
        SatelliteFrequenceIndexMap.tryCreate(G, len(frequencyPool)),
        frequencyPool
    )

    satelliteIndexes = list(sorted(freqMap.keys()))

    logging.info('frequencies of first 20 satellites')
    print('satIdx freqIdx  freq')
    for satIdx in satelliteIndexes[:20]:
        frequency = freqMap[satIdx]
        freqIdx = freqMap.indexMap[satIdx]
        print(f'{satIdx:6d}  {freqIdx:6d}  {frequency:6f}')
