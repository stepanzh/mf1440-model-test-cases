#!/usr/bin/env python3

import constellation
import task

import pathlib
import matplotlib.colors
import matplotlib.pyplot as plt
import mpl_toolkits
import numpy as np


FIG_PARAMS = {'dpi': 300}


def plotSphere(ax, r, **kwargs):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, **kwargs)

    return ax


def plotEarth(ax, **kwargs):
    r = constellation.Const.earthRadius
    nx, ny, nz = 0, 0, r
    sx, sy, sz = 0, 0, -r
    
    ax.scatter([nx], [ny], [nz], color='blue', s=2)
    ax.text(nx, ny, nz, 'N', color='blue')
    
    ax.scatter([sx], [sy], [sz], color='red', s=2)
    ax.text(sx, sy, sz, 'S', color='red')
    
    return plotSphere(ax, r, alpha=0.25, color='gray')


def createPrettyFigure(fignum):
    ax = plt.figure(figsize=(10, 10), num=fignum).add_subplot(projection='3d')
    return ax


def prettifyAxes(ax):
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()

    return ax


def plotTrajectories(stateEci, satellites):
    ax = createPrettyFigure('trajectories')
    plotEarth(ax)

    for satIdx in satellites:
        ax.plot(
            stateEci[satIdx, 0, :],
            stateEci[satIdx, 1, :],
            stateEci[satIdx, 2, :],
        )

    prettifyAxes(ax)
    ax.set_title('trajectories')

    return ax


def plotPositions(ax, xyz, **kwargs):
    ax.scatter(
        xyz[:, 0],
        xyz[:, 1],
        xyz[:, 2],
        s=1,
        **kwargs,
    )
    return ax


def plotSatelliteProjections(cosmosPositions):

    ax = createPrettyFigure('projection-on-earth')
    ax = plotEarth(ax)
    ax = plotPositions(ax, cosmosPositions, label='position in cosmos')

    shadows = task.projectOnSphere(cosmosPositions, constellation.Const.earthRadius)
    ax = plotPositions(ax, shadows, color='red', label='\'shadow\' on Earth')

    prettifyAxes(ax)
    ax.set_title('positions of satellities and their \'shadows\' on Earth')

    return ax

def plotVoronoiTesselation(voronoi):
    ax = createPrettyFigure('voronoi')

    Pview = voronoi.points # P[P[:, 1] < -2e6]
    Vview = voronoi.vertices # V[V[:, 1] < -2e6]

    ax.computed_zorder = True 
    for (regionIdx, region) in enumerate(voronoi.regions):
        n = len(region)

        poly = mpl_toolkits.mplot3d.art3d.Poly3DCollection(
            [voronoi.vertices[region]],
            facecolors='gray',
            edgecolors='y',
            linewidths=0.5,
            alpha=1.0,
        )

        ax.add_collection(poly)

    ax.scatter(Pview[:, 0], Pview[:, 1], Pview[:, 2], s=0.5, color='k', label='satellite \'shadow\'')
    ax.scatter(Vview[:, 0], Vview[:, 1], Vview[:, 2], s=0.5, color='w', label='voronoi vertice')

    prettifyAxes(ax)
    ax.set_title('Satellite regions: Convex hull of Voronoi tesselation on Earth surface')

    return ax


def plotSatelliteRegions(voronoi, coloring):
    ax = createPrettyFigure('voronoi-neighbors')

    chromaticNumber = len(set(coloring.values()))

    assert len(matplotlib.colors.TABLEAU_COLORS.keys()) >= chromaticNumber
    regionColors = list(matplotlib.colors.TABLEAU_COLORS.keys())[0:chromaticNumber]

    for (regionIdx, region) in enumerate(voronoi.regions):
        n = len(region)

        ax.add_collection(
            mpl_toolkits.mplot3d.art3d.Poly3DCollection(
                [voronoi.vertices[region]],
                facecolors=regionColors[coloring[regionIdx]],
                alpha=1.0,
            ),
        )

    P = voronoi.points
    Pview = P[P[:, 1] < -4e6]

    ax.scatter(Pview[:, 0], Pview[:, 1], Pview[:, 2],
         s=1,
         color='k',
         label='satellite shadow (some hidden)',
    )
    ax.set_title('Satellite regions: Non-overlapped')

    ymin, ymax = ax.get_ylim()
    ylim = 1.15 * max(abs(ymin), abs(ymax))

    ax.set_xlim(-ylim, ylim)
    ax.set_ylim(-ylim, ylim)
    ax.set_zlim(-ylim, ylim)

    prettifyAxes(ax)

    return ax


def main():
    constellation_ = constellation.Constellation.createFromJson(pathlib.Path('../constellationsTest.json'), 'Starlink')
    epochs = list(range(1002))
    epochIdx = 0
    frequencyPool = [11, 12, 13, 14, 15, 16]

    stateEci = constellation_.propagateJ2(epochs)
    cosmosPositionsAtEpoch = stateEci[:, :, epochIdx]
    G = task.SatelliteShadowsGraph.createFromCosmosPositions(cosmosPositionsAtEpoch, earthRadius=constellation.Const.earthRadius)
    indexMap = task.SatelliteFrequenceIndexMap.tryCreate(G, len(frequencyPool))


    save = lambda ax, p: ax.figure.savefig(p, **FIG_PARAMS)

    ax = plotTrajectories(stateEci, range(0, stateEci.shape[0]))
    save(ax, pathlib.Path('figs/01-trajectories.png'))

    ax = plotSatelliteProjections(cosmosPositionsAtEpoch)
    save(ax, pathlib.Path('figs/02-projected.png'))

    ax = plotVoronoiTesselation(G.voronoi)
    save(ax, pathlib.Path('figs/03-voronoi.png'))

    ax = plotSatelliteRegions(G.voronoi, indexMap)
    save(ax, pathlib.Path('figs/04-voronoi-non-overlapped.png'))

    return None


if __name__ == '__main__':
    main()
