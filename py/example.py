#!/usr/bin/env python3

USAGE = """
usage:    ./example.py satellite-index epoch-index
example:  ./example.py 87 912
"""

from constellation import Constellation
import numpy as np
import pathlib
from random import randint
import sys


def main(argv):
    # TODO. Use argparse.
    try:
        satIdx = int(sys.argv[1])
        epochIdx = int(sys.argv[2])
    except ValueError:
        print('indexes must be ints', file=sys.stderr)
        print(USAGE)
        return None
    except IndexError:
        print('miss parameter(s)', file=sys.stderr)
        print(USAGE)
        return None

    constellation = Constellation.createFromJson(pathlib.Path('../constellationsTest.json'), 'Starlink')

    # определение точек на оси времени, в которые будут проихзводиться расчёты
    epochs = list(range(1002))

    stateEci = constellation.propagateJ2(epochs)

    print('Положение КА-' + str(satIdx) + ' на эпоху ' + str(epochs[epochIdx]) + ':')
    print(stateEci[satIdx, :, epochIdx])


if __name__ == '__main__':
    main(sys.argv)
