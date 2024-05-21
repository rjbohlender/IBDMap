import gzip
from typing import List
import numpy as np


class GeneticMap:
    """Utility class for handling finding genetic distances.

    Args:
        gmaps: The collection of genetic map paths

    Attributes:
        gmap (Dict[int, Dict[int, float]]): The container for the recombination map.
    """

    def __init__(self, gmaps: List[str]):
        recom_map = {}
        for path in gmaps:
            gmap_file = gzip.open(path, 'rt')
            for line in gmap_file:
                line = line.strip().split()
                if line[0] == 'pos':
                    continue
                if line[1] == 'X' or line[1] == 'Y' or line[1] == 'MT':
                    break
                if line[1].startswith('chr'):  # Sanitize chromosome labels
                    line[1] = line[1][3:]
                if int(line[1]) not in recom_map:
                    recom_map[int(line[1])] = {}
                recom_map[int(line[1])][int(line[0])] = float(line[2])
        self.gmap = recom_map

    def find_nearest(self, chrom: int, pos: int) -> List[int]:
        """Find nearest element.

        Args:
            chrom: The chromosome number.
            pos: The position to search for.

        Returns:
            Pair of positions nearest the target.
        """
        gmap_list = np.array(list(self.gmap[chrom].keys()))
        nearest = np.abs(gmap_list - pos).argmin()
        if gmap_list[nearest] > pos:
            if nearest > 0:
                return [gmap_list[nearest - 1], gmap_list[nearest]]
            return [gmap_list[nearest], gmap_list[nearest + 1]]
        elif nearest == len(gmap_list) - 1:
            return [gmap_list[nearest - 1], gmap_list[nearest]]
        else:
            return [gmap_list[nearest], gmap_list[nearest + 1]]


