from collections.abc import MutableSequence

import math
import scipy.special

TURNS_FIRST_LAYER = 45
NUM_LAYERS = 26
WIRE_DIA = (1.0 / 1000.0)
WIRE_DIA_INSULATION = (1.09 / 1000.0)
LOOP_INNER_DIA = (76 / 1000.0 + WIRE_DIA_INSULATION)

WIRE_RADIUS = (WIRE_DIA / 2.0)
LOOP_INNER_RADIUS = (LOOP_INNER_DIA / 2.0)

LAYER_OFFSET = (WIRE_DIA_INSULATION * math.sqrt(3) / 2)

MU_0 = (1.25663706212 * 1.0E-6)

class Coil(MutableSequence):

    def __init__(self, loops):
        """init of loop.
        
        Parameters
        ----------
        loops : int
            Number of windings.
        """
        super().__init__()
        self._list = []
        
        for i in range(loops):
            self._list.append(self._make_loop())

    def __getitem__(self, index):
        """Return element(s) at index or slice."""
        return self._list[index]

    def __setitem__(self, index, value):
        """Set element(s) at index or slice."""
        self._list[index] = value

    def __delitem__(self, index):
        """Delete element(s) at index or slice."""
        del self._list[index]

    def __len__(self):
        """Return the length of self."""
        return len(self._list)

    def insert(self, index, value):
        """Insert one element at index."""
        self._list.insert(index, value)
        
    def _make_loop(self):
        """Add single loop to coil."""
        tmp = {
            "radius" : 0.0,
            "offset" : 0.0,
            "self_inductance" : 0.0,
            "mutual_inductance" : 0.0,
            }

        return tmp


def initialize_coil_geometry(loop):
    """Set up the position and radius of each conductor loop.
    
    Parameters
    ----------
    loop : Coil
        Container for loop specifications.
        ''
    Returns
    -------
    None.
    """
    
    # Initialize even-numbered layers 0, 2, 4, ...
    for i in range(0, NUM_LAYERS, 2):
        for j in range(0, TURNS_FIRST_LAYER):
            loop[int(j + i * TURNS_FIRST_LAYER - i / 2)]["radius"] = LOOP_INNER_RADIUS + i * LAYER_OFFSET
            loop[int(j + i * TURNS_FIRST_LAYER - i / 2)]["offset"] = - ((TURNS_FIRST_LAYER - 1) * WIRE_DIA_INSULATION) / 2 + j * WIRE_DIA_INSULATION
            
    # Initialize even-numbered layers 1, 3, 5, ...
    for i in range(1, NUM_LAYERS, 2):
        for j in range(0, TURNS_FIRST_LAYER - 1):
            loop[int(j + i * TURNS_FIRST_LAYER - int(i / 2))]["radius"] = (LOOP_INNER_RADIUS + LAYER_OFFSET) + (i - 1) * LAYER_OFFSET
            loop[int(j + i * TURNS_FIRST_LAYER - int(i / 2))]["offset"] = - (TURNS_FIRST_LAYER * WIRE_DIA_INSULATION) / 2 + (j + 1) * WIRE_DIA_INSULATION
            

def self_inductance(loop):
    """Calculate self inductance of each turn
    
    Parameters
    ----------
    loop : Coil
        Container for loop specifications.
        
    Returns
    -------
    None.
    """

    for i in range(0, loop_count):
        loop[i]["self_inductance"] = MU_0 * loop[i]["radius"] * (math.log(8 * loop[i]["radius"] / WIRE_RADIUS) - 2 + 0.25)


def mutual_inductance(loop):
    """Calculate mutual inductance between each pair of turns
    
    Parameters
    ----------
    loop : Coil
        Container for loop specifications.
        
    Returns
    -------
    None.
    """

    # Elliptic modulus 'k': argument to complete elliptic integrals K and E
    
    for i in range(0, loop_count):
        for j in range(0, loop_count):
            if i == j: continue
            
            center_dist = abs(loop[i]["offset"] - loop[j]["offset"])
            
            k = 2 * math.sqrt(loop[i]["radius"] * loop[j]["radius"] / ((loop[i]["radius"] + loop[j]["radius"]) * (loop[i]["radius"] + loop[j]["radius"]) + center_dist * center_dist))

            loop[i]["mutual_inductance"] += MU_0 * math.sqrt(loop[i]["radius"] * loop[j]["radius"]) * ((2 / k - k) * scipy.special.ellipk(k*k) - 2 / k * scipy.special.ellipe(k*k))


def total_inductance(loop):
    """Coil inductance = sum of self and mutual inductance
    
    Parameters
    ----------
    loop : Coil
        Container for loop specifications.
        
    Returns
    -------
    None.
    """

    si=0
    mi=0

    for i in range(0, loop_count):
        si += loop[i]["self_inductance"]

    for i in range(0, loop_count):
        mi += loop[i]["mutual_inductance"]

    print(f"Self inductance (sum):    {si} H")
    print(f"Mutual inductance (sum):   {mi} H")
    print(f"========================================")
    print(f"Total inductance:    {si + mi} H")


def get_loop_count():
    """Determine number of turns

    Returns
    -------
    loop_count : int
        Number of windings.
    """

    loop_count = 0

    for i in range (0,  NUM_LAYERS):
        if i % 2 == 0:
            loop_count += TURNS_FIRST_LAYER
        else:
            loop_count += TURNS_FIRST_LAYER - 1
    return loop_count


if __name__ == '__main__':

    loop_count = get_loop_count()

    loop = Coil(loop_count)

    initialize_coil_geometry(loop)
    
    self_inductance(loop)

    mutual_inductance(loop)
    
    total_inductance(loop)