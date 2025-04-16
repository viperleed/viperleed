from collections.abc import MutableSequence

import math


TURNS_FIRST_LAYER = 45
NUM_LAYERS = 26
WIRE_DIA = (1.0 / 1000.0)
WIRE_DIA_INSULATION = (1.09 / 1000.0)

LOOP_INNER_DIA = (76 / 1000.0 + WIRE_DIA_INSULATION)
CENTER_DISTANCE = (0.5)
COIL_CURRENT = 2.5

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
            "flux_density" : 0.0,
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



def flux_density(loop):
    """Calculate flux density of each turn.
    
    Parameters
    ----------
    loop : Coil
        Container for loop specifications.
        ''
    Returns
    -------
    None.
    """

    for i in range (0, loop_count):
        loop[i]["flux_density"] = (MU_0 * COIL_CURRENT) / 2 * pow(loop[i]["radius"], 2) / (pow(pow(loop[i]["radius"], 2) + pow(CENTER_DISTANCE - loop[i]["offset"], 2), 1.5))

             


def total_flux_density(loop):
    """Coil flux density = superposition of each turn.
    
    Parameters
    ----------
    loop : Coil
        Container for loop specifications.
        ''
    Returns
    -------
    None.
    """
    
    B = 0.0
    for i in range (0, loop_count):
        B += loop[i]["flux_density"]
    print(f"Magnetic flux density B({CENTER_DISTANCE}): {B} T\n")



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

    flux_density(loop)
    
    total_flux_density(loop)
    
    