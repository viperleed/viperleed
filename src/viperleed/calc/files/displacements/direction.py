import numpy as np
import re


class Direction:
    """Class to parse and handle direction information in 3D space, supporting
    either Cartesian or fractional coordinates."""

    def __init__(self, direction_str):
        self._fractional = (
            "a" in direction_str or "b" in direction_str or "c" in direction_str
        )
        self.direction_str = direction_str
        self._vectors, self.num_free_directions = self._parse_direction(
            direction_str
        )

    def _parse_direction(self, direction_str):
        """Parse the direction string and return normalized 3D vectors and
        number of free directions."""
        dir_labels = "abc" if self._fractional else "xyz"

        if "[" in direction_str:  # Handle vector cases like 'xy[1 1]'
            vector_match = re.match(
                rf"^(?P<dir>[{dir_labels}]+)\[(?P<vec>[\d\s]+)\]$",
                direction_str,
            )
            if vector_match:
                directions = list(vector_match.group("dir"))
                vector_components = list(
                    map(float, vector_match.group("vec").split())
                )
                if len(directions) != len(vector_components):
                    raise ValueError(
                        "Mismatch between directions and vector components "
                        f"in {direction_str}"
                    )

                # Embed the 2D/1D vector into 3D space
                vector = self._embed_vector_in_3d(
                    directions, vector_components, dir_labels
                )
                vector = self._normalize_vectors([vector])
                return np.array(vector), 1  # Only 1 free direction
        else:
            # Handle simple cases like 'x', 'xy', 'abc'
            if not all(c in dir_labels for c in direction_str):
                raise ValueError(
                    "Mixing of fractional and Cartesian coordinates is not "
                    "allowed."
                )
            directions = list(direction_str)
            vectors = [
                self._get_basis_vector(d, dir_labels) for d in directions
            ]
            vectors = self._normalize_vectors(vectors)
            return np.array(vectors), len(vectors)

        raise ValueError(f"Invalid direction format: {direction_str}")

    def _get_basis_vector(self, direction, dir_labels):
        """Return the 3D basis vector or a placeholder for fractional
        coordinates 'a', 'b', 'c'."""
        basis_vectors = {
            "x": [1, 0, 0],
            "y": [0, 1, 0],
            "z": [0, 0, 1],
            "a": "a",
            "b": "b",
            "c": "c",
        }
        return basis_vectors[direction]

    def _embed_vector_in_3d(self, directions, vector_components, dir_labels):
        """Embed the 2D/1D vector into a 3D space."""
        embedded_vector = np.zeros(3)
        dir_map = {"x": 0, "y": 1, "z": 2, "a": 0, "b": 1, "c": 2}
        for direction, component in zip(directions, vector_components):
            embedded_vector[dir_map[direction]] = component
        return embedded_vector

    def _normalize_vectors(self, vectors):
        """Normalize vectors and raise error if zero-length vector is
        detected."""
        normalized_vectors = []
        for vec in vectors:
            if isinstance(vec, str):
                normalized_vectors.append(
                    vec
                )  # Skip normalizing for fractional placeholders
            else:
                norm = np.linalg.norm(vec)
                if norm == 0:
                    raise ValueError(f"Zero-length vector found: {vec}")
                normalized_vectors.append(np.array(vec) / norm)
        return normalized_vectors

    def get_frac_vectors(self, ucell):
        """Return the fractional vectors, scaled by the unit cell."""
        if not self._fractional:
            raise ValueError("No fractional coordinates in this direction")
        return self.get_cart_vectors(
            ucell
        )  # The fractional coordinates in this case correspond to the scaled ones

    def get_cart_vectors(self, ucell):
        """Scale fractional coordinates by the given unit cell matrix (3x3).

        ucell[:,0] is the x vector, ucell[:,1] is the y vector, and ucell[:,2]
        is the z vector.
        Only applies to 'a', 'b', and 'c' directions, returns cartesian
        coordinates."""
        scaled_vectors = []
        for vec in self._vectors:
            if isinstance(vec, str):
                # Map 'a', 'b', 'c' to corresponding unit cell vectors
                direction_map = {
                    "a": ucell[:, 0],
                    "b": ucell[:, 1],
                    "c": ucell[:, 2],
                }
                scaled_vectors.append(direction_map[vec])
            else:
                scaled_vectors.append(vec)
        # Normalize after scaling with the unit cell
        return np.array(
            [
                v / np.linalg.norm(v)
                for v in scaled_vectors
                if np.linalg.norm(v) > 0
            ]
        )

    def __eq__(self, other):
        """Compare two Direction objects for equality."""
        if not isinstance(other, Direction):
            return False
        return (
            np.allclose(self._vectors, other._vectors)
            and self._fractional == other._fractional
            and self.num_free_directions == other.num_free_directions
        )

    def __repr__(self):
        return (
            f"Direction(vectors={self._vectors}, "
            f"num_free_directions={self.num_free_directions}, "
            f"fractional={self._fractional})"
        )
