__authors__ = ("Alexander M. Imre (@amimre)",)
__created__ = "2024-10-15"


class DisplacementsRange:
    _EPS = 1e-6

    def __init__(self, start, stop, step=None):
        self.start = start
        self.stop = stop
        self.step = step
        self.has_step = step is not None

    def __eq__(self, other):
        if not isinstance(other, DisplacementsRange):
            return False
        if self.has_step != other.has_step:
            return False
        if self.has_step:
            return (
                abs(self.start - other.start) < self._EPS
                and abs(self.stop - other.stop) < self._EPS
                and abs(self.step - other.step) < self._EPS
            )
        # otherwise, the step should be None
        if other.step is not None:
            return False
        return (
            abs(self.start - other.start) < self._EPS
            and abs(self.stop - other.stop) < self._EPS
        )

    def __repr__(self):
        if self.has_step:
            return f"(start={self.start}, stop={self.stop}, step={self.step})"
        return f"(start={self.start}, stop={self.stop})"
