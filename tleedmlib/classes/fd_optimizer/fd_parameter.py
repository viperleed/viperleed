class FDParameter():
    """Base class for parameters accessible via full dynamic calculation."""

    def __init__(self, name, start_value, bounds, transform_func):
        self.name = name
        self.start_value = start_value
        # make sure start value is within bounds
        if start_value < bounds[0] or start_value > bounds[1]:
            raise ValueError(f"Start value for {start_value} is not within "
                             "bounds {bounds}")
        self.bounds = bounds

        self.transform_func = transform_func

    def apply(self, rparams, slab, val):
        """Apply the parameter to the given rparams and slab."""
        self.transform_func(rparams, slab, val)
