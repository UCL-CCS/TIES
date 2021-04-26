


class Hybrid:
    """
    The class for the user interface that hides the complexity of SuperimposedTopology.

    Functionalities:
     - generate the metadata (.json) from the suptop
     - load itself from the metadata
     - allow for the analysis once the results are available
    """

    def __init__(self, config, morph, suptop):
        self.config = config
        self.morph = morph
        self.suptop = suptop
