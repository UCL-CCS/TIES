
import ties.generator


class Protein:
    """
    A helper tool for the protein file. It calculates the number of ions needed to neutralise it
    (using ambertools for now).

    :param filename: filepath to the protein
    :type filename: string
    :param config: Optional configuration for the protein
    :type config: :class:`Config`
    """
    def __init__(self, filename, config):
        self.file = filename
        self.config = config

        # fixme - check if the file exists at this point, throw an exception otherwise

        # calculate the charges of the protein (using ambertools)
        # fixme - turn this into a method? stage2: use propka or some other tool, not this workaround
        self.protein_net_charge = ties.generator.get_protein_net_charge(config.workdir, config.protein.absolute(),
                                                                   config.ambertools_tleap, config.tleap_check_protein,
                                                                   config.protein_ff)

        print(f'Protein net charge: {self.protein_net_charge}')

    def get_path(self):
        """
        Get a path to the protein.

        :return: the protein filename
        :rtype: string
        """
        return self.file