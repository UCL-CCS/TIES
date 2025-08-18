"""
Updated OOP approach to data analysis.
"""


class Replica:
    """
    Replica loads in and parses the actual information.
    Multiple replicas can work on the lambda.
    So replica is defined by lambda and by its directory path.
    This representation should contain all the details necessary.
    """

    pass


class Lambda:
    """
    Reflect the real lambda rather than the "general" lambda.
    However, stores the information about the "general" lambda as well.
    This class contains at least 1 replica.
    """

    pass


class Contribution:
    """
    Reflects one type of interactions. For example, appearing electrostatics, or dissapearing VDW.
    This class contains lambdas with their replicas.
    It can calculate the integral and plot different information relevant to each contribution.
    """

    pass


class DGSystem:
    """
    Single step dG system.

    Contains 4 contributions: disappearing and appearing electrostatics and vdw
    Uses contributions to calculate dG.
    Contains lots of dG analysis and plotting.
    """

    pass


class TCSystem:
    """
    Thermodynamics Cycle System.

    Contains 2 Systems, each providing one dG.
    This way it can provide the ddG.
    Contains lots of ddG analysis and plotting.
    """

    pass
