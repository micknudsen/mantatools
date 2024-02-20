class MantaToolsException(Exception):
    pass


class InfoFieldNotFound(MantaToolsException):
    pass


class GenotypeFieldNotFound(MantaToolsException):
    pass


class MissingMate(MantaToolsException):
    pass
