class MantaToolsException(Exception):
    pass


class FieldNotFound(MantaToolsException):
    pass


class InfoFieldNotFound(MantaToolsException):
    pass


class GenotypeFieldNotFound(MantaToolsException):
    pass


class MissingMate(MantaToolsException):
    pass
