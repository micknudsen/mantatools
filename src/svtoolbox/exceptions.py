class SVToolBoxException(Exception):
    pass


class FieldNotFound(SVToolBoxException):
    pass


class InfoFieldNotFound(SVToolBoxException):
    pass


class GenotypeFieldNotFound(SVToolBoxException):
    pass


class MissingMate(SVToolBoxException):
    pass
