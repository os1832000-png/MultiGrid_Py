"""Custom warning classes emitted by MARFA."""


class MarfaWarning(UserWarning):
    """Base class for all MARFA warnings."""


class LineOutOfRangeWarning(MarfaWarning):
    """Requested spectral range contains no lines from the database."""


class InputFileWarning(MarfaWarning):
    """A HITRAN input file is missing or empty."""