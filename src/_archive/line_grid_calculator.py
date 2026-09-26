# ==============================================================================
# LINE GRID CALCULATOR  (direct port of Fortran LineGridCalc)
# ==============================================================================

class LineGridCalc:
    """
    Python port of the Fortran LineGridCalc module.

    Implements leftLBL_full, centerLBL_full, rightLBL_full which route each
    spectral line's shape-function sample to the correct multi-resolution
    stencil cell, exactly as the Fortran GOTO-cascade does.

    The FSHAPE callable plays the role of Fortran's procedure(shape) pointer:
    it receives a single float offset and returns a float lineshape value.
    """

    def __init__(self, gs: 'GridState'):
        self.gs = gs
