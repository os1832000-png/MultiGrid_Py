# ==============================================================================
# ATMOSPHERIC PROFILE
# ==============================================================================

@dataclass
class AtmosphericProfile:
    z:   np.ndarray
    P:   np.ndarray
    T:   np.ndarray
    vmr: Dict[int, np.ndarray]

    @classmethod
    def us_standard(cls):
        z   = np.array([0, 10, 20, 30, 50, 70, 100], dtype=float)
        P   = np.array([1.0, 0.265, 0.055, 0.012, 0.001, 0.00005, 0.000003])
        T   = np.array([288, 223, 217, 227, 271, 220, 210], dtype=float)
        vmr = {2: 400e-6 * np.ones_like(z)}
        return cls(z, P, T, vmr)

