#!/usr/bin/env python3
"""MarfaPython: simplified 3-grid MARFA implementation.

CLI mirrors `marfa_simple`:
    marfa_python <parfile> <maxlines> <molec_id> <wv_min> <wv_max> <cutoff>
                 <temp> <p_self> <p_foreign> <TIPS_ref> <TIPS_target> <outfile>
"""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass

import numpy as np
from numba import njit


PI = np.float64(3.1415926)
SQLN2 = np.sqrt(np.log(2.0))
BOLSGS = np.float64(1.3806503e-16)
SPL = np.float64(2.99792458e10)
AVOGADRO = np.float64(6.02214076e23)
C2 = np.float64(1.438777)
REF_TEMPERATURE = np.float64(296.0)
DOPPLER_CONST = np.sqrt(2.0 * AVOGADRO * BOLSGS * np.log(2.0)) / SPL

DELTA_WV = np.float64(10.0)
H_FINE = np.float64(1.0 / 2048.0)  # 0.00048828125, same output resolution as MarfaSimple.


@dataclass
class ParsedArgs:
    parfile: str
    maxlines: int
    molecule_int_code: int
    wv_min: float
    wv_max: float
    cutoff: float
    temperature: float
    p_self: float
    p_foreign: float
    tips_ref: float
    tips_target: float
    outfile: str


def parse_command_line(argv: list[str]) -> ParsedArgs:
    if len(argv) != 13:
        print(
            "Usage: marfa_python <parfile> <maxlines> <molec_id> <wv_min> <wv_max> "
            "<cutoff> <temp> <p_self> <p_foreign> <TIPS_ref> <TIPS_target> <outfile>"
        )
        raise SystemExit(1)
    return ParsedArgs(
        parfile=argv[1],
        maxlines=int(argv[2]),
        molecule_int_code=int(argv[3]),
        wv_min=float(argv[4]),
        wv_max=float(argv[5]),
        cutoff=float(argv[6]),
        temperature=float(argv[7]),
        p_self=float(argv[8]),
        p_foreign=float(argv[9]),
        tips_ref=float(argv[10]),
        tips_target=float(argv[11]),
        outfile=argv[12],
    )


def convert_iso(iso_char: str) -> int:
    if iso_char.isdigit():
        return 10 if iso_char == "0" else int(iso_char)
    if iso_char == "A":
        return 11
    return 1


def read_hitran_file(
    filename: str,
    max_lines: int,
    start_wv: float,
    end_wv: float,
    cutoff: float,
) -> tuple[np.ndarray, ...]:
    ext_min = start_wv - cutoff
    ext_max = end_wv + cutoff

    wv = []
    intensity_ref = []
    gamma_foreign = []
    gamma_self = []
    lower_state = []
    temp_coeff = []
    delta_foreign = []
    mol_iso = []

    with open(filename, "r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            if len(raw) < 67:
                continue
            try:
                line_wv = float(raw[3:15])
            except ValueError:
                continue
            if line_wv < ext_min or line_wv > ext_max:
                continue
            try:
                mo = int(raw[0:2])
                iso = raw[2:3]
                wv.append(line_wv)
                intensity_ref.append(float(raw[15:25]))
                gamma_foreign.append(float(raw[35:40]))
                gamma_self.append(float(raw[40:45]))
                lower_state.append(float(raw[45:55]))
                temp_coeff.append(float(raw[55:59]))
                delta_foreign.append(float(raw[59:67]))
                mol_iso.append(mo * 100 + convert_iso(iso))
            except ValueError:
                continue
            if len(wv) >= max_lines:
                break

    if not wv:
        empty_f = np.empty(0, dtype=np.float64)
        empty_i = np.empty(0, dtype=np.int32)
        return empty_f, empty_f, empty_f, empty_f, empty_f, empty_f, empty_f, empty_i

    order = np.argsort(np.asarray(wv, dtype=np.float64))
    return (
        np.asarray(wv, dtype=np.float64)[order],
        np.asarray(intensity_ref, dtype=np.float64)[order],
        np.asarray(gamma_foreign, dtype=np.float64)[order],
        np.asarray(gamma_self, dtype=np.float64)[order],
        np.asarray(lower_state, dtype=np.float64)[order],
        np.asarray(temp_coeff, dtype=np.float64)[order],
        np.asarray(delta_foreign, dtype=np.float64)[order],
        np.asarray(mol_iso, dtype=np.int32)[order],
    )


@njit(cache=True)
def get_molar_mass(mol_num: int) -> float:
    if mol_num == 1:
        return 18.015
    if mol_num == 2:
        return 44.01
    if mol_num == 3:
        return 47.998
    if mol_num == 4:
        return 44.013
    if mol_num == 5:
        return 28.01
    if mol_num == 6:
        return 16.04
    if mol_num == 7:
        return 31.999
    if mol_num == 8:
        return 30.01
    if mol_num == 9:
        return 64.066
    return 44.01


@njit(cache=True)
def lorentz(x: float, lor_hwhm: float, line_intensity: float) -> float:
    return lor_hwhm / (PI * (x * x + lor_hwhm * lor_hwhm)) * line_intensity


@njit(cache=True)
def doppler(x: float, dop_hwhm: float, line_intensity: float) -> float:
    return SQLN2 / (np.sqrt(PI) * dop_hwhm) * np.exp(-(x / dop_hwhm) ** 2 * np.log(2.0)) * line_intensity


@njit(cache=True)
def voigt_asymptotic1(x: float, lor_hwhm: float, vx: float, line_intensity: float) -> float:
    return lorentz(x, lor_hwhm, line_intensity) * (1.0 + 1.5 / (vx * vx))


@njit(cache=True)
def voigt_asymptotic2(x: float, lor_hwhm: float, vx: float, dop_hwhm: float, line_intensity: float) -> float:
    u = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0], dtype=np.float64)
    w = np.array([-0.688, 0.2667, 0.6338, 0.4405, 0.2529, 0.1601, 0.1131, 0.0853, 0.068], dtype=np.float64)
    idx = int(vx / 0.5 - 1.00001)
    if idx < 0:
        idx = 0
    if idx > 7:
        idx = 7
    f = 2.0 * (w[idx] * (u[idx + 1] - vx) + w[idx + 1] * (vx - u[idx]))
    return doppler(x, dop_hwhm, line_intensity) + (lor_hwhm / (PI * x * x) * (1.0 + f)) * line_intensity


@njit(cache=True)
def voigt_profile(x: float, dop_hwhm: float, lor_hwhm: float, line_intensity: float) -> float:
    vx = abs(SQLN2 * x / dop_hwhm)
    vy = SQLN2 * lor_hwhm / dop_hwhm
    if vx >= 15.0:
        return lorentz(x, lor_hwhm, line_intensity)

    vx2 = vx * vx

    if vx + vy >= 15.0:
        y1 = vy
        y2 = y1 * y1
        a1 = (0.2820948 + 0.5641896 * y2) * y1
        b1 = 0.5641896 * y1
        a2 = 0.25 + y2 + y2 * y2
        b2 = y2 + y2 - 1.0
        return (a1 + b1 * vx2) / (a2 + b2 * vx2 + vx2 * vx2) / np.sqrt(PI) / (dop_hwhm / SQLN2) * line_intensity

    if vx + vy >= 5.5:
        y = vy
        y2 = y * y
        a3 = y * (((0.56419 * y2 + 3.10304) * y2 + 4.65456) * y2 + 1.05786)
        b3 = y * ((1.69257 * y2 + 0.56419) * y2 + 2.962)
        c3 = y * (1.69257 * y2 - 2.53885)
        d3 = y * 0.56419
        a4 = (((y2 + 6.0) * y2 + 10.5) * y2 + 4.5) * y2 + 0.5625
        b4 = ((4.0 * y2 + 6.0) * y2 + 9.0) * y2 - 4.5
        c4 = 10.5 + 6.0 * (y2 - 1.0) * y2
        d4 = 4.0 * y2 - 6.0
        return (
            (((d3 * vx2 + c3) * vx2 + b3) * vx2 + a3)
            / ((((vx2 + d4) * vx2 + c4) * vx2 + b4) * vx2 + a4)
            / np.sqrt(PI)
            / (dop_hwhm / SQLN2)
            * line_intensity
        )

    if vx <= 1.0 or vy >= 0.02:
        y = vy
        a5 = ((((((((0.564224 * y + 7.55895) * y + 49.5213) * y + 204.510) * y + 581.746) * y + 1174.8) * y + 1678.33) * y + 1629.76) * y + 973.778) * y + 272.102
        b5 = ((((((2.25689 * y + 22.6778) * y + 100.705) * y + 247.198) * y + 336.364) * y + 220.843) * y - 2.34403) * y - 60.5644
        c5 = ((((3.38534 * y + 22.6798) * y + 52.8454) * y + 42.5683) * y + 18.546) * y + 4.58029
        d5 = ((2.25689 * y + 7.56186) * y + 1.66203) * y - 0.128922
        e5 = 0.000971457 + 0.564224 * y
        a6 = (((((((((y + 13.3988) * y + 88.2674) * y + 369.199) * y + 1074.41) * y + 2256.98) * y + 3447.63) * y + 3764.97) * y + 2802.87) * y + 1280.83) * y + 272.102
        b6 = (((((((5.0 * y + 53.5952) * y + 266.299) * y + 793.427) * y + 1549.68) * y + 2037.31) * y + 1758.34) * y + 902.306) * y + 211.678
        c6 = (((((10.0 * y + 80.3928) * y + 269.292) * y + 479.258) * y + 497.302) * y + 308.186) * y + 78.866
        d6 = (((10.0 * y + 53.5952) * y + 92.7586) * y + 55.0293) * y + 22.0353
        e6 = (5.0 * y + 13.3988) * y + 1.49645
        return (
            ((((e5 * vx2 + d5) * vx2 + c5) * vx2 + b5) * vx2 + a5)
            / (((((vx2 + e6) * vx2 + d6) * vx2 + c6) * vx2 + b6) * vx2 + a6)
            / np.sqrt(PI)
            / (dop_hwhm / SQLN2)
            * line_intensity
        )

    if vx > 5.0:
        return voigt_asymptotic1(x, lor_hwhm, vx, line_intensity)
    if vx > np.sqrt(1.4):
        return voigt_asymptotic2(x, lor_hwhm, vx, dop_hwhm, line_intensity)
    return doppler(x, dop_hwhm, line_intensity)


@njit(cache=True)
def line_parameters(
    wv: np.ndarray,
    intensity_ref: np.ndarray,
    gamma_foreign: np.ndarray,
    gamma_self: np.ndarray,
    lower_state: np.ndarray,
    temp_coeff: np.ndarray,
    delta_foreign: np.ndarray,
    mol_iso: np.ndarray,
    temperature: float,
    p_self: float,
    p_foreign: float,
    tips_ref: float,
    tips_target: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = wv.shape[0]
    shifted = np.empty(n, dtype=np.float64)
    dop = np.empty(n, dtype=np.float64)
    lor = np.empty(n, dtype=np.float64)
    intens = np.empty(n, dtype=np.float64)

    tips_factor = tips_ref / tips_target
    pressure = p_self + p_foreign
    for i in range(n):
        shifted[i] = wv[i] + delta_foreign[i] * pressure
        mol = mol_iso[i] // 100
        mass = get_molar_mass(mol)
        dop[i] = DOPPLER_CONST * shifted[i] * np.sqrt(temperature / mass)
        lor[i] = ((REF_TEMPERATURE / temperature) ** temp_coeff[i]) * (gamma_foreign[i] * p_foreign + gamma_self[i] * p_self)
        boltz = np.exp(-C2 * lower_state[i] / temperature) / np.exp(-C2 * lower_state[i] / REF_TEMPERATURE)
        stim = (1.0 - np.exp(-C2 * wv[i] / temperature)) / (1.0 - np.exp(-C2 * wv[i] / REF_TEMPERATURE))
        intens[i] = intensity_ref[i] * tips_factor * boltz * stim
    return shifted, dop, lor, intens


@njit(cache=True)
def _clamp_index(lo: int, hi: int, size: int) -> tuple[int, int]:
    if lo < 0:
        lo = 0
    if hi >= size:
        hi = size - 1
    return lo, hi


@njit(cache=True)
def _accumulate_line_on_grid(
    acc: np.ndarray,
    step: float,
    start_wv: float,
    center: float,
    left: float,
    right: float,
    dop_hwhm: float,
    lor_hwhm: float,
    line_intensity: float,
) -> None:
    if right <= left:
        return
    lo = int(np.ceil((center - right - start_wv) / step))
    hi = int(np.floor((center + right - start_wv) / step))
    lo, hi = _clamp_index(lo, hi, acc.shape[0])
    for i in range(lo, hi + 1):
        nu = start_wv + step * i
        x = nu - center
        ax = abs(x)
        if ax <= left or ax > right:
            continue
        acc[i] += voigt_profile(x, dop_hwhm, lor_hwhm, line_intensity)


@njit(cache=True)
def _accumulate_line_on_fine(
    acc: np.ndarray,
    step: float,
    start_wv: float,
    center: float,
    radius: float,
    dop_hwhm: float,
    lor_hwhm: float,
    line_intensity: float,
) -> None:
    if radius <= 0.0:
        return
    lo = int(np.ceil((center - radius - start_wv) / step))
    hi = int(np.floor((center + radius - start_wv) / step))
    lo, hi = _clamp_index(lo, hi, acc.shape[0])
    for i in range(lo, hi + 1):
        nu = start_wv + step * i
        x = nu - center
        if abs(x) <= radius:
            acc[i] += voigt_profile(x, dop_hwhm, lor_hwhm, line_intensity)


@njit(cache=True)
def _cascade_add_factor2(coarse: np.ndarray, fine: np.ndarray) -> None:
    n_coarse = coarse.shape[0]
    n_fine = fine.shape[0]
    for j in range(n_coarse):
        i = 2 * j
        if i < n_fine:
            fine[i] += coarse[j]
        i_mid = i + 1
        if i_mid < n_fine:
            if j + 1 < n_coarse:
                fine[i_mid] += 0.5 * (coarse[j] + coarse[j + 1])
            else:
                fine[i_mid] += coarse[j]


@njit(cache=True)
def marfa_3grid_numba(
    shifted: np.ndarray,
    dop: np.ndarray,
    lor: np.ndarray,
    intens: np.ndarray,
    start_wv: float,
    end_wv: float,
    cutoff: float,
    h_fine: float,
) -> tuple[np.ndarray, np.ndarray]:
    n0 = int(np.floor((end_wv - start_wv) / h_fine + 1e-12)) + 1
    h1 = h_fine * 2.0
    h2 = h_fine * 4.0
    n1 = int(np.floor((end_wv - start_wv) / h1 + 1e-12)) + 1
    n2 = int(np.floor((end_wv - start_wv) / h2 + 1e-12)) + 1

    rk0 = np.zeros(n0, dtype=np.float64)  # fine
    rk1 = np.zeros(n1, dtype=np.float64)  # medium
    rk2 = np.zeros(n2, dtype=np.float64)  # coarse

    # In MarfaSimple, lines inside a 10 cm-1 subinterval are spread over that full subinterval.
    # So cutoff=0 still implies a non-zero profile support for "center" lines.
    wing_limit = cutoff
    if wing_limit < DELTA_WV:
        wing_limit = DELTA_WV

    # Three-grid partition radii.
    r0 = 2.0
    r1 = 8.0
    if r1 > wing_limit:
        r1 = wing_limit
    if r0 > r1:
        r0 = r1

    for i in range(shifted.shape[0]):
        line_intensity = intens[i]
        if line_intensity <= 0.0:
            continue
        center = shifted[i]
        if center < start_wv - wing_limit or center > end_wv + wing_limit:
            continue

        _accumulate_line_on_fine(rk0, h_fine, start_wv, center, min(r0, wing_limit), dop[i], lor[i], line_intensity)
        _accumulate_line_on_grid(rk1, h1, start_wv, center, min(r0, wing_limit), min(r1, wing_limit), dop[i], lor[i], line_intensity)
        _accumulate_line_on_grid(rk2, h2, start_wv, center, min(r1, wing_limit), wing_limit, dop[i], lor[i], line_intensity)

    # Cascade interpolation from coarse -> medium -> fine.
    _cascade_add_factor2(rk2, rk1)
    _cascade_add_factor2(rk1, rk0)

    wavenumbers = np.empty(n0, dtype=np.float64)
    for i in range(n0):
        wavenumbers[i] = start_wv + h_fine * i
    return wavenumbers, rk0


def run(args: ParsedArgs) -> None:
    total_cpu_start = time.process_time()
    total_wall_start = time.perf_counter()

    section_cpu_start = time.process_time()
    section_wall_start = time.perf_counter()
    (
        wv,
        intensity_ref,
        gamma_foreign,
        gamma_self,
        lower_state,
        temp_coeff,
        delta_foreign,
        mol_iso,
    ) = read_hitran_file(args.parfile, args.maxlines, args.wv_min, args.wv_max, args.cutoff)
    section_cpu_end = time.process_time()
    section_wall_end = time.perf_counter()
    print("================================")
    print("=== TIMING: readHITRANFile ===")
    print(f"CPU time:  {section_cpu_end - section_cpu_start:.6f} seconds")
    print(f"Wall clock time:  {section_wall_end - section_wall_start:.6f} seconds")
    print(f"Lines read:  {wv.shape[0]}")
    print("================================")

    section_cpu_start = time.process_time()
    section_wall_start = time.perf_counter()
    shifted, dop, lor, intens = line_parameters(
        wv,
        intensity_ref,
        gamma_foreign,
        gamma_self,
        lower_state,
        temp_coeff,
        delta_foreign,
        mol_iso,
        args.temperature,
        args.p_self,
        args.p_foreign,
        args.tips_ref,
        args.tips_target,
    )
    section_cpu_end = time.process_time()
    section_wall_end = time.perf_counter()
    print("================================")
    print("=== TIMING: lineParameters ===")
    print(f"CPU time:  {section_cpu_end - section_cpu_start:.6f} seconds")
    print(f"Wall clock time:  {section_wall_end - section_wall_start:.6f} seconds")
    print("================================")

    section_cpu_start = time.process_time()
    section_wall_start = time.perf_counter()
    nu, xsc = marfa_3grid_numba(shifted, dop, lor, intens, args.wv_min, args.wv_max, args.cutoff, H_FINE)
    section_cpu_end = time.process_time()
    section_wall_end = time.perf_counter()
    print("================================")
    print("=== TIMING: marfa_3grid_numba ===")
    print(f"CPU time:  {section_cpu_end - section_cpu_start:.6f} seconds")
    print(f"Wall clock time:  {section_wall_end - section_wall_start:.6f} seconds")
    print("================================")

    section_cpu_start = time.process_time()
    section_wall_start = time.perf_counter()
    with open(args.outfile, "w", encoding="utf-8") as handle:
        handle.write("# Wavenumber [cm-1]    Absorption cross-section [cm2/molec]\n")
        for n, s in zip(nu, xsc):
            handle.write(f"{n:12.6f} {s:15.6E}\n")
    section_cpu_end = time.process_time()
    section_wall_end = time.perf_counter()
    print("================================")
    print("=== TIMING: writeOutput ===")
    print(f"CPU time:  {section_cpu_end - section_cpu_start:.6f} seconds")
    print(f"Wall clock time:  {section_wall_end - section_wall_start:.6f} seconds")
    print("================================")

    total_cpu_end = time.process_time()
    total_wall_end = time.perf_counter()
    print("================================")
    print("=== TIMING: MarfaPython total ===")
    print(f"CPU time:  {total_cpu_end - total_cpu_start:.6f} seconds")
    print(f"Wall clock time:  {total_wall_end - total_wall_start:.6f} seconds")
    print("================================")
    print("NOTE: first run includes Numba JIT compilation overhead.")


if __name__ == "__main__":
    run(parse_command_line(sys.argv))
