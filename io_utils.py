# io_utils.py
import re
from presests import TORICSURFACE_PRESETS


def parse_rays(input_str):
    """Parse a string of the form '(1, 0), (0, 1), (-1, -1)' into a list of tuples."""
    matches = re.findall(r"\(([^)]+)\)", input_str)
    rays = []
    for m in matches:
        parts = m.split(",")
        if len(parts) != 2:
            raise ValueError(f"Invalid ray: {m}")
        rays.append(tuple(float(p.strip()) for p in parts))
    return rays


def get_surface_rays(surface_key):
    """Return the raw ray string for a preset key, e.g. 'P2Sigma0'."""
    return TORICSURFACE_PRESETS[surface_key]