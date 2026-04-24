# io_utils.py
import re

def parse_rays(input_str):
    """
    Parse user input of the form:
    (1, 0), (0, 1), (-1, -1)
    """

    # extract all tuples
    matches = re.findall(r"\(([^)]+)\)", input_str)

    rays = []

    for m in matches:
        parts = m.split(",")

        if len(parts) != 2:
            raise ValueError(f"Invalid ray: {m}")

        ray = tuple(float(p.strip()) for p in parts)
        rays.append(ray)

    return rays
