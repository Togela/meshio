"""
I/O for the LS-DYNA .k file format.
"""
import datetime

import numpy as np

from ..__about__ import __version__
from .._exceptions import WriteError
from .._files import open_file
from .._helpers import register_format
from .._mesh import CellBlock, Mesh

def read(filename):
    with open_file(filename, "r") as f:
        mesh = read_buffer(f)
    return mesh

def read_buffer(f):
    # Initialize the optional data fields
    points = []
    cells = []
    cell_ids = []
    field_data = {}
    cell_data = {}
    point_data = {}
    point_ids = None

    line = f.readline()
    while True:
        if not line:  # EOF
            break

        # Comments
        if line.startswith("$"):
            line = f.readline()
            continue

        # NODE
        if line.startswith("*NODE"):
            points, point_ids, line = _read_nodes(f)
            continue
        elif line.startswith == "*ELEMENT_BEAM":
            # TODO continue
            # cell_type, cells_data, ids, sets, line = _read_element_beam(f)
            cells.append(CellBlock(cell_type, cells_data))
            cell_ids.append(ids)
        elif line.startswith == "*ELEMENT_SHELL":
            # TODO continue
            cell_type, cells_data, ids, sets, line = _read_element_beam(f)
        elif line.startswith == "*ELEMENT_SOLID":
            # TODO continue
            

def _read_nodes(f):
    points = []
    point_ids = {}
    counter = 0
    while True:
        line = f.readline()
        if not line or line.startswith("*"):
            break

        if "," in line: # , structure
            line = line.strip().split(",")
        else: # standard structure
            line = [line[:8],line[8:24], line[24:40], line[40:56] ]
        point_id, coords = line[0], line[1:]
        point_ids[int(point_id)] = counter
        points.append([float(x) for x in coords])
        counter += 1

    return np.array(points, dtype=float), point_ids, line

def _read_element_beam(f):
    beam = []
    beam_ids = {}
    counter = 0

    while True:
        line = f.readline()
        if not line or line.startswith("*"):
            break

        if "," in line: # , structure
            line = line.strip().split(",")
        else: # standard structure
            # TODO: two node & three node structure
            line = [line[:8],line[8:24], line[24:40], line[40:56] ]
        point_id, coords = line[0], line[1:]
        point_ids[int(point_id)] = counter
        points.append([float(x) for x in coords])
        counter += 1

    return np.array(points, dtype=float), point_ids, line
