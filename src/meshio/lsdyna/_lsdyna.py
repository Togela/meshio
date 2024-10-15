"""
I/O for the LS-DYNA .k file format.
"""
import datetime

import numpy as np

from ..__about__ import __version__
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._helpers import register_format
from .._mesh import CellBlock, Mesh
from .._common import warn

lineskip_kw_element_beam = {'*ELEMENT_BEAM': 1, 
                            '*ELEMENT_BEAM_ELBOW': 2, 
                            '*ELEMENT_BEAM_OFFSET': 2, 
                            '*ELEMENT_BEAM_ORIENTATION_OFFSET': 3,
                            '*ELEMENT_BEAM_ORIENTATION': 2,
                            '*ELEMENT_BEAM_PID': 2,
                            '*ELEMENT_BEAM_PID_OFFSET': 3,
                            '*ELEMENT_BEAM_PID_ORIENTATION': 3,
                            '*ELEMENT_BEAM_SCALAR': 2,
                            '*ELEMENT_BEAM_SCALAR_PID': 3,
                            '*ELEMENT_BEAM_SCALAR_OFFSET': 3,
                            '*ELEMENT_BEAM_SCALAR_ORIENTATION': 3,
                            '*ELEMENT_BEAM_SECTION': 2,
                            '*ELEMENT_BEAM_SECTION_SCALAR': 3,
                            '*ELEMENT_BEAM_SECTION_PID': 3,
                            '*ELEMENT_BEAM_SECTION_OFFSET': 3,
                            '*ELEMENT_BEAM_SECTION_ORIENTATION': 3,
                            '*ELEMENT_BEAM_THICKNESS': 2,
                            '*ELEMENT_BEAM_THICKNESS_SCALAR': 3,
                            '*ELEMENT_BEAM_THICKNESS_PID': 3,
                            '*ELEMENT_BEAM_THICKNESS_OFFSET': 3,
                            '*ELEMENT_BEAM_THICKNESS_ORIENTATION': 3,
                            '*ELEMENT_BEAM_WARPAGE': 2
                            }

lineskip_kw_element_shell = {'*ELEMENT_SHELL': 1,
                             '*ELEMENT_SHELL_SHL4_TO_SHL8': 1,
                             '*ELEMENT_SHELL_BETA': 3,
                             '*ELEMENT_SHELL_BETA_OFFSET': 4,
                             '*ELEMENT_SHELL_BEXT_PATCH': 8,
                             '*ELEMENT_SHELL_COMPOSITE': 2,
                             '*ELEMENT_SHELL_COMPOSITE_LONG': 2,
                             '*ELEMENT_SHELL_BETA_COMPOSITE': 4,
                             '*ELEMENT_SHELL_BETA_COMPOSITE_LONG': 4,
                             '*ELEMENT_SHELL_DOF': 2,
                             '*ELEMENT_SHELL_MCID': 3,
                             '*ELEMENT_SHELL_MCID_OFFSET': 4,
                             '*ELEMENT_SHELL_OFFSET': 2,
                             '*ELEMENT_SHELL_OFFSET_COMPOSITE': 3,
                             '*ELEMENT_SHELL_OFFSET_COMPOSITE_LONG': 3,
                             '*ELEMENT_SHELL_BETA_OFFSET_COMPOSITE': 3,
                             '*ELEMENT_SHELL_BETA_OFFSET_COMPOSITE_LONG': 3,
                             '*ELEMENT_SHELL_THICKNESS': 3,
                             '*ELEMENT_SHELL_THICKNESS_BETA': 3,
                             '*ELEMENT_SHELL_THICKNESS_BETA_OFFSET': 4,
                             '*ELEMENT_SHELL_THICKNESS_MCID': 3,
                             '*ELEMENT_SHELL_THICKNESS_MCID_OFFSET': 3,
                             '*ELEMENT_SHELL_THICKNESS_OFFSET': 4
                             }

lineskip_kw_element_solid = {'*ELEMENT_SOLID': 1,
                             '*ELEMENT_SOLID_DOF': 2,
                             '*ELEMENT_SOLID_PERI': 1,
                             '*ELEMENT_SOLID_NURBS_PATCH': 3,
                             '*ELEMENT_SOLID_TET4TOTET10': 1,
                             '*ELEMENT_SOLID_TET4TOTET10_ORTHO': 3,
                             '*ELEMENT_SOLID_TET4TOTET10_ORTHO_DOF': 4,
                             '*ELEMENT_SOLID_TET4TOTET10_FORMAT2': 2,
                             '*ELEMENT_SOLID_TET4TOTET10_DOF_FORMAT2': 3,
                             '*ELEMENT_SOLID_ORTHO': 3,
                             '*ELEMENT_SOLID_ORTHO_DOF': 4,
                             '*ELEMENT_SOLID (ten nodes format)': 2,
                             '*ELEMENT_SOLID_DOF (ten nodes format)': 3,
                             '*ELEMENT_SOLID_ORTHO (ten nodes format)': 4,
                             '*ELEMENT_SOLID_ORTHO_DOF (ten nodes format)': 5,
                             '*ELEMENT_SOLID (ten nodes format)': 2
                             # '*ELEMENT_SOLID_H8TOH20': 1 are skipped
                             }

lsdyna_to_meshio_type = {
    '*NODE': "vertex",
    '*ELEMENT_BEAM': "line",
    '*ELEMENT_SHELL': ("triangle", "quad"),
    '*ELEMENT_SOLID': ("tetra", "hexahedron")
    }

lsdyna_to_meshio_nnids = {
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "triangle6": 6,
    "quad8": 8,
    "tetra": 4,
    "hexahedron": 8
    }

meshio_to_lsdyna_type = {}
for k, v in lsdyna_to_meshio_type.items():
    if type(v) is str:
        meshio_to_lsdyna_type[v] = k
    else:
        for _v in v:
            meshio_to_lsdyna_type[_v] = k


def read(filename):
    with open_file(filename, "r") as f:
        mesh = read_buffer(f)
    return mesh

def read_buffer(f):
    # Initialize the optional data fields
    points = []
    cells = []
    cell_ids = []
    cell_pids = []
    cell_blocks = []
    field_data = {}
    cell_data = {}
    point_data = {}
    point_ids = []
    point_sets = []
    cell_sets = []
    
    line = f.readline()
    while True:
        if not line:  # EOF
            break

        # Comments
        if line.startswith("$"):
            line = f.readline()
            continue

        if line.startswith("*NODE"):
            # Node
            _points, _point_ids, line = _read_nodes(f)
            points.extend(_points)
            point_ids.extend(_point_ids)
            continue
        elif line.startswith("*ELEMENT"):
            # Beam element
            # TODO continue               
            kw = line.rstrip()         
            _cell_block, _cell_ids, _element_ids, _element_pids, line = _read_element(f, kw=kw)
            for k, v in _cell_block.items():
                cell_blocks.append([k, np.array(v)])
            # cell_ids.append(_cell_ids)
            cell_ids.append(_element_ids)
            cell_pids.append(_element_pids)
        elif line.startswith("*END"):
            break
        else:
            line = f.readline()
            continue

    # Convert to natural point ordering
    # https://stackoverflow.com/questions/16992713/translate-every-element-in-numpy-array-according-to-key
    point_ids_ext_to_int = {v:k for k,v in enumerate(point_ids)}
    conversion_point_ids = np.vectorize(point_ids_ext_to_int.__getitem__)
    for block in cell_blocks:
        block_natural = conversion_point_ids(block[1])
        cells.append(CellBlock(block[0], block_natural))
        
    point_data['lsdyna:nid'] = point_ids
    cell_data['lsdyna:eid'] = cell_ids
    cell_data['lsdyna:pid'] = cell_pids
    mesh =  Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
        point_sets=point_sets,
        cell_sets=cell_sets,
    )
    return mesh

def _read_nodes(f):
    points = []
    point_ids = []
    counter = 0
    while True:
        line = f.readline()
        
        if line.startswith("$"): # skip comments
            continue
        
        if not line or line.startswith("*"): 
            break

        # read node data
        if "," in line: # , structure
            data = line.strip().split(",")
        elif len(line) <= 80: # l10 structure
            data = [line[:8],line[8:28],line[28:48],line[48:68]]
        elif len(line) <= 120: # long structure
            data = [line[:20], line[20:40], line[40:60], line[60:80]]
        else:
            raise ReadError('Unknown line lenght')
        point_id, coords = data[0], data[1:4]
        point_ids.append(int(point_id))
        points.append([float(x) for x in coords])
        counter += 1

    return np.array(points, dtype=float), point_ids, line

def _read_element(f:str, kw:str):#
    cell_block = {}
    cell_ids = {}
    
    # element_nids_type1 = []
    # element_nids_type2 = []
    # element_pids_type1 = []
    # element_pids_type2 = []
    # element_types = []
    element_ids = []
    element_pids = []
    counter = 0

    if kw.startswith('*ELEMENT_BEAM'):
        lineskip = lineskip_kw_element_beam.get(kw, 1)
        nnids = [3, 2] # [with 3rd beam node, without] 
        cell_type = ["line", "line"]
    elif kw.startswith('*ELEMENT_SHELL'):
        lineskip = lineskip_kw_element_shell.get(kw, 1)
        nnids = [4, 3] # [quad, triangle]
        cell_type = ["quad", "triangle"]
    if kw.startswith('*ELEMENT_SOLID'):
        lineskip = lineskip_kw_element_solid.get(kw, 1)
        nnids = [8, 4] # [hexahedron, tetra]
        cell_type = ["hexahedron", "tetra"]
        
    iskip = 1
    while True:
        line = f.readline()
        if not line or line.startswith("*"):
            break
        
        if line.startswith("$"): # skip comments
            continue

        if iskip != 1: # skip lines for certain keywords
            iskip += 1
            if iskip == lineskip:
                iskip = 1
            continue

        line = line.rstrip()
        if "," in line: # , fixed format        
            # only two node are taken from two node & three node structures
            data = line.strip().split(",")             # elid, pid, nid1, nid2 
        elif len(line) <= 80: # l10 structure
            try: 
                data = [line[(idata*8):((idata+1)*8)] for idata in range(nnids[0]+2)] # elid, pid, nid1, nid2, nid3 ...
            except IndexError:
                data = [line[(idata*8):((idata+1)*8)] for idata in range(nnids[1]+2)] # elid, pid, nid1, nid2, nid3 ...                
        elif len(line) <= 200: # long structure
            try: 
                data = [line[(idata*20):((idata+1)*20)] for idata in range(nnids[0]+2)] # elid, pid, nid1, nid2, nid3 ...
            except IndexError:
                data = [line[(idata*20):((idata+1)*20)] for idata in range(nnids[1]+2)] # elid, pid, nid1, nid2, nid3 ...
        else:
            raise ReadError('Unknown line lenght')

        data = [d for d in data if d != ""]
        element_id, part_id = int(data[0]), int(data[1])
        element_point_ids = [int(d) for d in data[2:]]
        if len(element_point_ids) == nnids[1]:
            element_type = cell_type[1]  
        elif all(np.equal(element_point_ids[-(nnids[0]-nnids[1]):],0)): # ls dyna might fill the empty cells with 0s
            element_type = cell_type[1] 
            element_point_ids = element_point_ids[:nnids[1]]
        else:
            element_type = cell_type[0]
                   
        if element_type in cell_block:
            cell_block[element_type].append(element_point_ids)
            cell_ids[element_type].append(counter)
        else:
            cell_block[element_type] = [element_point_ids]
            cell_ids[element_type] = [counter]
        element_ids.append(element_id)
        element_pids.append(part_id)

        counter += 1

    return cell_block, cell_ids, element_ids, element_pids, line

def floattochars(number, nchars, exp_digits):
    if isinstance(number, str):
        return number[:16]
    string_fixed = str(number)[:nchars]
    string_scientific = np.format_float_scientific(number, nchars-5-exp_digits, exp_digits=exp_digits)
    error_fixed = abs(float(string_fixed) - number)
    error_scientific = abs(float(string_scientific) - number)

    if error_fixed < error_scientific:
        return string_fixed
    else:
        return string_scientific

def write(filename, mesh, format="fixed"):
    
    if format == "fixed":
        blocksize = 8
        nchars = 16
        elementfmt = {"vertex": "{:>" + str(blocksize) + "}"+"{:>20}"*3+"\n",
                      "line": ("{:>"+str(blocksize)+"}")*4+"\n",
                      "triangle":("{:>"+str(blocksize)+"}")*5+"\n",
                      "quad":("{:>"+str(blocksize)+"}")*6+"\n",
                      "tetra":("{:>"+str(blocksize)+"}")*6+"\n",
                      "hexahedron":("{:>"+str(blocksize)+"}")*10+"\n",
                      }
        setfmt = ["{:>"+str(blocksize)+"}",
                  ("{:>"+str(blocksize)+"}")*10+"\n"]
    elif format == "fixedlong":
        blocksize = 20
        nchars = 16
        elementfmt = {"vertex": "{:>" + str(blocksize) + "}"+"{:>20}"*3+"\n",
                      "line": ("{:>"+str(blocksize)+"}")*4+"\n",
                      "triangle":("{:>"+str(blocksize)+"}")*5+"\n",
                      "quad":("{:>"+str(blocksize)+"}")*6+"\n",
                      "tetra":("{:>"+str(blocksize)+"}")*6+"\n",
                      "hexahedron":("{:>"+str(blocksize)+"}")*10+"\n",
                      }
        setfmt = ["{:>"+str(blocksize)+"}",
                  ("{:>"+str(blocksize)+"}")*10+"\n"]
    elif format == "free":
        nchars = 12
        elementfmt = {"vertex": ",".join(["{:}"]*4)+"\n",
                      "line": ",".join(["{:}"]*4)+"\n",
                      "triangle": ",".join(["{:}"]*5)+"\n",
                      "quad": ",".join(["{:}"]*6)+"\n",
                      "tetra": ",".join(["{:}"]*6)+"\n",
                      "hexahedron": ",".join(["{:}"]*10)+"\n",
                      }
        setfmt = ["{:}",
                  ",".join(["{:}"]*10)+"\n"]
    else:
        raise RuntimeError(f'unknown "{format}" format')

    def floatfmt(number): return floattochars(number=number, nchars=nchars, exp_digits=2) 

    with open_file(filename, "w") as f:
        f.write(f"$ LS-DYNA keyword file written by meshio v{__version__}\n")
        f.write("*KEYWORD\n")
        
        # Write nodes
        point_refs = mesh.point_data.get("lsdyna:nid", None)
        if len(mesh.points) > 0:
            f.write('*NODE\n')

            if mesh.points.shape[1] == 2:
                warn("2D mesh is converted to 3D (z-value is 0)")
                points = np.pad(mesh.points, [[0,0],[0,1]], mode='constant', constant_values=0)
            else: 
                points = mesh.points
            for point_id, x in enumerate(points, 1):
                fx = [floatfmt(k) for k in x]
                pref = str(point_refs[point_id]) if point_refs is not None else point_id
                string = elementfmt['vertex'].format(pref, *fx)
                f.write(string)

        # Write elements
        element_eids = mesh.cell_data.get("lsdyna:eid", None)  
        element_pids = mesh.cell_data.get("lsdyna:pid", np.ones(len(mesh.cells), dtype=int))
        eid = 1
        for ict, cell_block in enumerate(mesh.cells):
            cell_type = cell_block.type
            cells = cell_block.data
            assert cell_type in meshio_to_lsdyna_type, f"Cell type {cell_type} currently not supported"    
            keyword_name = meshio_to_lsdyna_type[cell_type]
            elfmt = elementfmt[cell_type]
                
            cell_pid = str(element_pids[ict])
            f.write(keyword_name+"\n")
            for row in cells:
                if element_eids is None:
                    cell_elid = str(eid)
                    eid += 1
                else:                
                    cell_elid = str(element_eids[ict])
                nids_strs = (str(nid + 1) for nid in row)
                f.write(elfmt.format(cell_elid, cell_pid, *nids_strs))
    
        for k, v in mesh.point_sets.items():
            nds = [str(i + 1) for i in v]
            f.write("*SET_NODE_LIST\n")
            f.write(setfmt[0].format(k))
            while len(nds > 10):
                _nds = nds.pop(10)
                f.write(setfmt[1].format(*_nds))
            _nds = np.pad(_nds, (0, 10-len(_nds)), mode='constant', constant_values="")
            f.write(setfmt[1].format(*_nds))
    
        for ic in range(len(mesh.cells)):
            cell_type = meshio_to_lsdyna_type[mesh.cells[ic].type]
            if cell_type == '*ELEMENT_BEAM':
                setstr = "*SET_BEAM\n" 
            elif cell_type == '*ELEMENT_SHELL':
                setstr = "*SET_SHELL\n" 
            elif cell_type == '*ELEMENT_SOLID':
                setstr = "*SET_SOLID\n" 
            else:
                continue
            for k, v in mesh.cell_sets.items():
                if len(v[ic]) > 0:
                    els = [str(i + 1) for i in v[ic]]
                    f.write(setstr)
                    f.write(setfmt[0].format(k))
                    while len(els > 10):
                        _els = els.pop(10)
                        f.write(setfmt[1].format(*_els))
                    _els = np.pad(_els, (0, 10-len(_els)), mode='constant', constant_values="")
                    f.write(setfmt[1].format(*_els))
                                
        f.write("*END\n")
        

if __name__ == "__main__":
    import meshio
    
    empty_mesh = meshio.Mesh(np.empty((0, 3)), [])

    line_mesh = meshio.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        [("line", [[0, 1], [0, 2], [0, 3], [1, 2], [2, 3]])],
)
    
    write('/tmp/test_mesh0_0', empty_mesh)
    write('/tmp/test_mesh0_1', empty_mesh)
