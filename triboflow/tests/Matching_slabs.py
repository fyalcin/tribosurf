#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 12:56:41 2021

@author: wolloch
"""
from mpinterfaces.transformations import get_aligned_lattices
from pymatgen.core.surface import Slab
from triboflow.phys.interface_matcher import InterfaceMatcher

pt111_dict = {'@module': 'pymatgen.core.surface',
 '@class': 'Slab',
 'charge': None,
 'lattice': {'matrix': [[2.8120010342192265, 0.0, 1.7218540328778114e-16],
   [-1.4060005171096142, 2.435264331101964, 1.7218540328778114e-16],
   [0.0, 0.0, 41.32780614009411]],
  'a': 2.8120010342192265,
  'b': 2.8120010342192265,
  'c': 41.32780614009411,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 120.00000000000001,
  'volume': 283.01140376607475},
 'sites': [{'species': [{'element': 'Pt', 'occu': 1}],
   'abc': [0.3333333333333336, 0.6666666666666666, 0.36111111111111116],
   'xyz': [1.6305477402018132e-16, 1.623509554067976, 14.923929995033987],
   'label': 'Pt',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Pt', 'occu': 1}],
   'abc': [4.440892098500626e-16, 5.551115123125783e-16, 0.4722222222222222],
   'xyz': [4.682922440189902e-16, 1.3518432657188906e-15, 19.515908455044443],
   'label': 'Pt',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Pt', 'occu': 1}],
   'abc': [0.6666666666666667, 0.3333333333333338, 0.41666666666666663],
   'xyz': [1.4060005171096126, 0.8117547770339892, 17.219919225039213],
   'label': 'Pt',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Pt', 'occu': 1}],
   'abc': [0.3333333333333339, 0.6666666666666669, 0.5277777777777778],
   'xyz': [8.500606668368282e-16, 1.6235095540679765, 21.811897685049672],
   'label': 'Pt',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Pt', 'occu': 1}],
   'abc': [0.9999999999999999, 0.9999999999999997, 0.6388888888888888],
   'xyz': [1.4060005171096124, 2.4352643311019633, 26.403876145060124],
   'label': 'Pt',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Pt', 'occu': 1}],
   'abc': [0.6666666666666661, 0.3333333333333327, 0.5833333333333333],
   'xyz': [1.4060005171096122, 0.8117547770339866, 24.107886915054895],
   'label': 'Pt',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}}],
 'oriented_unit_cell': {'@module': 'pymatgen.core.structure',
  '@class': 'Structure',
  'charge': None,
  'lattice': {'matrix': [[1.9883849999999998, -1.988385, -0.0],
    [2.220446049250313e-16, 1.988385, -1.9883849999999998],
    [3.9767700000000006, 3.97677, 3.9767700000000006]],
   'a': 2.8120010342192265,
   'b': 2.8120010342192265,
   'c': 6.887967690015685,
   'alpha': 90.0,
   'beta': 90.0,
   'gamma': 120.00000000000001,
   'volume': 47.168567294345806},
  'sites': [{'species': [{'element': 'Pt', 'occu': 1}],
    'abc': [0.33333333333333326, 0.6666666666666667, 0.33333333333333326],
    'xyz': [1.9883849999999998, 1.988385, -1.4631037122588476e-16],
    'label': 'Pt',
    'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
   {'species': [{'element': 'Pt', 'occu': 1}],
    'abc': [0.9999999999999998, 0.0, 6.195649072435257e-17],
    'xyz': [1.9883849999999996, -1.9883849999999994, 2.463867136178836e-16],
    'label': 'Pt',
    'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
   {'species': [{'element': 'Pt', 'occu': 1}],
    'abc': [0.6666666666666663, 0.3333333333333336, 0.6666666666666666],
    'xyz': [3.9767699999999992, 1.9883850000000012, 1.9883849999999998],
    'label': 'Pt',
    'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}}]},
 'miller_index': (1, 1, 1),
 'shift': 0.16666666666666677,
 'scale_factor': [[-1, 1, 0], [-1, 0, 1], [-1, -1, -1]],
 'reconstruction': None,
 'energy': 0.0}
ag111_dict = {'@module': 'pymatgen.core.surface',
 '@class': 'Slab',
 'charge': None,
 'lattice': {'matrix': [[2.9419517042521277, 0.0, 1.8014258689292344e-16],
   [-1.470975852126065, 2.547804912589266, 1.8014258689292347e-16],
   [0.0, 0.0, 43.23768313997447]],
  'a': 2.9419517042521277,
  'b': 2.941951704252128,
  'c': 43.23768313997447,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 120.00000000000001,
  'volume': 324.0888756946132},
 'sites': [{'species': [{'element': 'Ag', 'occu': 1}],
   'abc': [0.33333333333333315, 0.6666666666666673, 0.36111111111111116],
   'xyz': [-2.1836824750706918e-15, 1.6985366083928457, 15.613607800546339],
   'label': 'Ag',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Ag', 'occu': 1}],
   'abc': [0.0, 3.3306690738754696e-16, 0.4722222222222222],
   'xyz': [-4.899333779093901e-16, 8.485895028629063e-16, 20.417794816099054],
   'label': 'Ag',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Ag', 'occu': 1}],
   'abc': [0.6666666666666665, 0.3333333333333339, 0.41666666666666674],
   'xyz': [1.4709758521260623, 0.8492683041964235, 18.0157013083227],
   'label': 'Ag',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Ag', 'occu': 1}],
   'abc': [0.33333333333333337, 0.6666666666666672, 0.5277777777777778],
   'xyz': [-1.3542375343258014e-15, 1.6985366083928455, 22.819888323875418],
   'label': 'Ag',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Ag', 'occu': 1}],
   'abc': [2.220446049250313e-16, 5.551115123125783e-16, 0.638888888888889],
   'xyz': [-1.6331112596979718e-16, 1.4143158381048438e-15, 27.62407533942814],
   'label': 'Ag',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
  {'species': [{'element': 'Ag', 'occu': 1}],
   'abc': [0.666666666666667, 0.33333333333333404, 0.5833333333333334],
   'xyz': [1.4709758521260634, 0.8492683041964239, 25.221981831651778],
   'label': 'Ag',
   'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}}],
 'oriented_unit_cell': {'@module': 'pymatgen.core.structure',
  '@class': 'Structure',
  'charge': None,
  'lattice': {'matrix': [[2.0802739999999997, -2.080274, -0.0],
    [4.440892098500626e-16, 2.080274, -2.080274],
    [4.160548000000001, 4.160548, 4.160548000000001]],
   'a': 2.9419517042521277,
   'b': 2.941951704252128,
   'c': 7.206280523329078,
   'alpha': 90.0,
   'beta': 90.0,
   'gamma': 120.00000000000001,
   'volume': 54.01481261576888},
  'sites': [{'species': [{'element': 'Ag', 'occu': 1}],
    'abc': [0.33333333333333326, 0.6666666666666667, 0.33333333333333326],
    'xyz': [2.080274, 2.080274, -1.599126756938556e-16],
    'label': 'Ag',
    'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
   {'species': [{'element': 'Ag', 'occu': 1}],
    'abc': [2.220446049250313e-16, 2.220446049250313e-16, 1.0],
    'xyz': [4.160548000000002, 4.160548, 4.160548],
    'label': 'Ag',
    'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}},
   {'species': [{'element': 'Ag', 'occu': 1}],
    'abc': [0.6666666666666665, 0.3333333333333335, 0.6666666666666667],
    'xyz': [4.160548000000001, 2.080274000000001, 2.0802740000000006],
    'label': 'Ag',
    'properties': {'bulk_equivalent': 0, 'bulk_wyckoff': 'a'}}]},
 'miller_index': (1, 1, 1),
 'shift': 0.16666666666666666,
 'scale_factor': [[-1, 1, 0], [-1, 0, 1], [-1, -1, -1]],
 'reconstruction': None,
 'energy': None}
alc111_dict = {'@module': 'pymatgen.core.surface',
 '@class': 'Slab',
 'charge': None,
 'lattice': {'matrix': [[4.953950000000101, 0.0, 3.0334195053180773e-16],
   [2.476975866025646, 4.290247499999913, 3.033420274724064e-16],
   [0.0, 0.0, 32.005944]],
  'a': 4.953950000000101,
  'b': 4.953951256535434,
  'c': 32.005944,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 59.99999682476928,
  'volume': 680.2438231080062},
 'sites': [{'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.666667, 0.666667, 0.387696],
   'xyz': [4.953953054325787, 2.860166430082442, 12.408576465024],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.333333, 0.333333, 0.387696],
   'xyz': [2.47697281169996, 1.430081069917471, 12.408576465024],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.0, 0.0, 0.387696],
   'xyz': [0.0, 0.0, 12.408576465024],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.333333, 0.666667, 0.462678],
   'xyz': [3.302633085025753, 2.860166430082442, 14.808446158032],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [1.0, 0.333334, 0.462678],
   'xyz': [5.779610273325893, 1.4300853601649712, 14.808446158032],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.666666, 0.0, 0.462678],
   'xyz': [3.302630030700067, 0.0, 14.808446158032],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.666667, 0.333333, 0.537322],
   'xyz': [4.1282927809999945, 1.430081069917471, 17.197497841967998],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.333334, 1.0, 0.537322],
   'xyz': [4.128295835325679, 4.290247499999913, 17.197497841967998],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.0, 0.666666, 0.537322],
   'xyz': [1.651315592699853, 2.860162139834942, 17.197497841967998],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.666667, 0.666667, 0.612304],
   'xyz': [4.953953054325787, 2.860166430082442, 19.597367534975998],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.333333, 0.333333, 0.612304],
   'xyz': [2.47697281169996, 1.430081069917471, 19.597367534975998],
   'label': 'Al',
   'properties': {}},
  {'species': [{'element': 'C', 'occu': 1}],
   'abc': [0.0, 0.0, 0.612304],
   'xyz': [0.0, 0.0, 19.597367534975998],
   'label': 'C',
   'properties': {}}],
 'oriented_unit_cell': {'@module': 'pymatgen.core.structure',
  '@class': 'Structure',
  'charge': None,
  'lattice': {'matrix': [[1e-06, 4.95395, 0.0],
    [4.290248, 2.476975, 0.0],
    [-0.0, -0.0, -32.005944]],
   'a': 4.953950000000101,
   'b': 4.953951256535434,
   'c': 32.005944,
   'alpha': 90.0,
   'beta': 90.0,
   'gamma': 59.99999682476928,
   'volume': 680.243823108006},
  'sites': [{'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.666667, 0.666667, 0.387696],
    'xyz': [2.860167430083, 4.953952476975, -12.408576465024],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.333333, 0.333333, 0.387696],
    'xyz': [1.430081569917, 2.4769725230249997, -12.408576465024],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.0, 0.0, 0.387696],
    'xyz': [0.0, 0.0, -12.408576465024],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.333333, 0.666667, 0.462678],
    'xyz': [2.860167096749, 3.302632507675, -14.808446158032],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [1.0, 0.333334, 0.462678],
    'xyz': [1.430086526832, 5.7796099846499995, -14.808446158032],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.666666, 0.0, 0.462678],
    'xyz': [6.666659999999999e-07, 3.3026300306999996, -14.808446158032],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.666667, 0.333333, 0.537322],
    'xyz': [1.430081903251, 4.128292492325, -17.197497841967998],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.333334, 1.0, 0.537322],
    'xyz': [4.2902483333340005, 4.1282949693, -17.197497841967998],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.0, 0.666666, 0.537322],
    'xyz': [2.8601624731679998, 1.6513150153499998, -17.197497841967998],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.666667, 0.666667, 0.612304],
    'xyz': [2.860167430083, 4.953952476975, -19.597367534975998],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.333333, 0.333333, 0.612304],
    'xyz': [1.430081569917, 2.4769725230249997, -19.597367534975998],
    'label': 'Al',
    'properties': {}},
   {'species': [{'element': 'C', 'occu': 1}],
    'abc': [0.0, 0.0, 0.612304],
    'xyz': [0.0, 0.0, -19.597367534975998],
    'label': 'C',
    'properties': {}}]},
 'miller_index': (1, 1, 1),
 'shift': 0,
 'scale_factor': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
 'reconstruction': None,
 'energy': None}

cuau111_dict = {'@module': 'pymatgen.core.surface',
 '@class': 'Slab',
 'charge': None,
 'lattice': {'matrix': [[5.134393814334074, 0.0, 3.143909475143097e-16],
   [2.5671960928329307, 4.446516160791427, 3.143909588830688e-16],
   [0.0, 0.0, 31.288323]],
  'a': 5.134393814334074,
  'b': 5.134394,
  'c': 31.288323,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 60.000011689341854,
  'volume': 714.3175788942775},
 'sites': [{'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.5, 0.5, 0.401016],
   'xyz': [3.8507949535835024, 2.2232580803957136, 12.547118136167999],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.5, 1.0, 0.401016],
   'xyz': [5.134392999999967, 4.446516160791427, 12.547118136167999],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.0, 0.5, 0.401016],
   'xyz': [1.2835980464164654, 2.2232580803957136, 12.547118136167999],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.0, 0.0, 0.401016],
   'xyz': [0.0, 0.0, 12.547118136167999],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.166666, 0.666667, 0.466994],
   'xyz': [2.567193797080454, 2.9643455893663386, 14.611459111062],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.166666, 0.166667, 0.466994],
   'xyz': [1.2835957506639888, 0.7410875089706248, 14.611459111062],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.666666, 0.666667, 0.466994],
   'xyz': [5.134390704247491, 2.9643455893663386, 14.611459111062],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.666666, 0.166667, 0.466994],
   'xyz': [3.850792657831026, 0.7410875089706248, 14.611459111062],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.333334, 0.333333, 0.533006],
   'xyz': [2.5671992029195136, 1.4821705714250888, 16.676863888937998],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.333334, 0.833333, 0.533006],
   'xyz': [3.850797249335979, 3.705428651820802, 16.676863888937998],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.833333, 0.333333, 0.533006],
   'xyz': [5.134390975692736, 1.4821705714250888, 16.676863888937998],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.833333, 0.833333, 0.533006],
   'xyz': [6.417989022109202, 3.705428651820802, 16.676863888937998],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.5, 0.5, 0.598984],
   'xyz': [3.8507949535835024, 2.2232580803957136, 18.741204863832],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.5, 1.0, 0.598984],
   'xyz': [5.134392999999967, 4.446516160791427, 18.741204863832],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.0, 0.5, 0.598984],
   'xyz': [1.2835980464164654, 2.2232580803957136, 18.741204863832],
   'label': 'Cu',
   'properties': {}},
  {'species': [{'element': 'Au', 'occu': 1}],
   'abc': [0.0, 0.0, 0.598984],
   'xyz': [0.0, 0.0, 18.741204863832],
   'label': 'Au',
   'properties': {}}],
 'oriented_unit_cell': {'@module': 'pymatgen.core.structure',
  '@class': 'Structure',
  'charge': None,
  'lattice': {'matrix': [[2.567196, 4.446516, 0.0],
    [5.134394, -0.0, 0.0],
    [-0.0, -0.0, -31.288323]],
   'a': 5.134393814334074,
   'b': 5.134394,
   'c': 31.288323,
   'alpha': 90.0,
   'beta': 90.0,
   'gamma': 60.000011689341854,
   'volume': 714.3175788942776},
  'sites': [{'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.5, 0.5, 0.401016],
    'xyz': [3.850795, 2.223258, -12.547118136167999],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.5, 1.0, 0.401016],
    'xyz': [6.417992, 2.223258, -12.547118136167999],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.0, 0.5, 0.401016],
    'xyz': [2.567197, 0.0, -12.547118136167999],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.0, 0.0, 0.401016],
    'xyz': [0.0, 0.0, -12.547118136167999],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.166666, 0.666667, 0.466994],
    'xyz': [3.850795333334, 0.741083035656, -14.611459111062],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.166666, 0.166667, 0.466994],
    'xyz': [1.2835983333340002, 0.741083035656, -14.611459111062],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.666666, 0.666667, 0.466994],
    'xyz': [5.134393333334001, 2.964341035656, -14.611459111062],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.666666, 0.166667, 0.466994],
    'xyz': [2.567196333334, 2.964341035656, -14.611459111062],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.333334, 0.333333, 0.533006],
    'xyz': [2.567196666666, 1.482174964344, -16.676863888937998],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.333334, 0.833333, 0.533006],
    'xyz': [5.134393666666, 1.482174964344, -16.676863888937998],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.833333, 0.333333, 0.533006],
    'xyz': [3.85079209947, 3.705428517828, -16.676863888937998],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.833333, 0.833333, 0.533006],
    'xyz': [6.41798909947, 3.705428517828, -16.676863888937998],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.5, 0.5, 0.598984],
    'xyz': [3.850795, 2.223258, -18.741204863832],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.5, 1.0, 0.598984],
    'xyz': [6.417992, 2.223258, -18.741204863832],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.0, 0.5, 0.598984],
    'xyz': [2.567197, 0.0, -18.741204863832],
    'label': 'Cu',
    'properties': {}},
   {'species': [{'element': 'Au', 'occu': 1}],
    'abc': [0.0, 0.0, 0.598984],
    'xyz': [0.0, 0.0, -18.741204863832],
    'label': 'Au',
    'properties': {}}]},
 'miller_index': (1, 1, 1),
 'shift': 0,
 'scale_factor': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
 'reconstruction': None,
 'energy': None}
al111_dict = {'@module': 'pymatgen.core.surface',
 '@class': 'Slab',
 'charge': None,
 'lattice': {'matrix': [[2.860165, 0.0, 1.7513459561416446e-16],
   [-1.430082, 2.476975, 1.751345511924877e-16],
   [0.0, 0.0, 32.005944]],
  'a': 2.860165,
  'b': 2.860164274538964,
  'c': 32.005944,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 119.99999682477132,
  'volume': 226.74794103600198},
 'sites': [{'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.0, 0.0, 0.38769624],
   'xyz': [0.0, 0.0, 12.408584146450561],
   'label': 'Al',
   'properties': {'magmom': 0.0}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.333333, 0.666667, 0.46267833],
   'xyz': [-1.0967490001730144e-06, 1.651317492325, 14.80845671999352],
   'label': 'Al',
   'properties': {'magmom': -0.0}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.0, 0.0, 0.61230376],
   'xyz': [0.0, 0.0, 19.597359853549438],
   'label': 'Al',
   'properties': {'magmom': 0.0}},
  {'species': [{'element': 'Al', 'occu': 1}],
   'abc': [0.666667, 0.333333, 0.53732167],
   'xyz': [1.4300840967489998, 0.8256575076749999, 17.197487280006477],
   'label': 'Al',
   'properties': {'magmom': -0.0}}],
 'oriented_unit_cell': {'@module': 'pymatgen.core.structure',
  '@class': 'Structure',
  'charge': None,
  'lattice': {'matrix': [[2.860165, 0.0, 0.0],
    [-1.430082, 2.476975, 0.0],
    [0.0, 0.0, 32.005944]],
   'a': 2.860165,
   'b': 2.860164274538964,
   'c': 32.005944,
   'alpha': 90.0,
   'beta': 90.0,
   'gamma': 119.99999682477132,
   'volume': 226.74794103600198},
  'sites': [{'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.0, 0.0, 0.38769624],
    'xyz': [0.0, 0.0, 12.408584146450561],
    'label': 'Al',
    'properties': {'magmom': 0.0}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.333333, 0.666667, 0.46267833],
    'xyz': [-1.0967490001730144e-06, 1.651317492325, 14.80845671999352],
    'label': 'Al',
    'properties': {'magmom': -0.0}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.0, 0.0, 0.61230376],
    'xyz': [0.0, 0.0, 19.597359853549438],
    'label': 'Al',
    'properties': {'magmom': 0.0}},
   {'species': [{'element': 'Al', 'occu': 1}],
    'abc': [0.666667, 0.333333, 0.53732167],
    'xyz': [1.4300840967489998, 0.8256575076749999, 17.197487280006477],
    'label': 'Al',
    'properties': {'magmom': -0.0}}]},
 'miller_index': (1, 1, 1),
 'shift': 0,
 'scale_factor': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
 'reconstruction': None,
 'energy': None}

cu111_dict = {'@module': 'pymatgen.core.surface',
 '@class': 'Slab',
 'charge': None,
 'lattice': {'matrix': [[2.567197, 0.0, 1.571954794415344e-16],
   [-1.2835990000000004, 2.2232580000000004, 1.5719550437332295e-16],
   [0.0, 0.0, 31.288323]],
  'a': 2.567197,
  'b': 2.567197407167007,
  'c': 31.288323,
  'alpha': 90.0,
  'beta': 90.0,
  'gamma': 120.00000763897559,
  'volume': 178.57939472356944},
 'sites': [{'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.0, -0.0, 0.40101582],
   'xyz': [0.0, 0.0, 12.54711250426986],
   'label': 'Cu',
   'properties': {'magmom': -0.0}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.333333, 0.666667, 0.46699376],
   'xyz': [-1.6169320002745272e-06, 1.4821727410860004, 14.611451601864479],
   'label': 'Cu',
   'properties': {'magmom': -0.0}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.0, -0.0, 0.59898418],
   'xyz': [0.0, 0.0, 18.741210495730137],
   'label': 'Cu',
   'properties': {'magmom': -0.0}},
  {'species': [{'element': 'Cu', 'occu': 1}],
   'abc': [0.666667, 0.333333, 0.53300624],
   'xyz': [1.2835996169319999, 0.7410852589140001, 16.676871398135518],
   'label': 'Cu',
   'properties': {'magmom': -0.0}}],
 'oriented_unit_cell': {'@module': 'pymatgen.core.structure',
  '@class': 'Structure',
  'charge': None,
  'lattice': {'matrix': [[2.567197, 0.0, 0.0],
    [-1.283599, 2.223258, 0.0],
    [0.0, 0.0, 31.288323]],
   'a': 2.567197,
   'b': 2.567197407167006,
   'c': 31.288323,
   'alpha': 90.0,
   'beta': 90.0,
   'gamma': 120.00000763897559,
   'volume': 178.5793947235694},
  'sites': [{'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.0, -0.0, 0.40101582],
    'xyz': [0.0, 0.0, 12.54711250426986],
    'label': 'Cu',
    'properties': {'magmom': -0.0}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.333333, 0.666667, 0.46699376],
    'xyz': [-1.6169319999784675e-06, 1.482172741086, 14.611451601864479],
    'label': 'Cu',
    'properties': {'magmom': -0.0}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.0, -0.0, 0.59898418],
    'xyz': [0.0, 0.0, 18.741210495730137],
    'label': 'Cu',
    'properties': {'magmom': -0.0}},
   {'species': [{'element': 'Cu', 'occu': 1}],
    'abc': [0.666667, 0.333333, 0.53300624],
    'xyz': [1.283599616932, 0.741085258914, 16.676871398135518],
    'label': 'Cu',
    'properties': {'magmom': -0.0}}]},
 'miller_index': (1, 1, 1),
 'shift': 0,
 'scale_factor': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
 'reconstruction': None,
 'energy': None}

slab1 = Slab.from_dict(alc111_dict)
al_bm = 76.50705929207275
slab2 = Slab.from_dict(cu111_dict)
cu_bm = 139.67909156673758
match_params = {'max_area': 200.0,
                'max_mismatch': 0.02,
                'max_angle_diff': 1.5,
                'r1r2_tol': 1,
                'best_match': 'area',
                'interface_distance': 'auto',
                'vacuum': 12
                }

Matcher = InterfaceMatcher(slab_2=slab1,
                         slab_1=slab2,
                         strain_weight_1=1,
                         strain_weight_2=1,
                         **match_params)


#al111_new, cu111_new = Matcher.get_aligned_slabs()
top_aligned, bottom_aligned = Matcher.get_aligned_slabs()
top_centered, bottom_centered = Matcher.get_centered_slabs()
interface = Matcher.get_interface()


# al111_old, cu111_old = get_aligned_lattices(
#                 al111,
#                 cu111,
#                 **match_params)