"""
Parse the true MuSIC_v2 detector geometry (dd4hep/lcgeo compact XML,
key4hep/k4geo repo: MuColl/MuSIC/compact/MuSIC_v2/) to compute the real
sensitive area of each (system, side, layer) region, replacing the
hit-inferred area estimate used for the mean-density calculation in
bib_common.region_table().

Source files (as included by MuSIC_v2.xml):
    MuSIC_v2.xml               - top-level constants (DetID_* mapping)
    Vertex_o2_v06_01.xml       - VXD barrel + endcap
    InnerTracker_o2_v07_01.xml - IT barrel + endcap
    OuterTracker_o2_v07_01.xml - OT barrel + endcap

Method:
  - Barrel layers (VXD/IT/OT): exact. VXD ladders are single silicon
    pieces (nLadders x width x length). IT/OT barrels tile each layer
    with square "module_envelope" tiles on an (nphi x nz) grid
    (<rphi_layout nphi=.../><z_layout nz=.../>) - area = nphi*nz*tile_area.
  - Endcap disks (VXD/IT/OT): each disk is built from 1+ concentric
    rings of trapezoidal wedge modules (<ring r=... nmodules=... module=.../>),
    each ring's module type giving its trapezoid dimensions
    (<trd x1= x2= z=/> for VXD, <trd x= y=/> approximated as a
    rectangle for IT/OT). Area = sum over rings of nmodules * module_area.
    This is the true tiled area (not an approximation).

A small recursive constant resolver handles the fact that XML constants
reference each other symbolically (e.g. "Vertex_outer_radius - VertexEndcap_offset")
and carry unit suffixes (mm, cm, deg, rad) - it is NOT a general dd4hep
XML/geometry-driver reimplementation, just enough arithmetic to resolve
the specific constants and shape formulas used in this geometry file set.
"""
import math
import re
from pathlib import Path

_UNIT_MULT = {
    "mm": 1.0, "cm": 10.0, "m": 1000.0,
    "deg": math.pi / 180.0, "rad": 1.0,
    "tesla": 1.0, "ns": 1.0,
}

_SAFE_NAMES = {
    "pi": math.pi, "tan": math.tan, "sin": math.sin, "cos": math.cos,
    "sqrt": math.sqrt, "int": int,
}


def _strip_comments(xml_text):
    """Remove <!-- ... --> blocks so commented-out (e.g. historically
    disabled double-layer) geometry elements aren't picked up by the
    regex-based extraction below, which isn't comment-aware."""
    return re.sub(r"<!--.*?-->", "", xml_text, flags=re.S)


def _parse_constants(xml_text):
    """Return {name: raw_expression_string} for every <constant> tag."""
    out = {}
    for m in re.finditer(
        r'<constant\s+name="([^"]+)"[^>]*\bvalue="([^"]+)"', xml_text
    ):
        out[m.group(1)] = m.group(2)
    return out


def _strip_units(expr):
    # "30*mm" -> "30*1.0", "45*deg" -> "45*(pi/180.0)" handled via _UNIT_MULT
    def repl(m):
        return f"({_UNIT_MULT[m.group(1)]})"
    return re.sub(r"\*(mm|cm|m|deg|rad|tesla|ns)\b", lambda m: f"*{_UNIT_MULT[m.group(1)]}", expr)


class GeometryResolver:
    """Resolves symbolic constant expressions to floats, on demand, with
    memoization. Also resolves one-off formula strings (e.g. shape
    dimensions) via `eval_expr()`."""

    def __init__(self):
        self.raw = {}      # name -> raw expression string
        self._cache = {}    # name -> resolved float

    def load(self, xml_text):
        self.raw.update(_parse_constants(xml_text))

    def resolve(self, name):
        if name in self._cache:
            return self._cache[name]
        if name not in self.raw:
            raise KeyError(f"Unknown constant: {name}")
        val = self.eval_expr(self.raw[name])
        self._cache[name] = val
        return val

    def eval_expr(self, expr):
        expr = _strip_units(expr)
        # substitute known constant names, longest first to avoid
        # partial-name collisions (e.g. "Vertex_outer_radius" inside
        # "VertexEndcap_something")
        names = sorted(self.raw.keys(), key=len, reverse=True)
        pattern = re.compile(r"\b(" + "|".join(re.escape(n) for n in names) + r")\b")

        def sub(m):
            return f"({self.resolve(m.group(1))})"

        # iterate substitution until stable (handles nested references)
        prev = None
        cur = expr
        for _ in range(20):
            if cur == prev:
                break
            prev = cur
            cur = pattern.sub(sub, cur)
        return eval(cur, {"__builtins__": {}}, _SAFE_NAMES)


def _find_modules_trd_xy(xml_text):
    """module name=... with <trd x=.. y=..> (rectangle approx) -> {name: (x,y)}"""
    out = {}
    for m in re.finditer(
        r'<module\s+name="([^"]+)"[^>]*>\s*<trd\s+x="([^"]+)"\s+y="([^"]+)"',
        xml_text,
    ):
        out[m.group(1)] = (m.group(2), m.group(3))
    return out


def _find_modules_trd_x1x2z(xml_text):
    """module name=... with <trd x1=.. x2=.. z=..> (trapezoid) -> {name: (x1,x2,z)}"""
    out = {}
    for m in re.finditer(
        r'<module\s+name="([^"]+)"[^>]*>\s*<trd\s+x1="([^"]+)"\s+x2="([^"]+)"\s+z="([^"]+)"',
        xml_text,
    ):
        out[m.group(1)] = (m.group(2), m.group(3), m.group(4))
    return out


def _find_rings(xml_text):
    """<layer id=K> ... <ring r=.. nmodules=.. module="name"/> ... </layer>
    -> {layer_id: [(r_expr, nmodules_expr, module_name), ...]}"""
    out = {}
    for lm in re.finditer(r'<layer\s+id="(\d+)"[^>]*>(.*?)</layer>', xml_text, re.S):
        layer_id = int(lm.group(1))
        body = lm.group(2)
        rings = []
        for rm in re.finditer(
            r'<ring\s+r="([^"]+)"[^>]*\bnmodules="([^"]+)"[^>]*\bmodule="([^"]+)"',
            body,
        ):
            rings.append((rm.group(1), rm.group(2), rm.group(3)))
        if rings:
            out[layer_id] = rings
    return out


def _find_barrel_layers(xml_text, layer_id_re=r'<layer\s+(?:module="[^"]*"\s+)?id="(\d+)"[^>]*>(.*?)</layer>'):
    """<layer id=K [module=..]> <rphi_layout nphi=.. rc=../> <z_layout nz=../> </layer>
    -> {layer_id: (nphi_expr, nz_expr)} ; nz missing (VXD) -> nz_expr=None"""
    out = {}
    for lm in re.finditer(layer_id_re, xml_text, re.S):
        layer_id = int(lm.group(1))
        body = lm.group(2)
        m_phi = re.search(r'<rphi_layout\s+[^>]*\bnphi="([^"]+)"', body)
        m_z = re.search(r'<z_layout\s+[^>]*\bnz="([^"]+)"', body)
        if m_phi:
            out[layer_id] = (m_phi.group(1), m_z.group(1) if m_z else None)
    return out


def build_area_lookup(geom_dir):
    """
    Parse the 4 geometry XML files under `geom_dir` and return
    {(system, side, layer): area_mm2} for every region we can compute.

    side convention matches bib_common's decode_id0(): 0 = barrel,
    1/3 = endcap +z/-z (endcap area is the same for both sides, by
    symmetry, so both are filled with the same value).
    """
    geom_dir = Path(geom_dir)
    main_xml = _strip_comments((geom_dir / "MuSIC_v2.xml").read_text())
    vtx_xml = _strip_comments((geom_dir / "Vertex_o2_v06_01.xml").read_text())
    it_xml = _strip_comments((geom_dir / "InnerTracker_o2_v07_01.xml").read_text())
    ot_xml = _strip_comments((geom_dir / "OuterTracker_o2_v07_01.xml").read_text())

    R = GeometryResolver()
    for txt in (main_xml, vtx_xml, it_xml, ot_xml):
        R.load(txt)

    area = {}

    # ---------------- VXD barrel (system 1): single-piece ladders --------
    for m in re.finditer(
        r'<layer\s+nLadders="([^"]+)"\s+phi0="[^"]*"\s+id="(\d+)"\s*>(.*?)</layer>',
        vtx_xml, re.S,
    ):
        n_ladders_expr, layer_id, body = m.group(1), int(m.group(2)), m.group(3)
        wm = re.search(r'width="([^"]+)"', body)
        lm = re.search(r'length="([^"]+)"', body)
        if not (wm and lm):
            continue
        n_ladders = R.eval_expr(n_ladders_expr)
        width = R.eval_expr(wm.group(1))
        length = R.eval_expr(lm.group(1))
        area[(1, 0, layer_id)] = n_ladders * width * length

    # ---------------- VXD endcap (system 2): ring of trapezoid modules ----
    vxd_modules = _find_modules_trd_x1x2z(vtx_xml)
    vxd_rings = _find_rings(vtx_xml)
    for layer_id, rings in vxd_rings.items():
        total = 0.0
        for r_expr, nmod_expr, mod_name in rings:
            if mod_name not in vxd_modules:
                continue
            x1e, x2e, ze = vxd_modules[mod_name]
            x1, x2, zlen = R.eval_expr(x1e), R.eval_expr(x2e), R.eval_expr(ze)
            nmod = R.eval_expr(nmod_expr)
            # G4Trd/dd4hep <trd x1 x2 z> convention: x1, x2, z are all
            # HALF-widths/half-length (like G4Trd's dx1,dx2,dz) -> full
            # parallel widths = 2*x1, 2*x2; full radial extent = 2*z.
            # Trapezoid area = average_full_width * full_height.
            avg_full_width = x1 + x2  # = (2*x1 + 2*x2) / 2
            full_height = 2.0 * zlen
            total += nmod * avg_full_width * full_height
        for side in (1, 3):
            area[(2, side, layer_id)] = total

    # ---------------- IT / OT barrels (systems 3, 5): tiled square modules
    for system, xml_text, tile_key in ((3, it_xml, None), (5, ot_xml, None)):
        # module_envelope width/length (assume same for all layers in file)
        me = re.search(
            r'<module_envelope\s+width="([^"]+)"\s+length="([^"]+)"', xml_text
        )
        tile_w = R.eval_expr(me.group(1))
        tile_l = R.eval_expr(me.group(2))
        tile_area = tile_w * tile_l
        layers = _find_barrel_layers(xml_text)
        for layer_id, (nphi_e, nz_e) in layers.items():
            nphi = R.eval_expr(nphi_e)
            nz = R.eval_expr(nz_e) if nz_e else 1
            area[(system, 0, layer_id)] = nphi * nz * tile_area

    # ---------------- IT / OT endcaps (systems 4, 6): ring of rect modules
    for system, xml_text in ((4, it_xml), (6, ot_xml)):
        modules = _find_modules_trd_xy(xml_text)
        rings = _find_rings(xml_text)
        for layer_id, ring_list in rings.items():
            total = 0.0
            for r_expr, nmod_expr, mod_name in ring_list:
                if mod_name not in modules:
                    continue
                xe, ye = modules[mod_name]
                mx, my = R.eval_expr(xe), R.eval_expr(ye)
                nmod = R.eval_expr(nmod_expr)
                total += nmod * mx * my
            for side in (1, 3):
                area[(system, side, layer_id)] = total

    return area


def annotate_rows_with_geometry_area(rows, area_lookup):
    """
    Add 'area_mm2_geometry' and recompute 'density_hits_per_mm2' from it
    for every row whose (system, side, layer) is found in area_lookup.
    The original hit-inferred area/density are kept as
    'area_mm2_hit_inferred' / 'density_hits_per_mm2_hit_inferred' for
    comparison. Rows with no geometry match keep the hit-inferred values
    as their primary density (and area_mm2_geometry is left absent).
    """
    for row in rows:
        key = (row["system"], row["side"], row["layer"])
        row["area_mm2_hit_inferred"] = row["area_mm2"]
        row["density_hits_per_mm2_hit_inferred"] = row["density_hits_per_mm2"]
        if key in area_lookup:
            geo_area = area_lookup[key]
            row["area_mm2_geometry"] = geo_area
            row["area_mm2"] = geo_area
            row["density_hits_per_mm2"] = (
                row["n_hits"] / geo_area if geo_area > 0 else float("nan")
            )
    return rows
