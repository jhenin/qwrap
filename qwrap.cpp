#include <tcl.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>
#include <limits>
#include <vector>

/*
 * qwrap is a Tcl routine for VMD, with a C++(ish) implementation
 * it is equivalent to "pbc wrap -all -center origin"
 * with several differences:
 * 1) the center of each "wrapping block" is the reference point, rather than one given atom
 * 2) it's faster (up to 30 times in my tests)
 * 3) some options are hard-coded right now, they're the most likely ones for a
 *    trajectory from a NAMD biomolecular simulation.
 * 4) triclinic wrapping is available with the nearest option
 *
 * Jérôme Hénin <jerome.henin@cnrs.fr> 2013-2026
 */

extern "C" {
int Qwrap_Init(Tcl_Interp *interp);
}
static int obj_qwrap(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[]);

// Calculate the shift as an integer multiple of vector b
// such as a is between -b/2 and +b/2

void calc_shift(float *a, const std::vector<float> &b)
{
  for (int c=0; c<3; c++)
    a[c] = floor (a[c] / b[c] + 0.5) * b[c];
  return;
}

struct TriclinicCell {
  double vectors[3][3];
  double reciprocal[3][3];
};

bool build_triclinic_cell(const std::vector<float> &pbc, TriclinicCell &cell)
{
  for (int i = 0; i < 6; i++)
    if (!std::isfinite(pbc[i])) return false;

  const double radians = 3.14159265358979323846 / 180.0;
  double alpha = pbc[3] * radians;
  double beta = pbc[4] * radians;
  double gamma = pbc[5] * radians;
  double cos_alpha = cos(alpha);
  double cos_beta = cos(beta);
  double cos_gamma = cos(gamma);
  double sin_gamma = sin(gamma);

  if (pbc[0] <= 0.0f || pbc[1] <= 0.0f || pbc[2] <= 0.0f || fabs(sin_gamma) < 1e-12)
    return false;

  cell.vectors[0][0] = pbc[0];
  cell.vectors[0][1] = 0.0;
  cell.vectors[0][2] = 0.0;
  cell.vectors[1][0] = pbc[1] * cos_gamma;
  cell.vectors[1][1] = pbc[1] * sin_gamma;
  cell.vectors[1][2] = 0.0;
  cell.vectors[2][0] = pbc[2] * cos_beta;
  cell.vectors[2][1] = pbc[2] * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
  double cz_squared = double(pbc[2]) * pbc[2] -
      cell.vectors[2][0] * cell.vectors[2][0] -
      cell.vectors[2][1] * cell.vectors[2][1];
  if (cz_squared <= 0.0)
    return false;
  cell.vectors[2][2] = sqrt(cz_squared);

  double determinant = cell.vectors[0][0] * cell.vectors[1][1] * cell.vectors[2][2];
  if (determinant <= 0.0)
    return false;

  // Rows of the inverse basis are the reciprocal lattice vectors.
  cell.reciprocal[0][0] = cell.vectors[1][1] * cell.vectors[2][2] / determinant;
  cell.reciprocal[0][1] = -cell.vectors[1][0] * cell.vectors[2][2] / determinant;
  cell.reciprocal[0][2] = (cell.vectors[1][0] * cell.vectors[2][1] -
                           cell.vectors[1][1] * cell.vectors[2][0]) / determinant;
  cell.reciprocal[1][0] = 0.0;
  cell.reciprocal[1][1] = cell.vectors[0][0] * cell.vectors[2][2] / determinant;
  cell.reciprocal[1][2] = -cell.vectors[0][0] * cell.vectors[2][1] / determinant;
  cell.reciprocal[2][0] = 0.0;
  cell.reciprocal[2][1] = 0.0;
  cell.reciprocal[2][2] = cell.vectors[0][0] * cell.vectors[1][1] / determinant;
  return true;
}

bool calc_nearest_lattice_shift(const float position[3], const TriclinicCell &cell, double shift[3])
{
  double fractional[3];
  long long center[3];
  const long double min_index = static_cast<long double>(std::numeric_limits<long long>::min()) + 1.0L;
  const long double max_index = static_cast<long double>(std::numeric_limits<long long>::max()) - 1.0L;
  for (int i = 0; i < 3; i++) {
    fractional[i] = 0.0;
    for (int c = 0; c < 3; c++)
      fractional[i] += cell.reciprocal[i][c] * position[c];
    long double rounded = floor(static_cast<long double>(fractional[i]) + 0.5L);
    if (!std::isfinite(rounded) || rounded < min_index || rounded > max_index)
      return false;
    center[i] = static_cast<long long>(rounded);
  }

  double residual[3] = {position[0], position[1], position[2]};
  for (int i = 0; i < 3; i++)
    for (int c = 0; c < 3; c++)
      residual[c] -= center[i] * cell.vectors[i][c];
  double best_squared = residual[0] * residual[0] + residual[1] * residual[1] + residual[2] * residual[2];
  double best[3] = {double(center[0]), double(center[1]), double(center[2])};
  long long lower[3], upper[3];

  // A better candidate must lie within this reciprocal-space bound in every direction.
  for (int i = 0; i < 3; i++) {
    double reciprocal_norm = 0.0;
    for (int c = 0; c < 3; c++)
      reciprocal_norm += cell.reciprocal[i][c] * cell.reciprocal[i][c];
    long double bound = sqrt(static_cast<long double>(best_squared) * reciprocal_norm);
    long double epsilon = 1e-12L * (1.0L + fabs(fractional[i]) + bound);
    long double low = ceil(static_cast<long double>(fractional[i]) - bound - epsilon);
    long double high = floor(static_cast<long double>(fractional[i]) + bound + epsilon);
    if (!std::isfinite(low) || !std::isfinite(high) || low < min_index ||
      high > max_index || low > high)
      return false;
    lower[i] = static_cast<long long>(low);
    upper[i] = static_cast<long long>(high);
  }

  for (long long i = lower[0];; i++) {
    for (long long j = lower[1];; j++) {
      for (long long k = lower[2];; k++) {
        double candidate[3] = {position[0], position[1], position[2]};
        long long indices[3] = {i, j, k};
        for (int axis = 0; axis < 3; axis++)
          for (int c = 0; c < 3; c++)
            candidate[c] -= indices[axis] * cell.vectors[axis][c];
        double squared = candidate[0] * candidate[0] + candidate[1] * candidate[1] + candidate[2] * candidate[2];
        if (squared < best_squared) {
          best_squared = squared;
          best[0] = i;
          best[1] = j;
          best[2] = k;
        }
        if (k == upper[2]) break;
      }
      if (j == upper[1]) break;
    }
    if (i == upper[0]) break;
  }

  for (int c = 0; c < 3; c++) {
    shift[c] = 0.0;
    for (int axis = 0; axis < 3; axis++)
      shift[c] += best[axis] * cell.vectors[axis][c];
  }
  return true;
}

// This version (for LATtice-based unwrapping) gets the shift as integers
// (in units of PBC vectors)
void calc_shift_int(float *a, const std::vector<float> &b, int shift[3])
{
  for (int c=0; c<3; c++)
    shift[c] = int(floor (a[c] / b[c] + 0.5));
  return;
}

// This version (for TORoidal-view unwrapping) adds the shift to the shift accumulator
void add_shift(float *a, const std::vector<float> &b, double shift[3])
{
  for (int c=0; c<3; c++)
    shift[c] += floor (a[c] / b[c] + 0.5) * b[c];
  return;
}

// Parse a vector of floats from a Tcl object
int parse_vector (Tcl_Obj * const obj, std::vector<float> &vec, Tcl_Interp *interp)
{
  Tcl_Obj **data;
  int num;
  double d;

  if (Tcl_ListObjGetElements(interp, obj, &num, &data) != TCL_OK) {
    Tcl_SetResult(interp, (char *) "qwrap: error parsing arguments", TCL_STATIC);
    return -1;
  }

  vec.resize(num);

  for (int i = 0; i < num; i++) {
    if (Tcl_GetDoubleFromObj(interp, data[i], &d) != TCL_OK) {
      Tcl_SetResult(interp, (char *) "qwrap: error parsing vector element as floating-point", TCL_STATIC);
      return -1;
    }
    // Tcl gives us doubles, make them float
    vec[i] = float (d);
  }
  return num;
}

// Parse a vector of integers from a Tcl object
int parse_ivector (Tcl_Obj * const obj, std::vector<int> &vec, Tcl_Interp *interp, bool fromDouble)
{
  Tcl_Obj **data;
  int num, i;
  double d;

  if (Tcl_ListObjGetElements(interp, obj, &num, &data) != TCL_OK) {
    Tcl_SetResult(interp, (char *) "qwrap: error parsing arguments", TCL_STATIC);
    return -1;
  }

  vec.resize(num);

  if (fromDouble == false) {
    for (i = 0; i < num; i++) {
      if (Tcl_GetIntFromObj(interp, data[i], &vec[i]) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "qwrap: error parsing vector element as integer", TCL_STATIC);
        return -1;
      }
    }
  } else {
    // do a double-to-int conversion first
    for (i = 0; i < num; i++) {
      if (Tcl_GetDoubleFromObj(interp, data[i], &d) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "qwrap: error parsing vector element as integer", TCL_STATIC);
        return -1;
      }
      vec[i] = int (d);
    }
  }
  return num;
}

// ***************************************************************************************
// ***************************************************************************************


static int do_qwrap(ClientData data, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[], bool unwrap, bool TOR)
{
  Tcl_Obj *atomselect, *object, *bytes, *centersel = NULL;
  Tcl_Obj *axiscentersel[3];
  int ncoords, result, length, ncenter, nsel;
  int naxiscenter[3];
  int num_frames, first_frame, last_frame;
  int num_atoms;
  int nearest_value;
  enum { NONE, RES, BETA, FRAGMENT } compound;
  bool refatoms, nearest;
  const char *sel_text;
  int i, c;

  std::vector<int> blockID;
  std::vector<int> centerID;
  std::vector<int> axisCenterID[3];
  std::vector<int> selID;
  std::vector<float> prev_ref_pos;
  std::vector<int> shifts;
  std::vector<double> shifts_TOR;
  std::vector<float> PBC;
  std::vector<int> is_ref;
  float *coords;

  const char* usage =  "Available options:\n"
"  [first <n>] [last <n>]: frame range to process (default: all)\n"
"  [sel <seltext>]: selection of atoms to wrap (default: all)\n"
"  [compound none|res|beta|fragment]: define wrapping blocks\n"
"  [center <seltext>]: center the atoms in <seltext>\n"
"  [centerx <seltext>] [centery <seltext>] [centerz <seltext>]: center the atoms in <seltext> along each axis\n"
"  [nearest yes|no]: wrap to the nearest image (for triclinic cells)";

  if (argc % 2 != 1) {
    Tcl_WrongNumArgs(interp, 1, objv, usage);
    return TCL_ERROR;
  }

  // Default Values

  compound = FRAGMENT; // Correct default for NAMD style wrapping
  first_frame = 0;
  last_frame = -1;
  ncenter = 0;
  for (c = 0; c < 3; c++) {
    axiscentersel[c] = NULL;
    naxiscenter[c] = 0;
  }
  sel_text = NULL;
  refatoms = false;
  nearest = false;

  // end default values

  for (i = 1; i + 1 < argc; i += 2) {
    const char *cmd = Tcl_GetString(objv[i]);
    if (!strncmp(cmd, "first", 4)) {
      if (Tcl_GetIntFromObj(interp, objv[i+1], &first_frame) != TCL_OK) { return TCL_ERROR; }

    } else if (!strncmp(cmd, "last", 4)) {
      if (Tcl_GetIntFromObj(interp, objv[i+1], &last_frame) != TCL_OK) { return TCL_ERROR; }
    } else if (!strncmp(cmd, "compound", 4)) {

      const char *comp = Tcl_GetString(objv[i+1]);
      if (!strncmp(comp, "res", 4)) compound = RES;
      else if (!strncmp(comp, "none", 4)) compound = NONE;
      else if (!strncmp(comp, "beta", 4)) compound = BETA;
      else if (!strncmp(comp, "fragment", 4)) compound = FRAGMENT;
      else {
        Tcl_SetResult(interp, (char *) "qwrap: unknown compound type", TCL_STATIC);
        return TCL_ERROR;
      }

    } else if (!strncmp(cmd, "refatoms", 4)) {

      const char *ref = Tcl_GetString(objv[i+1]);
      if (!strncmp(ref, "occ", 4)) refatoms = true;
      else if (!strncmp(ref, "none", 4)) refatoms = false;
      else {
        Tcl_SetResult(interp, (char *) "qwrap: unknown value for refatoms (accepted: occ, none)", TCL_STATIC);
        return TCL_ERROR;
      }

    } else if (!strncmp(cmd, "nearest", 7)) {
      if (unwrap) {
        Tcl_SetResult(interp, (char *) "qwrap: nearest is only supported by qwrap", TCL_STATIC);
        return TCL_ERROR;
      }
      if (Tcl_GetBooleanFromObj(interp, objv[i+1], &nearest_value) != TCL_OK) { return TCL_ERROR; }
      nearest = nearest_value != 0;

    } else if (!strncmp(cmd, "sel", 4)) {
      sel_text = Tcl_GetString(objv[i+1]);

    } else if (!strcmp(cmd, "centerx") || !strcmp(cmd, "centery") || !strcmp(cmd, "centerz")) {
      if (unwrap) {
        Tcl_SetResult(interp, (char *) "qwrap: qunwrap does not support axis-specific center options", TCL_STATIC);
        return TCL_ERROR;
      }
      int axis = 0;
      if (cmd[6] == 'y') axis = 1;
      else if (cmd[6] == 'z') axis = 2;
      Tcl_Obj *axiscmd = Tcl_ObjPrintf("atomselect top \"%s\"", Tcl_GetString(objv[i+1]));
      result = Tcl_EvalObjEx(interp, axiscmd, TCL_EVAL_DIRECT);
      if (result != TCL_OK) {
        Tcl_SetResult(interp, (char *) "qwrap: error calling atomselect for axis-specific center", TCL_STATIC);
        return TCL_ERROR;
      }
      axiscentersel[axis] = Tcl_GetObjResult(interp);
      Tcl_IncrRefCount(axiscentersel[axis]);
      Tcl_Obj *script = Tcl_DuplicateObj(axiscentersel[axis]);
      Tcl_AppendToObj (script, " get index", -1);
      result = Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT);
      naxiscenter[axis] = parse_ivector(Tcl_GetObjResult(interp), axisCenterID[axis], interp, false );
      if (naxiscenter[axis] < 1) {
        Tcl_SetResult(interp, (char *) "qwrap: no atoms for axis-specific centering", TCL_STATIC);
        return TCL_ERROR;
      }
    } else if (!strncmp(cmd, "center", 4)) {
      if (unwrap) {
        Tcl_SetResult(interp, (char *) "qwrap: qunwrap does not support the center option", TCL_STATIC);
        return TCL_ERROR;
      }
      Tcl_Obj *cmd = Tcl_ObjPrintf("atomselect top \"%s\"", Tcl_GetString(objv[i+1]));
      result = Tcl_EvalObjEx(interp, cmd, TCL_EVAL_DIRECT);
      if (result != TCL_OK) {
        Tcl_SetResult(interp, (char *) "qwrap: error calling atomselect for center", TCL_STATIC);
        return TCL_ERROR;
      }
      centersel = Tcl_GetObjResult(interp);
      Tcl_IncrRefCount(centersel);
      Tcl_Obj *script = Tcl_DuplicateObj(centersel);
      Tcl_AppendToObj (script, " get index", -1);
      result = Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT);
      ncenter = parse_ivector(Tcl_GetObjResult(interp), centerID, interp, false );
      if (ncenter < 1) {
        Tcl_SetResult(interp, (char *) "qwrap: no atoms for centering", TCL_STATIC);
        return TCL_ERROR;
      }
    } else {
      Tcl_SetResult(interp, (char *) usage, TCL_STATIC);
      return TCL_ERROR;
    }
  }

  // Build main selection, which is 'all' unless otherwise specified with the sel option
  if (sel_text == NULL) {
    sel_text = "all";
  }
  Tcl_Obj *cmd =  Tcl_ObjPrintf("atomselect top \"%s\"", sel_text);
  result = Tcl_EvalObjEx(interp, cmd, TCL_EVAL_DIRECT);
  if (result != TCL_OK) {
    Tcl_SetResult(interp, (char *) "qwrap: error calling atomselect for complete selection", TCL_STATIC);
    return TCL_ERROR;
  }
  atomselect = Tcl_GetObjResult(interp);
  Tcl_IncrRefCount(atomselect);   // needed to retain the atomselect object beyond this point!

  // ********* atom IDs for whole selection *******

  Tcl_Obj *script = Tcl_DuplicateObj(atomselect);
  Tcl_AppendToObj (script, " get index", -1);
  result = Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT);
  nsel = parse_ivector(Tcl_GetObjResult(interp), selID, interp, false );
  if (nsel < 1) {
    Tcl_SetResult(interp, (char *) "qwrap: no atoms in selection", TCL_STATIC);
    return TCL_ERROR;
  }

  // ********* block IDs *******

  {
    Tcl_Obj *script = Tcl_DuplicateObj(atomselect);
    if ( compound == RES )
      Tcl_AppendToObj (script, " get residue", -1);
    else if ( compound == BETA )
      Tcl_AppendToObj (script, " get beta", -1);
    else if ( compound == FRAGMENT )
      Tcl_AppendToObj (script, " get fragment", -1);
    else // this case is just to find out how many atoms we have
      Tcl_AppendToObj (script, " get occupancy", -1);
    result = Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT);
    if (result != TCL_OK) {
      Tcl_SetResult(interp, (char *) "qwrap: error calling atomselect", TCL_STATIC);
      return TCL_ERROR;
    }
    ncoords = parse_ivector(Tcl_GetObjResult(interp), blockID, interp, (compound != RES) );
    if (ncoords == -1) {
      Tcl_SetResult(interp, (char *) "qwrap: error parsing atomselect result", TCL_STATIC);
      return TCL_ERROR;
    }

    if (unwrap) {
      // to unwrap, we need to store one set of previous coordinates per block
      // we also set blockIDs to consecutive integers
      int nblocks = 1;

      if (compound == NONE) {
        // We don't have meaningful blockIDs in this case
        nblocks = ncoords;
      } else {
        int current_block = blockID[0];
        blockID[0] = 0;
        for (i = 1; i < ncoords; i++) {
          if (blockID[i] != current_block) {
            current_block = blockID[i];
            nblocks++;
          }
          blockID[i] = nblocks - 1;
        }
      }
      prev_ref_pos.resize(3 * nblocks);
      if (TOR) {
        shifts_TOR.resize(3 * nblocks);
      } else {
        shifts.resize(3 * nblocks);
      }
    }
  }

  // ********* total number of atoms *******

  result = Tcl_EvalEx(interp, "molinfo top get numatoms", -1, 0);
  if (result != TCL_OK) {
    Tcl_SetResult(interp, (char *) "qwrap: error calling molinfo", TCL_STATIC);
    return TCL_ERROR;
  }
  object = Tcl_GetObjResult(interp);
  if (Tcl_GetIntFromObj(interp, object, &num_atoms) != TCL_OK) {
    Tcl_SetResult(interp, (char *) "qwrap: error parsing number of atoms", TCL_STATIC);
    return TCL_ERROR;
  }

  // ********* reference atoms *******

  if (refatoms) {
    Tcl_Obj *script = Tcl_DuplicateObj(atomselect);
    Tcl_AppendToObj (script, " get occupancy", -1);
    result = Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT);
    if (result != TCL_OK) {
      Tcl_SetResult(interp, (char *) "qwrap: error calling atomselect", TCL_STATIC);
      return TCL_ERROR;
    }
    int r = parse_ivector(Tcl_GetObjResult(interp), is_ref, interp, true);
    if (r == -1) {
      Tcl_SetResult(interp, (char *) "qwrap: error parsing atomselect result", TCL_STATIC);
      return TCL_ERROR;
    }
  }

  result = Tcl_EvalEx(interp, "molinfo top get numframes", -1, 0);
  if (result != TCL_OK) {
    Tcl_SetResult(interp, (char *) "qwrap: error calling molinfo", TCL_STATIC);
    return TCL_ERROR;
  }
  object = Tcl_GetObjResult(interp);
  if (Tcl_GetIntFromObj(interp, object, &num_frames) != TCL_OK) {
    Tcl_SetResult(interp, (char *) "qwrap: error parsing number of frames", TCL_STATIC);
    return TCL_ERROR;
  }

  if ( first_frame < 0 || first_frame >= num_frames ) {
    Tcl_SetResult(interp, (char *) "qwrap: illegal value of first_frame", TCL_STATIC);
    return TCL_ERROR;
  }
  if ( last_frame == -1 || last_frame >= num_frames ) last_frame = num_frames - 1;
  int print = ((last_frame - first_frame) / 10);
  if (print < 10) print = 10;
  if (print > 100) print = 100;

  // Loop on frames
  for (int frame = first_frame; frame <= last_frame; frame++) {
    bool orthorhombic = true;
    TriclinicCell triclinic;

    if (frame % print == 0) {
      Tcl_Obj *msg = Tcl_ObjPrintf ("puts \"Frame %i\"", frame);
      result = Tcl_EvalObjEx(interp, msg, TCL_EVAL_DIRECT);
      if (result != TCL_OK) { return TCL_ERROR; }
    }

    Tcl_Obj *chgframe = Tcl_DuplicateObj(atomselect);
    Tcl_AppendPrintfToObj (chgframe, " frame %i", frame);
    result = Tcl_EvalObjEx(interp, chgframe, TCL_EVAL_DIRECT);
    if (result != TCL_OK) { return TCL_ERROR; }

    Tcl_Obj *mol_chgframe = Tcl_ObjPrintf ("molinfo top set frame %i", frame);
    result = Tcl_EvalObjEx(interp, mol_chgframe, TCL_EVAL_DIRECT);
    if (result != TCL_OK) { return TCL_ERROR; }

    // ********* get current PBC *******
    // except if unwrapping and at first frame, then it's not needed

    if (!(unwrap && frame == first_frame)) {
      Tcl_Obj *get_abc = Tcl_ObjPrintf ("molinfo top get {a b c alpha beta gamma}");
      result = Tcl_EvalObjEx(interp, get_abc, TCL_EVAL_DIRECT);
      if (result != TCL_OK) { return TCL_ERROR; }

      object = Tcl_GetObjResult(interp);
      {
      int num = parse_vector(object, PBC, interp);
        if (num != 6) {
          Tcl_SetResult(interp, (char *) "qwrap: error parsing PBC", TCL_STATIC);
          return TCL_ERROR;
        }
        if (PBC[0]*PBC[1]*PBC[2] == 0.0) {
          Tcl_SetResult(interp, (char *) "qwrap: error: at least one PBC box length is zero", TCL_STATIC);
          return TCL_ERROR;
        }
        orthorhombic = PBC[3] == 90. && PBC[4] == 90. && PBC[5] == 90.;
        if (!orthorhombic) {
          if (!nearest || unwrap) {
            Tcl_SetResult(interp, (char *) "qwrap: non-orthorhombic cell detected; use qwrap nearest yes to wrap to the nearest image", TCL_STATIC);
            return TCL_ERROR;
          }
          if (!build_triclinic_cell(PBC, triclinic)) {
            Tcl_SetResult(interp, (char *) "qwrap: invalid non-orthorhombic cell geometry", TCL_STATIC);
            return TCL_ERROR;
          }
        }
      }
    }

    // ********* get current coordinates *******

    Tcl_Obj *get_ts = Tcl_ObjPrintf ("gettimestep %s %i", "top", frame);
    result = Tcl_EvalObjEx(interp, get_ts,  TCL_EVAL_DIRECT);
    if (result != TCL_OK) {
      Tcl_SetResult(interp, (char *) "qwrap: error getting coordinates", TCL_STATIC);
      return TCL_ERROR;
    }

    bytes = Tcl_GetObjResult(interp);
    Tcl_IncrRefCount(bytes);
    Tcl_InvalidateStringRep (bytes);
    coords = reinterpret_cast<float *> (Tcl_GetByteArrayFromObj(bytes, &length));

    if (!unwrap) {
      // ******** centering *******
      if (ncenter != 0 || naxiscenter[0] != 0 || naxiscenter[1] != 0 || naxiscenter[2] != 0) {
        float shift[3] = {0.0f, 0.0f, 0.0f};

        for (c = 0; c < 3; c++) {
          int n;
          const std::vector<int> *ids;

          if (naxiscenter[c] != 0) {
            n   = naxiscenter[c];
            ids = &axisCenterID[c];
          } else if (ncenter != 0) {
            n   = ncenter;
            ids = &centerID;
          } else {
            // No selection defined for this axis: leave it unshifted.
            continue;
          }

          for (i = 0; i < n; i++)
            shift[c] += coords[3 * (*ids)[i] + c];
          shift[c] /= n;
        }

        for (i = 0; i < num_atoms; i++)
          for (c = 0; c < 3; c++)
            coords[3*i + c] -= shift[c];
      }
    }

    // ******** (un)wrapping *******
    float ref_pos[3];
    int current_block, n_ref, current_atom;
    int c;

    for (int start_atom = 0; start_atom < ncoords; ) {

      // First, get the reference position of the current wrapping block
      // (ie the center of its reference atoms)

      if ( compound != NONE ) { // ref position is the center of ref atoms of the block
        current_block = blockID[start_atom];
        n_ref = 0;  // number of reference atoms within current block
        for (c = 0; c < 3; c++) ref_pos[c] = 0.0;

        for ( current_atom = start_atom;
              current_atom < ncoords && blockID[current_atom] == current_block;
              current_atom++) {

          if (refatoms && is_ref[current_atom] == 0)  // skip non-ref atoms
            continue;

          for (c = 0; c < 3; c++) {
            ref_pos[c] += coords[3*selID[current_atom] + c];
          }
          n_ref++;
        }
        if (n_ref == 0) {
          Tcl_SetResult(interp, (char *) "qwrap: block contains no reference atoms", TCL_STATIC);
          return TCL_ERROR;
        }
        for (c = 0; c < 3; c++) ref_pos[c] /= n_ref;

      } else {  // compound == NONE: ref position is simply the atom position
        current_atom = start_atom;
        current_block = start_atom; // no need for blockID array
        for (c = 0; c < 3; c++) {
          ref_pos[c] = coords[3*selID[current_atom] + c];
        }
        current_atom++; // Will be the next start_atom
      }

      if (unwrap) {
        // On all frames, record the reference position
        for (c = 0; c < 3; c++) {
          float tmp = ref_pos[c]; // remember ref position
          ref_pos[c] -= prev_ref_pos[current_block * 3 + c]; // new refpos is the displacement from previous one
          prev_ref_pos[current_block * 3 + c] = tmp;  // save the refpos for next frame
        }

        // On frames after the first one, unwrap
        if (frame != first_frame) {

          if (TOR) {
            // ************ TORoidal-view Unwrapping **************
            // see Bullerjahn et al. for a comparison
            // https://doi.org/10.1021/acs.jctc.3c00308

            // Get the shift needed to unwrap the reference position, increment the floating-point shift
            add_shift(ref_pos, PBC, &(shifts_TOR[current_block * 3]));

            // Shift all atoms within the unwrapping block
            for (i = start_atom; i < current_atom; i++) {
              for (c = 0; c < 3; c++) coords[3*selID[i] + c] -= shifts_TOR[current_block * 3 + c];
            }

          } else {
            // ************ LATt-ice-based Unwrapping **************
            // see Bullerjahn et al. for a comparison
            // https://doi.org/10.1021/acs.jctc.3c00308
            // Note on the unwrapping algorithm:
            // We compare the current ref_pos to the *unmodified* reference position from the previous frame
            // so we detect PBC jumps one at a time and accumulate them in the shifts array
            // Other algorithms just wrap the current frame around the already unwrapped previous frame
            // we don't do it so we don't have to recompute the reference position after unwrapping atom-wise
            // this would be trivial for compound none, but not in other cases
            // TODO: for unwrapping, ref_pos should not be the average, but the position of a single reference atom
            // compounds should be kept whole by wrapping relative displacements within them

            int shift[3];
            // Get the shift needed to unwrap the reference position displacement, increment the (integer) shift counter
            calc_shift_int(ref_pos, PBC, shift);
            for (c = 0; c < 3; c++) {
              shifts[current_block * 3 + c] += shift[c];
            }

            // Shift all atoms within the unwrapping block
            for (i = start_atom; i < current_atom; i++) {
              for (c = 0; c < 3; c++) coords[3*selID[i] + c] -= shifts[current_block * 3 + c] * PBC[c];
            }
          }
        }
      } else {

        // Get the shift needed to wrap the reference position
        if (orthorhombic) {
          calc_shift(ref_pos, PBC);
        } else {
          double shift[3];
          if (!calc_nearest_lattice_shift(ref_pos, triclinic, shift)) {
            Tcl_SetResult(interp, (char *) "qwrap: could not calculate nearest periodic image", TCL_STATIC);
            return TCL_ERROR;
          }
          for (c = 0; c < 3; c++) ref_pos[c] = float(shift[c]);
        }

        // Actually shift all atoms within the wrapping block
        for (i = start_atom; i < current_atom; i++) {
          for (c = 0; c < 3; c++) coords[3*selID[i] + c] -= ref_pos[c];
        }
      }

      // Next wrapping block starts here
      start_atom = current_atom;
    }
    // ******** wrapping done *******

    // call rawtimestep to set VMD coords from byte array
    Tcl_Obj *set_ts[5];

    set_ts[0] = Tcl_NewStringObj("rawtimestep", -1);
    set_ts[1] = Tcl_NewStringObj("top", -1);
    set_ts[2] = bytes;
    set_ts[3] = Tcl_NewStringObj("-frame", -1);
    set_ts[4] = Tcl_NewIntObj(frame);

    result = Tcl_EvalObjv (interp, 5, set_ts, 0);
    if (result != TCL_OK) { return TCL_ERROR; }
    Tcl_DecrRefCount(bytes);
  } // end loop on frames

  Tcl_DecrRefCount(atomselect);
  if (ncenter != 0) Tcl_DecrRefCount(centersel);
  for (c = 0; c < 3; c++) {
    if (naxiscenter[c] != 0) Tcl_DecrRefCount(axiscentersel[c]);
  }
  Tcl_SetResult(interp, (char *) "", TCL_STATIC);
  return TCL_OK;
}


// "Wrapper" functions

static int obj_qwrap(ClientData data, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[])
{
    return do_qwrap(data, interp, argc, objv, false, false);
}

static int obj_qunwrap(ClientData data, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[])
{
    return do_qwrap(data, interp, argc, objv, true, false);
}

static int obj_qunwrap_tor(ClientData data, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[])
{
    return do_qwrap(data, interp, argc, objv, true, true);
}

extern "C" {
  int Qwrap_Init(Tcl_Interp *interp) {

#ifdef USE_TCL_STUBS
    const char *version = Tcl_InitStubs(interp, "8.5", 0);
#else
    const char *version = Tcl_PkgRequire(interp, "Tcl", "8.5", 0);
#endif

    if (version == NULL) {
      return TCL_ERROR;
    }

    Tcl_PkgProvide(interp, "qwrap", VERSION);
    Tcl_CreateObjCommand(interp, "qwrap", obj_qwrap,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "qunwrap", obj_qunwrap,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateObjCommand(interp, "qunwrap_tor", obj_qunwrap_tor,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    return TCL_OK;
  }
}
