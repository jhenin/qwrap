#include <tcl.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

/*
 * qwrap is a Tcl routine for VMD, with a C++(ish) implementation
 * it is equivalent to "pbc wrap -all -center origin"
 * there are two differences:
 * 1) the center of each "wrapping block" is the reference point, rather than one given atom
 * 2) it's faster (up to 30 times in my tests)
 * 3) some options are hard-coded right now, they're the most likely ones for a 
 *    trajectory from a NAMD biomolecular simulation.
 * 4) it only deals with orthorhombic boxes
 * 5) I can't count!
 *
 * Jerome Henin <jerome.henin@ibpc.fr> March 2013
 */

extern "C" {
int Qwrap_Init(Tcl_Interp *interp);
}
static int obj_qwrap(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[]);

// this is "truncating" a to integer number of times b
// such as a is between -b/2 and +b/2

void truncate(float *a, const std::vector<float> &b)
{
  for (int c=0; c<3; c++)
    a[c] = floor (a[c] / b[c] + 0.5) * b[c];
  return;
}

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

int parse_ivector (Tcl_Obj * const obj, std::vector<int> &vec, Tcl_Interp *interp, bool fromDouble)
{
  Tcl_Obj **data;
  int num;
  double d;

  if (Tcl_ListObjGetElements(interp, obj, &num, &data) != TCL_OK) {
    Tcl_SetResult(interp, (char *) "qwrap: error parsing arguments", TCL_STATIC);
    return -1;
  }

  vec.resize(num);

  if (fromDouble == false) {
    for (int i = 0; i < num; i++) {
      if (Tcl_GetIntFromObj(interp, data[i], &vec[i]) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "qwrap: error parsing vector element as integer", TCL_STATIC);
        return -1;
      }
    }
  } else {
    // do a double-to-int conversion first
    for (int i = 0; i < num; i++) {
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

static int obj_qwrap(ClientData data, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[])
{
  Tcl_Obj *atomselect, *object, *bytes, *centersel;
  int ncoords, result, length, ncenter, nsel;
  int num_frames, first_frame, last_frame;
  int num_atoms;
  enum { NONE, RES, BETA } compound;
  bool refatoms;
  const char *sel_text;

  std::vector<int> blockID;
  std::vector<int> centerID;
  std::vector<int> selID;
  float *coords;
  std::vector<float> PBC;
  std::vector<int> is_ref;

  if (argc % 2 != 1) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"[first <n>] [last <n>] [compound none|res|beta] [refatoms none|occ] [center <seltext>] [sel <seltext>]");
    return TCL_ERROR;
  }

  // Default Values

  compound = RES;
  first_frame = 0;
  last_frame = -1;
  ncenter = 0;
  sel_text = NULL;
  refatoms = false;

  // end default values

  for (int i = 1; i + 1 < argc; i += 2) { 
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

    } else if (!strncmp(cmd, "sel", 4)) {
      sel_text = Tcl_GetString(objv[i+1]);

    } else if (!strncmp(cmd, "center", 4)) {
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
      Tcl_SetResult(interp, (char *) "Usage: qwrap [first <n>] [last <n>] [compound none|res|beta] [center <seltext>]", TCL_STATIC);
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

  for (int frame = first_frame; frame <= last_frame; frame++) {

    if (frame % print == 0) {
      Tcl_Obj *msg = Tcl_ObjPrintf ("puts \"Wrapping frame %i\"", frame);
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

    // ********* PBC *******

    Tcl_Obj *get_abc = Tcl_ObjPrintf ("molinfo top get {a b c}");
    result = Tcl_EvalObjEx(interp, get_abc, TCL_EVAL_DIRECT);
    if (result != TCL_OK) { return TCL_ERROR; }

    object = Tcl_GetObjResult(interp);
    {
    int num = parse_vector(object, PBC, interp); 
      if (num != 3 || PBC[0]*PBC[1]*PBC[2] == 0.0) {
        Tcl_SetResult(interp, (char *) "qwrap: error parsing PBC", TCL_STATIC);
        return TCL_ERROR;
      }
    }

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

    // ******** centering *******
    float shift[3];
    for (int c = 0; c < 3; c++) shift[c] = 0.0;
    if ( ncenter != 0 ) {
      for (int i = 0; i < ncenter; i++) {
        for (int c = 0; c < 3; c++) shift[c] += coords[3 * centerID[i] + c]; 
      }
      for (int c = 0; c < 3; c++) shift[c] /= ncenter;
    }
    
    for (int i = 0; i < num_atoms; i++) {
      for (int c = 0; c < 3; c++) {
        coords[3*i + c] -= shift[c];
      }
    }

    // ******** wrapping *******
    float ref_pos[3];
    int current_block, n_ref, current_atom;

    for (int start_atom = 0; start_atom < ncoords; ) {

      if ( compound != NONE ) { // ref position is the center of ref atoms of the block
        current_block = blockID[start_atom];
        n_ref = 0;  // number of reference atoms within current block
        for (int c = 0; c < 3; c++) ref_pos[c] = 0.0;

        for ( current_atom = start_atom;
              current_atom < ncoords && blockID[current_atom] == current_block;
              current_atom++) {
          
          if (refatoms && is_ref[current_atom] == 0)  // skip non-ref atoms
            continue;

          // accumulating average position of block
          for (int c = 0; c < 3; c++) {
            ref_pos[c] += coords[3*selID[current_atom] + c];
          }
          n_ref++;
        }
        if (n_ref == 0) {
          Tcl_SetResult(interp, (char *) "qwrap: block contains no reference atoms", TCL_STATIC);
          return TCL_ERROR;
        }
        for (int c = 0; c < 3; c++) ref_pos[c] /= n_ref;

      } else {  // ref position is simply the atom position
        current_atom = start_atom;
        for (int c = 0; c < 3; c++) {
          ref_pos[c] = coords[3*selID[current_atom] + c];
        }
        current_atom++;
      }

      truncate (ref_pos, PBC);
      
      for (int i = start_atom; i < current_atom; i++) {
        for (int c = 0; c < 3; c++) coords[3*selID[i] + c] -= ref_pos[c];
      } 

      start_atom = current_atom;
    }
    // ******** wrapping done *******

    // call rawtimestep to set byte array
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
  Tcl_SetResult(interp, (char *) "", TCL_STATIC);
  return TCL_OK;
}

extern "C" {
  int Qwrap_Init(Tcl_Interp *interp) {
    Tcl_CreateObjCommand(interp, "qwrap", obj_qwrap,
                    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_EvalEx(interp, "package provide qwrap 1.0", -1, 0);
    return TCL_OK;
  }
}
