#pragma once
#include "Mesh_config.h"
#include <vector>

// =============================================================================
// ObstacleMask — shared helpers for immersed-obstacle-aware stencils.
//
// Used by MomentumSolver and DiffusionOperator to decide, at every U/V
// neighbor lookup, whether that neighbor sits on the outer domain edge, is
// blocked by a solid obstacle, or is a legitimate fluid node — so the same
// ghost-image method already used for outer walls extends naturally to an
// interior obstacle surface. Vacuous (always false) for scenarios that never
// call MeshConfig::mark_solid_block() (Lid-Driven Cavity, DHC).
// =============================================================================

// P-cell (i,j) is solid if marked so by mark_solid_block(); out-of-range
// indices are never queried here (callers only pass valid P-cell indices).
inline bool isSolidP(const GridData& pMesh, int i, int j) {
    const int idx = (j - 1) * pMesh.N_cells_x + (i - 1);
    return pMesh.cells[idx].is_solid;
}

// U-node (k,j), interior range k=2..N, straddles P-cells (k-1,j) and (k,j).
// "Blocked" = strictly inside the solid block (both straddling P-cells solid)
// — such a node is never a legitimate flow unknown, its stencil is skipped.
inline bool uNodeBlocked(const GridData& pMesh, int k, int j) {
    if (k < 2 || k > pMesh.N_cells_x) return false; // boundary nodes handled elsewhere
    return isSolidP(pMesh, k - 1, j) && isSolidP(pMesh, k, j);
}

// V-node (i,l), interior range l=2..M, straddles P-cells (i,l-1) and (i,l).
inline bool vNodeBlocked(const GridData& pMesh, int i, int l) {
    if (l < 2 || l > pMesh.M_cells_y) return false;
    return isSolidP(pMesh, i, l - 1) && isSolidP(pMesh, i, l);
}

// True if the U-node (k,j) sits exactly on the obstacle surface (normal
// velocity node): exactly one of the two straddling P-cells is solid.
inline bool uNodeOnObstacleSurface(const GridData& pMesh, int k, int j) {
    if (k < 2 || k > pMesh.N_cells_x) return false;
    return isSolidP(pMesh, k - 1, j) != isSolidP(pMesh, k, j);
}

// True if the V-node (i,l) sits exactly on the obstacle surface.
inline bool vNodeOnObstacleSurface(const GridData& pMesh, int i, int l) {
    if (l < 2 || l > pMesh.M_cells_y) return false;
    return isSolidP(pMesh, i, l - 1) != isSolidP(pMesh, i, l);
}
