#include "LaplaceOperator.hpp"

LaplaceOperator::LaplaceOperator(const GridData& pMesh)
    : pMesh_(pMesh)
{}

// ---------------------------------------------------------------------------
// assembleToTriplets
//
// FVM Laplaciano sobre malla P.
// Convención positivo-definida: diagonal = +(a_E+a_W+a_N+a_S),
// fuera-de-diagonal = −a_nb.
// ---------------------------------------------------------------------------
std::vector<Triplet> LaplaceOperator::assembleToTriplets() const
{
    const int n_cells = pMesh_.num_active_nodes;
    std::vector<Triplet> triplets;
    triplets.reserve(5 * n_cells);

    for (int k = 0; k < n_cells; ++k) {
        const InternalCell& cell = pMesh_.cells[k];

        // Solid (immersed obstacle) P-cells are decoupled from the system via
        // a trivial identity row — their pressure value is never used by any
        // real physics (fluid-side cells never reference a solid neighbor,
        // since bc_* is set true at the interface by mark_solid_block).
        if (cell.is_solid) {
            triplets.push_back({ k, k, 1.0 });
            continue;
        }

        const int i   = pMesh_.id_to_i[k];
        const int j   = pMesh_.id_to_j[k];
        const int row = k;

        double a_P = 0.0;

        // Outer-boundary faces are homogeneous Neumann (skip entirely) UNLESS
        // the face is a real OUTLET, in which case a Dirichlet pressure
        // condition (p = bound_*[...].value, typically 0 gauge pressure) is
        // used instead: the face still contributes to the diagonal (a_P +=
        // a_face) but gets NO off-diagonal entry, since the neighbor is a
        // prescribed value, not an unknown. With p_bc = 0 the corresponding
        // RHS term (b[k] -= a_face * p_bc) vanishes, so no RHS change is
        // needed here. Only EAST is exercised in practice (the Square
        // Cylinder's only OUTLET face); all four are handled for generality.
        // LDC/DHC never have an OUTLET-typed boundary, so this is dead code
        // for both (their bc_east/west/north/south faces always take the
        // existing Neumann-skip branch below, byte-identical to before).
        // bound_east/west/north/south entries describe the REAL outer domain
        // edge only (one entry per row/column, built once in generate_mesh());
        // they are meaningless for an interior obstacle interface, so the
        // Dirichlet-outlet check below is gated on actually being at the
        // outer edge (i==N_cells_x etc.), not merely on bc_east being true
        // (which mark_solid_block also sets at a fluid/solid interface).
        const bool at_east_edge  = (i == pMesh_.N_cells_x);
        const bool at_west_edge  = (i == 1);
        const bool at_north_edge = (j == pMesh_.M_cells_y);
        const bool at_south_edge = (j == 1);

        if (!cell.bc_east) {
            const double dist = pMesh_.x[i + 1] - pMesh_.x[i];
            const double a_E  = pMesh_.Se[i][j] / dist;
            triplets.push_back({ row, cell.n_east, -a_E });
            a_P += a_E;
        } else if (at_east_edge && !pMesh_.bound_east.empty() && pMesh_.bound_east[j-1].type == OUTLET) {
            const double dist = pMesh_.x[i + 1] - pMesh_.x[i];
            a_P += pMesh_.Se[i][j] / dist;
        }
        if (!cell.bc_west) {
            const double dist = pMesh_.x[i] - pMesh_.x[i - 1];
            const double a_W  = pMesh_.Sw[i][j] / dist;
            triplets.push_back({ row, cell.n_west, -a_W });
            a_P += a_W;
        } else if (at_west_edge && !pMesh_.bound_west.empty() && pMesh_.bound_west[j-1].type == OUTLET) {
            const double dist = pMesh_.x[i] - pMesh_.x[i - 1];
            a_P += pMesh_.Sw[i][j] / dist;
        }
        if (!cell.bc_north) {
            const double dist = pMesh_.y[j + 1] - pMesh_.y[j];
            const double a_N  = pMesh_.Sn[i][j] / dist;
            triplets.push_back({ row, cell.n_north, -a_N });
            a_P += a_N;
        } else if (at_north_edge && !pMesh_.bound_north.empty() && pMesh_.bound_north[i-1].type == OUTLET) {
            const double dist = pMesh_.y[j + 1] - pMesh_.y[j];
            a_P += pMesh_.Sn[i][j] / dist;
        }
        if (!cell.bc_south) {
            const double dist = pMesh_.y[j] - pMesh_.y[j - 1];
            const double a_S  = pMesh_.Ss[i][j] / dist;
            triplets.push_back({ row, cell.n_south, -a_S });
            a_P += a_S;
        } else if (at_south_edge && !pMesh_.bound_south.empty() && pMesh_.bound_south[i-1].type == OUTLET) {
            const double dist = pMesh_.y[j] - pMesh_.y[j - 1];
            a_P += pMesh_.Ss[i][j] / dist;
        }
        triplets.push_back({ row, row, a_P });
    }
    return triplets;
}

CSRMatrix LaplaceOperator::assemble() const
{
    auto trips = assembleToTriplets();
    return CSRMatrix::fromTriplets(pMesh_.num_active_nodes, trips);
}

Vector LaplaceOperator::apply(const Vector& p) const
{
    return multiply(assemble(), p);
}

void LaplaceOperator::pinReference(std::vector<Triplet>& triplets, Vector& b) const
{
    bool found_diag = false;
    for (auto& t : triplets) {
        if (t.row == 0 && t.col != 0) { t.value = 0.0; continue; }
        if (t.col == 0 && t.row != 0) { t.value = 0.0; continue; }
        if (t.row == 0 && t.col == 0) { t.value = 1.0; found_diag = true; }
    }
    if (!found_diag) triplets.push_back({ 0, 0, 1.0 });
    b[0] = 0.0;
}
