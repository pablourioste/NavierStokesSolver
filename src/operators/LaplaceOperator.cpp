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
        const int i   = pMesh_.id_to_i[k];
        const int j   = pMesh_.id_to_j[k];
        const int row = k;

        double a_P = 0.0;

        if (!cell.bc_east) {
            const double dist = pMesh_.x[i + 1] - pMesh_.x[i];
            const double a_E  = pMesh_.Se[i][j] / dist;
            triplets.push_back({ row, cell.n_east, -a_E });
            a_P += a_E;
        }
        if (!cell.bc_west) {
            const double dist = pMesh_.x[i] - pMesh_.x[i - 1];
            const double a_W  = pMesh_.Sw[i][j] / dist;
            triplets.push_back({ row, cell.n_west, -a_W });
            a_P += a_W;
        }
        if (!cell.bc_north) {
            const double dist = pMesh_.y[j + 1] - pMesh_.y[j];
            const double a_N  = pMesh_.Sn[i][j] / dist;
            triplets.push_back({ row, cell.n_north, -a_N });
            a_P += a_N;
        }
        if (!cell.bc_south) {
            const double dist = pMesh_.y[j] - pMesh_.y[j - 1];
            const double a_S  = pMesh_.Ss[i][j] / dist;
            triplets.push_back({ row, cell.n_south, -a_S });
            a_P += a_S;
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
