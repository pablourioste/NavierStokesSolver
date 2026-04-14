#include "GradientOperator.hpp"

GradientOperator::GradientOperator(const GridData& pMesh,
                                   const GridData& uMesh,
                                   const GridData& vMesh)
    : pMesh_(pMesh), uMesh_(uMesh), vMesh_(vMesh)
{}

void GradientOperator::apply(const Vector& p,
                              Vector& grad_u,
                              Vector& grad_v) const
{
    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;

    grad_u.assign(uMesh_.num_active_nodes, 0.0);
    grad_v.assign(vMesh_.num_active_nodes, 0.0);

    // ---- ∂p/∂x en nodos U interiores (k=2..N, j=1..M) ----
    // El nodo U (k,j) se encuentra entre las celdas P (k-1,j) y (k,j).
    // La distancia entre centros es pMesh_.x[k] - pMesh_.x[k-1].
    for (int j = 1; j <= M; ++j) {
        for (int k = 2; k <= N; ++k) {
            const double P_east = p[ pIdx(k,     j) ];
            const double P_west = p[ pIdx(k - 1, j) ];
            const double dist   = pMesh_.x[k] - pMesh_.x[k - 1];
            grad_u[ uIdx(k, j) ] = (P_east - P_west) / dist;
        }
    }

    // ---- ∂p/∂y en nodos V interiores (i=1..N, l=2..M) ----
    // El nodo V (i,l) se encuentra entre las celdas P (i,l-1) y (i,l).
    for (int l = 2; l <= M; ++l) {
        for (int i = 1; i <= N; ++i) {
            const double P_north = p[ pIdx(i, l    ) ];
            const double P_south = p[ pIdx(i, l - 1) ];
            const double dist    = pMesh_.y[l] - pMesh_.y[l - 1];
            grad_v[ vIdx(i, l) ] = (P_north - P_south) / dist;
        }
    }
}
