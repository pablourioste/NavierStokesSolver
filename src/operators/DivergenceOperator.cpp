#include "DivergenceOperator.hpp"

DivergenceOperator::DivergenceOperator(const GridData& pMesh,
                                       const GridData& uMesh,
                                       const GridData& vMesh)
    : pMesh_(pMesh), uMesh_(uMesh), vMesh_(vMesh)
{}

Vector DivergenceOperator::apply(const Vector& u, const Vector& v) const
{
    const int n_cells = pMesh_.num_active_nodes;
    Vector div(n_cells, 0.0);

    for (int k = 0; k < n_cells; ++k) {
        // Solid (immersed obstacle) P-cells have no fluid volume — their
        // divergence is meaningless and must not be computed from the
        // surrounding no-penetration-fixed U/V nodes.
        if (pMesh_.cells[k].is_solid) { div[k] = 0.0; continue; }

        const int i = pMesh_.id_to_i[k];   // 1-based logical index
        const int j = pMesh_.id_to_j[k];

        // Velocities on the four faces of P cell (i,j)
        //   East face  ↔ U node (k_u = i+1, j)
        //   West face  ↔ U node (k_u = i,   j)
        //   North face ↔ V node (i,   l_v = j+1)
        //   South face ↔ V node (i,   l_v = j  )
        const double u_E = u[ uIdx(i + 1, j    ) ];
        const double u_W = u[ uIdx(i,     j    ) ];
        const double v_N = v[ vIdx(i,     j + 1) ];
        const double v_S = v[ vIdx(i,     j    ) ];

        // FVM integral of the net flux
        const double flux =
              u_E * pMesh_.Se[i][j]
            - u_W * pMesh_.Sw[i][j]
            + v_N * pMesh_.Sn[i][j]
            - v_S * pMesh_.Ss[i][j];

        // Divide by volume → divergence per unit volume
        div[k] = flux / pMesh_.V[i][j];
    }
    return div;
}
