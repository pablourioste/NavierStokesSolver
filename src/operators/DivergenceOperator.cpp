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
        const int i = pMesh_.id_to_i[k];   // índice lógico 1-based
        const int j = pMesh_.id_to_j[k];

        // Velocidades en las cuatro caras de la celda P (i,j)
        //   Cara Este  ↔ nodo U (k_u = i+1, j)
        //   Cara Oeste ↔ nodo U (k_u = i,   j)
        //   Cara Norte ↔ nodo V (i,   l_v = j+1)
        //   Cara Sur   ↔ nodo V (i,   l_v = j  )
        const double u_E = u[ uIdx(i + 1, j    ) ];
        const double u_W = u[ uIdx(i,     j    ) ];
        const double v_N = v[ vIdx(i,     j + 1) ];
        const double v_S = v[ vIdx(i,     j    ) ];

        // Integral FVM del flujo neto
        const double flux =
              u_E * pMesh_.Se[i][j]
            - u_W * pMesh_.Sw[i][j]
            + v_N * pMesh_.Sn[i][j]
            - v_S * pMesh_.Ss[i][j];

        // Dividir por volumen → divergencia por unidad de volumen
        div[k] = flux / pMesh_.V[i][j];
    }
    return div;
}
