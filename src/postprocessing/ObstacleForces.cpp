#include "ObstacleForces.hpp"
#include "mesh/ObstacleMask.hpp"
#include <cmath>

ObstacleForces::ObstacleForces(const GridData& pMesh, const GridData& uMesh,
                               const GridData& vMesh, double nu, double rho)
    : pMesh_(pMesh), uMesh_(uMesh), vMesh_(vMesh), nu_(nu), rho_(rho)
{}

std::pair<double, double> ObstacleForces::computeForces(const Vector& p,
                                                          const Vector& u,
                                                          const Vector& v) const
{
    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;
    double Fx = 0.0, Fy = 0.0;

    for (int j = 1; j <= M; ++j) {
        for (int i = 1; i <= N; ++i) {
            if (pMesh_.cells[pIdx(i, j)].is_solid) continue;

            // East face of fluid cell (i,j) touches the obstacle's west face
            // when cell (i+1,j) is solid: outward normal (from the fluid's
            // perspective, pushing on the obstacle) is +x.
            if (i < N && pMesh_.cells[pIdx(i + 1, j)].is_solid) {
                const double p_f   = p[pIdx(i, j)];
                const double area  = pMesh_.Se[i][j];
                Fx += p_f * area;  // pressure pushes +x on the obstacle's west face

                // Viscous shear: tangential velocity is v here (face normal
                // is x), no-slip at the obstacle surface -> one-sided diff
                // between the fluid-side V-node and the (zero) wall value,
                // same ghost distance as DiffusionOperator's west/east branch.
                const double v_S = v[vIdx(i, j)];
                const double v_N = v[vIdx(i, j + 1)];
                const double v_avg = 0.5 * (v_S + v_N);
                const double dist = pMesh_.obstacle_xmin - pMesh_.x[i];
                const double dv_dn = (0.0 - v_avg) / dist;
                Fx += -rho_ * nu_ * dv_dn * area;
            }
            // West face of fluid cell (i,j) touches the obstacle's east face.
            if (i > 1 && pMesh_.cells[pIdx(i - 1, j)].is_solid) {
                const double p_f  = p[pIdx(i, j)];
                const double area = pMesh_.Sw[i][j];
                Fx -= p_f * area;  // pressure pushes -x on the obstacle's east face

                const double v_S = v[vIdx(i, j)];
                const double v_N = v[vIdx(i, j + 1)];
                const double v_avg = 0.5 * (v_S + v_N);
                const double dist = pMesh_.x[i] - pMesh_.obstacle_xmax;
                const double dv_dn = (v_avg - 0.0) / dist;
                Fx += rho_ * nu_ * dv_dn * area;
            }
            // North face of fluid cell (i,j) touches the obstacle's south face.
            if (j < M && pMesh_.cells[pIdx(i, j + 1)].is_solid) {
                const double p_f  = p[pIdx(i, j)];
                const double area = pMesh_.Sn[i][j];
                Fy += p_f * area;

                const double u_W = u[uIdx(i, j)];
                const double u_E = u[uIdx(i + 1, j)];
                const double u_avg = 0.5 * (u_W + u_E);
                const double dist = pMesh_.obstacle_ymin - pMesh_.y[j];
                const double du_dn = (0.0 - u_avg) / dist;
                Fy += -rho_ * nu_ * du_dn * area;
            }
            // South face of fluid cell (i,j) touches the obstacle's north face.
            if (j > 1 && pMesh_.cells[pIdx(i, j - 1)].is_solid) {
                const double p_f  = p[pIdx(i, j)];
                const double area = pMesh_.Ss[i][j];
                Fy -= p_f * area;

                const double u_W = u[uIdx(i, j)];
                const double u_E = u[uIdx(i + 1, j)];
                const double u_avg = 0.5 * (u_W + u_E);
                const double dist = pMesh_.y[j] - pMesh_.obstacle_ymax;
                const double du_dn = (u_avg - 0.0) / dist;
                Fy += rho_ * nu_ * du_dn * area;
            }
        }
    }

    return { Fx, Fy };
}

double estimate_shedding_frequency(const std::vector<double>& time,
                                    const std::vector<double>& Cl)
{
    const std::size_t n = Cl.size();
    if (n < 4) return 0.0;

    // Discard the first half as the startup transient.
    const std::size_t start = n / 2;

    // Mean-subtract over the retained window so crossings are centered.
    double mean = 0.0;
    for (std::size_t k = start; k < n; ++k) mean += Cl[k];
    mean /= static_cast<double>(n - start);

    std::vector<double> crossing_times;
    for (std::size_t k = start + 1; k < n; ++k) {
        const double prev = Cl[k - 1] - mean;
        const double curr = Cl[k]     - mean;
        if (prev < 0.0 && curr >= 0.0) {
            // Linear interpolation for the zero-crossing time.
            const double frac = (curr == prev) ? 0.0 : (-prev) / (curr - prev);
            crossing_times.push_back(time[k - 1] + frac * (time[k] - time[k - 1]));
        }
    }

    if (crossing_times.size() < 2) return 0.0;

    double total_period = 0.0;
    for (std::size_t k = 1; k < crossing_times.size(); ++k)
        total_period += crossing_times[k] - crossing_times[k - 1];
    const double avg_period = total_period / static_cast<double>(crossing_times.size() - 1);

    return (avg_period > 0.0) ? 1.0 / avg_period : 0.0;
}
