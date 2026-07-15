#include "custom_elements/suspended_particle.h"
#include "DEM_application_variables.h"

namespace Kratos
{

void SuspendedParticle::ComputeAdditionalForces(
    array_1d<double, 3>& additionally_applied_force,
    array_1d<double, 3>& additionally_applied_moment,
    const ProcessInfo& r_current_process_info,
    const array_1d<double,3>& gravity)
{
    KRATOS_TRY

    SphericParticle::ComputeAdditionalForces(
        additionally_applied_force,
        additionally_applied_moment,
        r_current_process_info,
        gravity);

    const Properties& props = GetProperties();
    if (!props.Has(COUPLING_KERNEL_LENGTH) || !props.Has(COUPLING_SEARCH_RADIUS) || !props.Has(DPD_FLUID_VISCOSITY)) {
        return;
    }

    const double h  = props[COUPLING_KERNEL_LENGTH];
    const double rc = props[COUPLING_SEARCH_RADIUS];
    const double mu = props[DPD_FLUID_VISCOSITY];

    if (h <= 0.0 || rc <= 0.0 || mu <= 0.0) {
        return;
    }

    const auto& my_pos = GetGeometry()[0].Coordinates();
    const auto& my_vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

    array_1d<double, 3> uf = ZeroVector(3);
    double wsum = 0.0;

    for (SphericParticle* p_neighbour : mNeighbourElements)
    {
        if (!p_neighbour) continue;
        if (p_neighbour->Info() != "DPDParticle") continue;

        const auto& xj = p_neighbour->GetGeometry()[0].Coordinates();
        const double dx = my_pos[0] - xj[0];
        const double dy = my_pos[1] - xj[1];
        const double dz = my_pos[2] - xj[2];
        const double r = std::sqrt(dx*dx + dy*dy + dz*dz);

        if (r <= 1.0e-12 || r >= rc) continue;

        const double w = KernelWeight(r, h, rc);
        if (w <= 0.0) continue;

        const auto& vj = p_neighbour->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        uf[0] += w * vj[0];
        uf[1] += w * vj[1];
        uf[2] += w * vj[2];
        wsum += w;
    }

    if (wsum <= 0.0) {
        return;
    }

    uf[0] /= wsum;
    uf[1] /= wsum;
    uf[2] /= wsum;

    const double a = GetRadius();
    const double coeff = 6.0 * Globals::Pi * mu * a;

    additionally_applied_force[0] += coeff * (uf[0] - my_vel[0]);
    additionally_applied_force[1] += coeff * (uf[1] - my_vel[1]);
    additionally_applied_force[2] += coeff * (uf[2] - my_vel[2]);

    KRATOS_CATCH("")
}

} // namespace Kratos