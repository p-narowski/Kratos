#include "custom_elements/dpd_particle.h"

namespace Kratos {

void DPDParticle::ComputeBallToBallContactForceAndMoment(
    ParticleDataBuffer& data_buffer,
    const ProcessInfo& r_process_info,
    array_1d<double, 3>& rElasticForce,
    array_1d<double, 3>& rContactForce)
{
    // DPD range force: skip the indentation > 0 guard in the base class.
    // CalculateRelativePositionsOrSkipContact sets data_buffer.mIndentation
    // and mDistance. We call the base only when the pair is within rc.
    // The base SphericParticle::ComputeBallToBallContactForceAndMoment
    // exits early for indentation < 0, so we bypass it entirely and call
    // EvaluateBallToBallForcesForPositiveIndentiations directly.

    const int size = mNeighbourElements.size();
    for (int i = 0; i < size; i++) {
        if (mNeighbourElements[i] == nullptr) continue;
        if (!data_buffer.SetNextNeighbourOrExit(i)) break;
        if (!CalculateRelativePositionsOrSkipContact(data_buffer)) continue;

        // Allow non-overlapping pairs through (range force)
        // data_buffer.mIndentation may be <= 0 — that is fine for DPD

        double LocalElasticContactForce[3]      = {0.0, 0.0, 0.0};
        double LocalDeltDisp[3]                 = {0.0, 0.0, 0.0};
        double LocalRelVel[3]                   = {0.0, 0.0, 0.0};
        double DeltDisp[3]                      = {0.0, 0.0, 0.0};
        double ViscoDampingLocalContactForce[3] = {0.0, 0.0, 0.0};
        double cohesive_force                   = 0.0;
        bool   sliding                          = false;

        array_1d<double, 3> neighbour_elastic_contact_force = ZeroVector(3);

        EvaluateDeltaDisplacement(data_buffer, DeltDisp, LocalRelVel,
                                  data_buffer.mLocalCoordSystem,
                                  data_buffer.mOldLocalCoordSystem,
                                  GetGeometry()[0].FastGetSolutionStepValue(VELOCITY),
                                  GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT));

        EvaluateBallToBallForcesForPositiveIndentiations(
            data_buffer, r_process_info,
            LocalElasticContactForce, DeltDisp, LocalDeltDisp, LocalRelVel,
            data_buffer.mIndentation,
            ViscoDampingLocalContactForce, cohesive_force,
            data_buffer.mpOtherParticle, sliding,
            data_buffer.mLocalCoordSystem, data_buffer.mOldLocalCoordSystem,
            neighbour_elastic_contact_force);
    }
}
void DPDParticle::EvaluateBallToBallForcesForPositiveIndentiations(
    SphericParticle::ParticleDataBuffer& data_buffer,
    const ProcessInfo& r_process_info,
    double LocalElasticContactForce[3],
    double DeltDisp[3],
    double LocalDeltDisp[3],
    double RelVel[3],
    double indentation,
    double ViscoDampingLocalContactForce[3],
    double& cohesive_force,
    SphericParticle* element2,
    bool& sliding,
    double LocalCoordSystem[3][3],
    double OldLocalCoordSystem[3][3],
    array_1d<double, 3>& neighbour_elastic_contact_force)
{
    // Delegate to base — allows range force (indentation may be <= 0)
    SphericParticle::EvaluateBallToBallForcesForPositiveIndentiations(
        data_buffer, r_process_info,
        LocalElasticContactForce, DeltDisp, LocalDeltDisp, RelVel,
        indentation, ViscoDampingLocalContactForce, cohesive_force,
        element2, sliding, LocalCoordSystem, OldLocalCoordSystem,
        neighbour_elastic_contact_force);
}

} // namespace Kratos