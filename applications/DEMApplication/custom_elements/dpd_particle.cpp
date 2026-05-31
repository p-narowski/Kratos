#include "custom_elements/dpd_particle.h"
#include "DEM_application_variables.h"

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
    // Call the DPD constitutive law directly, bypassing the base class
    // indentation > 0 guard that would silently skip all range-force pairs.
    DEMDiscontinuumConstitutiveLaw& cl =
        *GetProperties()[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER];

    cl.InitializeContact(this, element2, indentation);

    cl.CalculateForces(
        r_process_info,
        data_buffer.mOldLocalElasticContactForce,
        LocalElasticContactForce,
        LocalDeltDisp,
        RelVel,
        indentation,
        data_buffer.mPreviousIndentation,
        ViscoDampingLocalContactForce,
        cohesive_force,
        this,
        element2,
        sliding,
        LocalCoordSystem);

    // Rotate local forces to global and accumulate
    double TotalContactForce[3];
    for (int k = 0; k < 3; k++) {
        TotalContactForce[k] = LocalElasticContactForce[k] + ViscoDampingLocalContactForce[k];
    }

    // Local-to-global rotation: F_global = R^T * F_local
    // LocalCoordSystem rows are the local axes expressed in global coords
    array_1d<double, 3>& elastic_force = GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
    array_1d<double, 3>& total_force   = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCES);

    for (int j = 0; j < 3; j++) {
        double elastic_global_j = 0.0;
        double total_global_j   = 0.0;
        for (int k = 0; k < 3; k++) {
            elastic_global_j += LocalCoordSystem[k][j] * LocalElasticContactForce[k];
            total_global_j   += LocalCoordSystem[k][j] * TotalContactForce[k];
        }
        elastic_force[j] += elastic_global_j;
        total_force[j]   += total_global_j;
        neighbour_elastic_contact_force[j] = -elastic_global_j;
    }
}

} // namespace Kratos
