#include "custom_elements/dpd_particle.h"
#include "DEM_application_variables.h"

namespace Kratos {

// ---------------------------------------------------------------------------
// ComputeBallToBallContactForceAndMoment
// Override to skip the base-class indentation>0 guard so that DPD range
// forces act on non-overlapping pairs too.
// ---------------------------------------------------------------------------
void DPDParticle::ComputeBallToBallContactForceAndMoment(
    ParticleDataBuffer& data_buffer,
    const ProcessInfo& r_process_info,
    array_1d<double, 3>& rElasticForce,
    array_1d<double, 3>& rContactForce)
{
    const int size = (int)mNeighbourElements.size();

    for (int i = 0; i < size; i++) {
        if (mNeighbourElements[i] == nullptr) continue;
        if (!data_buffer.SetNextNeighbourOrExit(i)) break;
        if (!CalculateRelativePositionsOrSkipContact(data_buffer)) continue;

        // data_buffer.mIndentation may be <= 0 for range-force pairs — that is OK

        double LocalElasticContactForce[3]      = {0.0, 0.0, 0.0};
        double DeltDisp[3]                      = {0.0, 0.0, 0.0};
        double LocalDeltDisp[3]                 = {0.0, 0.0, 0.0};
        double RelVel[3]                        = {0.0, 0.0, 0.0};
        double ViscoDampingLocalContactForce[3] = {0.0, 0.0, 0.0};
        double cohesive_force                   = 0.0;
        bool   sliding                          = false;

        array_1d<double, 3> neighbour_elastic_contact_force = ZeroVector(3);

        EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel,
                                  data_buffer.mLocalCoordSystem,
                                  data_buffer.mOldLocalCoordSystem,
                                  GetGeometry()[0].FastGetSolutionStepValue(VELOCITY),
                                  GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT));

        EvaluateBallToBallForcesForPositiveIndentiations(
            data_buffer, r_process_info,
            LocalElasticContactForce,
            DeltDisp, LocalDeltDisp, RelVel,
            data_buffer.mIndentation,
            ViscoDampingLocalContactForce, cohesive_force,
            data_buffer.mpOtherParticle, sliding,
            data_buffer.mLocalCoordSystem,
            data_buffer.mOldLocalCoordSystem,
            neighbour_elastic_contact_force);
    }
}

// ---------------------------------------------------------------------------
// EvaluateBallToBallForcesForPositiveIndentiations
// Calls the constitutive law directly — bypasses base indentation>0 guard.
// ---------------------------------------------------------------------------
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
    // Retrieve the per-neighbour index so we can read mNeighbourElasticContactForces[i]
    // The neighbour pointer lets us locate the index.
    int i_neighbour = -1;
    for (int k = 0; k < (int)mNeighbourElements.size(); k++) {
        if (mNeighbourElements[k] == element2) { i_neighbour = k; break; }
    }

    // Old local elastic contact force for this neighbour (history variable)
    double OldLocalElasticContactForce[3] = {0.0, 0.0, 0.0};
    if (i_neighbour >= 0 && i_neighbour < (int)mNeighbourElasticContactForces.size()) {
        // Rotate stored global force to current local frame
        const array_1d<double,3>& stored = mNeighbourElasticContactForces[i_neighbour];
        for (int k = 0; k < 3; k++) {
            OldLocalElasticContactForce[k] = LocalCoordSystem[k][0]*stored[0]
                                           + LocalCoordSystem[k][1]*stored[1]
                                           + LocalCoordSystem[k][2]*stored[2];
        }
    }

    // Previous indentation: not stored in ParticleDataBuffer; pass 0.0.
    // The DPD-SPH CL does not use it (Hertz spring = 0 for kn=0).
    const double previous_indentation = 0.0;

    DEMDiscontinuumConstitutiveLaw& cl =
        *GetProperties()[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER];

    cl.InitializeContact(this, element2, indentation);

    cl.CalculateForces(
        r_process_info,
        OldLocalElasticContactForce,
        LocalElasticContactForce,
        LocalDeltDisp,
        RelVel,
        indentation,
        previous_indentation,
        ViscoDampingLocalContactForce,
        cohesive_force,
        this,
        element2,
        sliding,
        LocalCoordSystem);

    // Rotate local forces to global and accumulate into node variables
    array_1d<double, 3>& elastic_force_node =
        GetGeometry()[0].FastGetSolutionStepValue(ELASTIC_FORCES);
    array_1d<double, 3>& contact_force_node =
        GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCES);

    for (int j = 0; j < 3; j++) {
        double elastic_global_j = 0.0;
        double damping_global_j = 0.0;
        for (int k = 0; k < 3; k++) {
            elastic_global_j += LocalCoordSystem[k][j] * LocalElasticContactForce[k];
            damping_global_j += LocalCoordSystem[k][j] * ViscoDampingLocalContactForce[k];
        }
        elastic_force_node[j]            += elastic_global_j;
        contact_force_node[j]            += elastic_global_j + damping_global_j;
        neighbour_elastic_contact_force[j] = -elastic_global_j;
    }

    // Update stored history for this neighbour
    if (i_neighbour >= 0 && i_neighbour < (int)mNeighbourElasticContactForces.size()) {
        mNeighbourElasticContactForces[i_neighbour] = neighbour_elastic_contact_force * (-1.0);
    }
}

} // namespace Kratos
