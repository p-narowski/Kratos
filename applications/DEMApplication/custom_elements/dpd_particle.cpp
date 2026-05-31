#include "custom_elements/dpd_particle.h"
#include "DEM_application_variables.h"

namespace Kratos {

// ---------------------------------------------------------------------------
// EvaluateBallToBallForcesForPositiveIndentiations
//
// Called by the base-class ComputeBallToBallContactForceAndMoment loop for
// every neighbour (after CalculateRelativePositionsOrSkipContact passes).
// The base class normally guards this call with  indentation > 0 ; we bypass
// that guard so that DPD/SPH range forces act on non-overlapping pairs too.
//
// We call mDiscontinuumConstitutiveLaw — the per-particle CLONED instance
// that has already been initialised by Kratos — NOT the raw prototype stored
// in Properties.  The base class then rotates LocalElasticContactForce and
// ViscoDampingLocalContactForce to global frame and accumulates them into
// rElasticForce / rContactForce / TOTAL_FORCES correctly.
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
    // Locate history index for this neighbour
    double OldLocalElasticContactForce[3] = {0.0, 0.0, 0.0};
    int i_neighbour = -1;
    for (int k = 0; k < (int)mNeighbourElements.size(); k++) {
        if (mNeighbourElements[k] == element2) { i_neighbour = k; break; }
    }
    if (i_neighbour >= 0 && i_neighbour < (int)mNeighbourElasticContactForces.size()) {
        const array_1d<double,3>& stored = mNeighbourElasticContactForces[i_neighbour];
        for (int k = 0; k < 3; k++) {
            OldLocalElasticContactForce[k] = LocalCoordSystem[k][0]*stored[0]
                                           + LocalCoordSystem[k][1]*stored[1]
                                           + LocalCoordSystem[k][2]*stored[2];
        }
    }

    // Use the per-particle cloned CL instance (NOT Properties prototype)
    mDiscontinuumConstitutiveLaw->InitializeContact(this, element2, indentation);

    mDiscontinuumConstitutiveLaw->CalculateForces(
        r_process_info,
        OldLocalElasticContactForce,
        LocalElasticContactForce,
        LocalDeltDisp,
        RelVel,
        indentation,
        0.0, // previous_indentation — DPD-SPH CL does not need it
        ViscoDampingLocalContactForce,
        cohesive_force,
        this,
        element2,
        sliding,
        LocalCoordSystem);
}

} // namespace Kratos
