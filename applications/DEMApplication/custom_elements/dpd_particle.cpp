#include "custom_elements/dpd_particle.h"
#include "DEM_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos {

// ---------------------------------------------------------------------------
// ComputeBallToBallContactForceAndMoment
//
// Full override of the base-class method.  We cannot rely on the base-class
// loop because it guards the constitutive-law call with (indentation > 0),
// which is never true for DPD/SPH soft-core pairs that do not geometrically
// overlap.
//
// Strategy:
//   1. Call the base class first so that any truly overlapping contacts
//      (e.g. wall-particle pairs using the Hertz law on Property 2) are
//      handled by the existing machinery without duplication.
//   2. Loop all search-radius neighbours ourselves and call the DPD CL
//      directly for every pair with dist < rc, where rc is read from the
//      contact sub-properties (material_relations in MaterialsDEM.json).
//      This makes the DPD interaction range purely a material-law parameter,
//      fully independent of the geometric particle radii.
// ---------------------------------------------------------------------------
void DPDParticle::ComputeBallToBallContactForceAndMoment(
    ParticleDataBuffer& data_buffer,
    const ProcessInfo& r_process_info,
    array_1d<double, 3>& rElasticForce,
    array_1d<double, 3>& rContactForce)
{
    // --- Step 1: base-class handles overlapping contacts (Hertz, walls, etc.) ---
    SphericParticle::ComputeBallToBallContactForceAndMoment(
        data_buffer, r_process_info, rElasticForce, rContactForce);

    // --- Step 2: DPD range forces for all neighbours within rc ---
    const array_1d<double, 3>& my_pos = GetGeometry()[0].Coordinates();
    const array_1d<double, 3>& my_vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

    for (unsigned int i = 0; i < mNeighbourElements.size(); ++i) {
        SphericParticle* neighbour = mNeighbourElements[i];
        if (!neighbour) continue;

        // DPD cut-off radius from the contact sub-properties (material_relations).
        // DPD_CUTOFF_RADIUS lives there, NOT on the element's own base properties.
        Properties& contact_props =
            GetProperties().GetSubProperties(neighbour->GetProperties().Id());
        if (!contact_props.Has(DPD_CUTOFF_RADIUS)) continue;
        const double rc = contact_props[DPD_CUTOFF_RADIUS];

        // Vector from neighbour to me (same sign convention as base class)
        array_1d<double, 3> r_ij;
        noalias(r_ij) = my_pos - neighbour->GetGeometry()[0].Coordinates();
        const double dist = norm_2(r_ij);
        if (dist < 1.0e-12) continue;

        // Pure DPD range condition — independent of geometric radii
        if (dist >= rc) continue;

        // "indentation" in the DPD sense: rc - dist > 0
        const double indentation = rc - dist;

        // Build a minimal local coordinate system (normal = r_ij / dist)
        double LocalCoordSystem[3][3] = {{0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0},
                                         {0.0, 0.0, 0.0}};
        for (int k = 0; k < 3; ++k)
            LocalCoordSystem[2][k] = r_ij[k] / dist; // normal (local z)

        // Relative velocity
        const array_1d<double, 3>& nb_vel =
            neighbour->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> rel_vel;
        noalias(rel_vel) = my_vel - nb_vel;

        double RelVel[3]          = {rel_vel[0], rel_vel[1], rel_vel[2]};
        double LocalDeltDisp[3]   = {0.0, 0.0, 0.0};

        // Recover previous elastic force in local frame from history
        double OldLocalElasticContactForce[3] = {0.0, 0.0, 0.0};
        if (i < mNeighbourElasticContactForces.size()) {
            const array_1d<double, 3>& stored = mNeighbourElasticContactForces[i];
            for (int k = 0; k < 3; ++k)
                OldLocalElasticContactForce[k] =
                    LocalCoordSystem[k][0] * stored[0] +
                    LocalCoordSystem[k][1] * stored[1] +
                    LocalCoordSystem[k][2] * stored[2];
        }

        double LocalElasticContactForce[3]      = {0.0, 0.0, 0.0};
        double ViscoDampingLocalContactForce[3]  = {0.0, 0.0, 0.0};
        double cohesive_force = 0.0;
        bool   sliding        = false;

        // Invoke the per-particle cloned DPD constitutive law
        mDiscontinuumConstitutiveLaw->InitializeContact(this, neighbour, indentation);

        mDiscontinuumConstitutiveLaw->CalculateForces(
            r_process_info,
            OldLocalElasticContactForce,
            LocalElasticContactForce,
            LocalDeltDisp,
            RelVel,
            indentation,
            0.0, // previous_indentation — not needed by DPD-SPH CL
            ViscoDampingLocalContactForce,
            cohesive_force,
            this,
            neighbour,
            sliding,
            LocalCoordSystem);

        // Rotate local forces back to global frame and accumulate
        for (int k = 0; k < 3; ++k) {
            const double f_elastic = LocalElasticContactForce[k];
            const double f_visco   = ViscoDampingLocalContactForce[k];

            double global_elastic[3] = {0.0, 0.0, 0.0};
            double global_visco[3]   = {0.0, 0.0, 0.0};
            for (int d = 0; d < 3; ++d) {
                global_elastic[d] += LocalCoordSystem[k][d] * f_elastic;
                global_visco[d]   += LocalCoordSystem[k][d] * f_visco;
            }
            for (int d = 0; d < 3; ++d) {
                rElasticForce[d] += global_elastic[d];
                rContactForce[d] += global_elastic[d] + global_visco[d];
            }
        }

        // Persist elastic force in history for next step
        if (i < mNeighbourElasticContactForces.size()) {
            for (int d = 0; d < 3; ++d) {
                mNeighbourElasticContactForces[i][d] = 0.0;
                for (int k = 0; k < 3; ++k)
                    mNeighbourElasticContactForces[i][d] +=
                        LocalCoordSystem[k][d] * LocalElasticContactForce[k];
            }
        }
    }
}

} // namespace Kratos
