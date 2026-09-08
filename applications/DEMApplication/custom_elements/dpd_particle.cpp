#include "custom_elements/dpd_particle.h"
#include "DEM_application_variables.h"
#include "utilities/math_utils.h"
#include <cmath>

namespace Kratos
{
    namespace
    {

        void BuildDPDLocalBasis(
            const array_1d<double, 3> &rNormal,
            double rLocalCoordSystem[3][3])
        {
            // Local direction 2: unit pair normal.
            array_1d<double, 3> normal = rNormal;

            const double normal_norm = norm_2(normal);

            KRATOS_ERROR_IF(normal_norm <= 1.0e-15)
                << "Cannot construct DPD basis from a zero normal.";

            for (int d = 0; d < 3; ++d)
            {
                normal[d] /= normal_norm;
            }

            // Select a reference direction safely non-parallel to normal.
            array_1d<double, 3> reference{0.0, 0.0, 0.0};

            if (std::abs(normal[0]) < 0.9)
            {
                reference[0] = 1.0;
            }
            else
            {
                reference[1] = 1.0;
            }

            // Local direction 0: t0 = normalized(reference x normal).
            array_1d<double, 3> tangent_0{0.0, 0.0, 0.0};

            tangent_0[0] =
                reference[1] * normal[2] -
                reference[2] * normal[1];

            tangent_0[1] =
                reference[2] * normal[0] -
                reference[0] * normal[2];

            tangent_0[2] =
                reference[0] * normal[1] -
                reference[1] * normal[0];

            const double tangent_0_norm = norm_2(tangent_0);

            KRATOS_ERROR_IF(tangent_0_norm <= 1.0e-15)
                << "Cannot construct DPD tangent direction.";

            for (int d = 0; d < 3; ++d)
            {
                tangent_0[d] /= tangent_0_norm;
            }

            // Local direction 1: t1 = normal x t0.
            array_1d<double, 3> tangent_1{0.0, 0.0, 0.0};

            tangent_1[0] =
                normal[1] * tangent_0[2] -
                normal[2] * tangent_0[1];

            tangent_1[1] =
                normal[2] * tangent_0[0] -
                normal[0] * tangent_0[2];

            tangent_1[2] =
                normal[0] * tangent_0[1] -
                normal[1] * tangent_0[0];

            for (int d = 0; d < 3; ++d)
            {
                rLocalCoordSystem[0][d] = tangent_0[d];
                rLocalCoordSystem[1][d] = tangent_1[d];
                rLocalCoordSystem[2][d] = normal[d];
            }
        }

        void ProjectGlobalToLocal(
            const array_1d<double, 3> &rGlobalVector,
            const double rLocalCoordSystem[3][3],
            double rLocalVector[3])
        {
            for (int local_direction = 0; local_direction < 3; ++local_direction)
            {
                rLocalVector[local_direction] = 0.0;

                for (int global_direction = 0; global_direction < 3; ++global_direction)
                {
                    rLocalVector[local_direction] +=
                        rLocalCoordSystem[local_direction][global_direction] *
                        rGlobalVector[global_direction];
                }
            }
        }

        void ProjectLocalToGlobal(
            const double rLocalVector[3],
            const double rLocalCoordSystem[3][3],
            array_1d<double, 3> &rGlobalVector)
        {
            rGlobalVector[0] = 0.0;
            rGlobalVector[1] = 0.0;
            rGlobalVector[2] = 0.0;

            for (int local_direction = 0; local_direction < 3; ++local_direction)
            {
                for (int global_direction = 0; global_direction < 3; ++global_direction)
                {
                    rGlobalVector[global_direction] +=
                        rLocalCoordSystem[local_direction][global_direction] *
                        rLocalVector[local_direction];
                }
            }
        }

    } // unnamed namespace
    void DPDParticle::ComputeBallToBallContactForceAndMoment(
        ParticleDataBuffer &data_buffer,
        const ProcessInfo &r_process_info,
        array_1d<double, 3> &rElasticForce,
        array_1d<double, 3> &rContactForce)
    {
        // Do not call the base ball-to-ball DEM contact routine for DPD particles.
        // --- Step 2: DPD range forces for all neighbours within rc ---
        const array_1d<double, 3> &my_pos = GetGeometry()[0].Coordinates();
        const array_1d<double, 3> &my_vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

        for (unsigned int i = 0; i < mNeighbourElements.size(); ++i)
        {
            SphericParticle *neighbour = mNeighbourElements[i];
            if (!neighbour)
                continue;

            Properties &contact_props =
                GetProperties().GetSubProperties(neighbour->GetProperties().Id());

            if (!contact_props.Has(DPD_SMOOTHING_LENGTH) ||
                !contact_props.Has(DPD_CUTOFF_RADIUS) ||
                !contact_props.Has(DPD_CONSERVATIVE_COEFF) ||
                !contact_props.Has(DPD_DISSIPATIVE_COEFF_NORMAL) ||
                !contact_props.Has(DPD_DISSIPATIVE_COEFF_TANGENTIAL))
            {
                continue;
            }

            const double rc = contact_props[DPD_CUTOFF_RADIUS];

            mDiscontinuumConstitutiveLaw =
                pCloneDiscontinuumConstitutiveLawWithNeighbour(neighbour);

            if (!mDiscontinuumConstitutiveLaw)
            {
                continue;
            }

            if (mDiscontinuumConstitutiveLaw->GetTypeOfLaw() != "DPD_SPH_LIKE")
            {
                continue;
            }

            // Vector from neighbour to me (same sign convention as base class)
            array_1d<double, 3> r_ij;
            noalias(r_ij) =
                my_pos - neighbour->GetGeometry()[0].Coordinates();

            const double dist = norm_2(r_ij);

            if (dist <= 1.0e-12)
            {
                continue;
            }

            if (dist >= rc)
            {
                continue;
            }

            const double indentation = rc - dist;

            array_1d<double, 3> normal{0.0, 0.0, 0.0};

            for (int d = 0; d < 3; ++d)
            {
                normal[d] = r_ij[d] / dist;
            }

            double LocalCoordSystem[3][3] = {};
            BuildDPDLocalBasis(normal, LocalCoordSystem);

            // Relative velocity
            const array_1d<double, 3> &nb_vel =
                neighbour->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3> rel_vel;
            noalias(rel_vel) = my_vel - nb_vel;

            double LocalRelVel[3] = {0.0, 0.0, 0.0};
            ProjectGlobalToLocal(rel_vel, LocalCoordSystem, LocalRelVel);
            double LocalDeltDisp[3] = {0.0, 0.0, 0.0};

            // Recover previous elastic force in local frame from history
            double OldLocalElasticContactForce[3] = {0.0, 0.0, 0.0};

            double LocalElasticContactForce[3] = {0.0, 0.0, 0.0};
            double ViscoDampingLocalContactForce[3] = {0.0, 0.0, 0.0};
            double cohesive_force = 0.0;
            bool sliding = false;

            // Ensure we have a valid CL instance for this neighbour
            mDiscontinuumConstitutiveLaw =
                pCloneDiscontinuumConstitutiveLawWithNeighbour(neighbour);
            if (!mDiscontinuumConstitutiveLaw)
            {
                continue;
            }
            // Guard against accidentally running the DPD branch with another law
            if (mDiscontinuumConstitutiveLaw->GetTypeOfLaw() != "DPD_SPH_LIKE")
            {
                continue;
            }
            mDiscontinuumConstitutiveLaw->InitializeContact(
                this,
                neighbour,
                indentation);
            mDiscontinuumConstitutiveLaw->CalculateForces(
                r_process_info,
                OldLocalElasticContactForce,
                LocalElasticContactForce,
                LocalDeltDisp,
                LocalRelVel,
                indentation,
                0.0, // previous_indentation — not needed by DPD-SPH CL
                ViscoDampingLocalContactForce,
                cohesive_force,
                this,
                neighbour,
                sliding,
                LocalCoordSystem);
            // Rotate local forces back to global frame and accumulate
            array_1d<double, 3> global_elastic_force{0.0, 0.0, 0.0};

            array_1d<double, 3> global_viscous_force{0.0, 0.0, 0.0};

            ProjectLocalToGlobal(
                LocalElasticContactForce,
                LocalCoordSystem,
                global_elastic_force);

            ProjectLocalToGlobal(
                ViscoDampingLocalContactForce,
                LocalCoordSystem,
                global_viscous_force);

            array_1d<double, 3> global_pair_force = global_elastic_force;
            noalias(global_pair_force) += global_viscous_force;

            noalias(rElasticForce) += global_elastic_force;
            noalias(rContactForce) += global_pair_force;
        }
    }

} // namespace Kratos
