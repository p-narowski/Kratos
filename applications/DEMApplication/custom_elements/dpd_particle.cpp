#include "custom_elements/dpd_particle.h"
#include "custom_elements/suspended_particle.h"

#include "DEM_application_variables.h"
#include "custom_utilities/GeometryFunctions.h"
#include "utilities/math_utils.h"

#include <algorithm>
#include <cmath>

namespace Kratos
{

    namespace
    {

        void BuildDPDLocalBasis(
            const array_1d<double, 3> &rNormalInput,
            double rLocalCoordSystem[3][3])
        {
            array_1d<double, 3> normal;
            normal.clear();

            normal[0] = rNormalInput[0];
            normal[1] = rNormalInput[1];
            normal[2] = rNormalInput[2];

            const double normal_norm = norm_2(normal);

            KRATOS_ERROR_IF(normal_norm <= 1.0e-15)
                << "Cannot build DPD local basis from zero normal."
                << std::endl;

            normal[0] /= normal_norm;
            normal[1] /= normal_norm;
            normal[2] /= normal_norm;

            array_1d<double, 3> reference;
            reference.clear();

            if (std::abs(normal[0]) < 0.9)
            {
                reference[0] = 1.0;
            }
            else
            {
                reference[1] = 1.0;
            }

            // Local tangent 0 = normalized(reference x normal).
            array_1d<double, 3> tangent_0;
            tangent_0.clear();

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
                << "Cannot build first DPD tangent direction."
                << std::endl;

            tangent_0[0] /= tangent_0_norm;
            tangent_0[1] /= tangent_0_norm;
            tangent_0[2] /= tangent_0_norm;

            // Local tangent 1 = normal x tangent_0.
            array_1d<double, 3> tangent_1;
            tangent_1.clear();

            tangent_1[0] =
                normal[1] * tangent_0[2] -
                normal[2] * tangent_0[1];

            tangent_1[1] =
                normal[2] * tangent_0[0] -
                normal[0] * tangent_0[2];

            tangent_1[2] =
                normal[0] * tangent_0[1] -
                normal[1] * tangent_0[0];

            rLocalCoordSystem[0][0] = tangent_0[0];
            rLocalCoordSystem[0][1] = tangent_0[1];
            rLocalCoordSystem[0][2] = tangent_0[2];

            rLocalCoordSystem[1][0] = tangent_1[0];
            rLocalCoordSystem[1][1] = tangent_1[1];
            rLocalCoordSystem[1][2] = tangent_1[2];

            rLocalCoordSystem[2][0] = normal[0];
            rLocalCoordSystem[2][1] = normal[1];
            rLocalCoordSystem[2][2] = normal[2];
        }

        void ProjectGlobalToLocal(
            const array_1d<double, 3> &rGlobalVector,
            const double rLocalCoordSystem[3][3],
            double rLocalVector[3])
        {
            for (int local_direction = 0; local_direction < 3; ++local_direction)
            {
                rLocalVector[local_direction] = 0.0;

                for (int global_direction = 0;
                     global_direction < 3;
                     ++global_direction)
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
            rGlobalVector.clear();

            for (int local_direction = 0; local_direction < 3; ++local_direction)
            {
                for (int global_direction = 0;
                     global_direction < 3;
                     ++global_direction)
                {
                    rGlobalVector[global_direction] +=
                        rLocalCoordSystem[local_direction][global_direction] *
                        rLocalVector[local_direction];
                }
            }
        }

        bool HasDPDProperties(
            const Properties &rContactProperties)
        {
            return rContactProperties.Has(DPD_SMOOTHING_LENGTH) &&
                   rContactProperties.Has(DPD_CUTOFF_RADIUS) &&
                   rContactProperties.Has(DPD_CONSERVATIVE_COEFF) &&
                   rContactProperties.Has(DPD_DISSIPATIVE_COEFF_NORMAL) &&
                   rContactProperties.Has(DPD_DISSIPATIVE_COEFF_TANGENTIAL);
        }

    } // unnamed namespace

    void DPDParticle::ComputeBallToBallContactForceAndMoment(
        ParticleDataBuffer &rDataBuffer,
        const ProcessInfo &rProcessInfo,
        array_1d<double, 3> &rElasticForce,
        array_1d<double, 3> &rContactForce)
    {
        const array_1d<double, 3> &rMyPosition =
            GetGeometry()[0].Coordinates();

        const array_1d<double, 3> &rMyVelocity =
            GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

        for (std::size_t i = 0; i < mNeighbourElements.size(); ++i)
        {
            SphericParticle *p_neighbour = mNeighbourElements[i];

            if (p_neighbour == nullptr)
            {
                continue;
            }

            Properties &rContactProperties =
                GetProperties().GetSubProperties(
                    p_neighbour->GetProperties().Id());

            if (!HasDPDProperties(rContactProperties))
            {
                continue;
            }

            const double rc =
                rContactProperties[DPD_CUTOFF_RADIUS];

            mDiscontinuumConstitutiveLaw =
                pCloneDiscontinuumConstitutiveLawWithNeighbour(
                    p_neighbour);

            KRATOS_ERROR_IF(mDiscontinuumConstitutiveLaw == nullptr)
                << "Null DPD constitutive law for particle "
                << Id()
                << " and neighbour "
                << p_neighbour->Id()
                << std::endl;

            if (mDiscontinuumConstitutiveLaw->GetTypeOfLaw() !=
                "DEM_DPD_SPH_LIKE")
            {
                continue;
            }

            array_1d<double, 3> r_ij;
            r_ij.clear();

            r_ij[0] =
                rMyPosition[0] -
                p_neighbour->GetGeometry()[0].Coordinates()[0];

            r_ij[1] =
                rMyPosition[1] -
                p_neighbour->GetGeometry()[0].Coordinates()[1];

            r_ij[2] =
                rMyPosition[2] -
                p_neighbour->GetGeometry()[0].Coordinates()[2];

            const double distance =
                norm_2(r_ij);

            if (distance <= 1.0e-12 || distance >= rc)
            {
                continue;
            }

            array_1d<double, 3> normal;
            normal.clear();

            normal[0] = r_ij[0] / distance;
            normal[1] = r_ij[1] / distance;
            normal[2] = r_ij[2] / distance;

            double local_coord_system[3][3] = {};
            BuildDPDLocalBasis(
                normal,
                local_coord_system);

            const array_1d<double, 3> &rNeighbourVelocity =
                p_neighbour->GetGeometry()[0]
                    .FastGetSolutionStepValue(VELOCITY);

            array_1d<double, 3> global_relative_velocity;
            global_relative_velocity.clear();

            global_relative_velocity[0] =
                rMyVelocity[0] - rNeighbourVelocity[0];

            global_relative_velocity[1] =
                rMyVelocity[1] - rNeighbourVelocity[1];

            global_relative_velocity[2] =
                rMyVelocity[2] - rNeighbourVelocity[2];

            double local_relative_velocity[3] =
                {0.0, 0.0, 0.0};

            ProjectGlobalToLocal(
                global_relative_velocity,
                local_coord_system,
                local_relative_velocity);

            const double dpd_indentation =
                rc - distance;

            double old_local_elastic_contact_force[3] =
                {0.0, 0.0, 0.0};

            double local_elastic_contact_force[3] =
                {0.0, 0.0, 0.0};

            double local_viscous_contact_force[3] =
                {0.0, 0.0, 0.0};

            double local_delta_displacement[3] =
                {0.0, 0.0, 0.0};

            double cohesive_force = 0.0;
            bool sliding = false;

            mDiscontinuumConstitutiveLaw->InitializeContact(
                this,
                p_neighbour,
                dpd_indentation);

            mDiscontinuumConstitutiveLaw->CalculateForces(
                rProcessInfo,
                old_local_elastic_contact_force,
                local_elastic_contact_force,
                local_delta_displacement,
                local_relative_velocity,
                dpd_indentation,
                dpd_indentation,
                local_viscous_contact_force,
                cohesive_force,
                this,
                p_neighbour,
                sliding,
                local_coord_system);

            array_1d<double, 3> global_elastic_force;
            global_elastic_force.clear();

            array_1d<double, 3> global_viscous_force;
            global_viscous_force.clear();

            ProjectLocalToGlobal(
                local_elastic_contact_force,
                local_coord_system,
                global_elastic_force);

            ProjectLocalToGlobal(
                local_viscous_contact_force,
                local_coord_system,
                global_viscous_force);

            array_1d<double, 3> global_pair_force(
                global_elastic_force);

            noalias(global_pair_force) +=
                global_viscous_force;

            noalias(rElasticForce) +=
                global_elastic_force;

            noalias(rContactForce) +=
                global_pair_force;
        }
    }

    void DPDParticle::ComputeBallToRigidFaceContactForceAndMoment(
        ParticleDataBuffer &rDataBuffer,
        array_1d<double, 3> &rElasticForce,
        array_1d<double, 3> &rContactForce,
        array_1d<double, 3> &rRigidElementForce,
        const ProcessInfo &rProcessInfo)
    {
        static int dpd_wall_override_counter = 0;

        const array_1d<double, 3> &rParticleVelocity =
            GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

        std::vector<DEMWall *> &rWallNeighbours =
            mNeighbourRigidFaces;

        static int dpd_wall_geometry_debug_counter = 0;
        static int dpd_wall_law_debug_counter = 0;
        static int wall_debug_counter = 0;

        for (std::size_t i = 0; i < rWallNeighbours.size(); ++i)
        {
            DEMWall *p_wall = rWallNeighbours[i];

            if (p_wall == nullptr)
            {
                continue;
            }

            Properties &rContactProperties =
                GetProperties().GetSubProperties(
                    p_wall->GetProperties().Id());

            if (!HasDPDProperties(rContactProperties))
            {
                continue;
            }

            const double rc =
                rContactProperties[DPD_CUTOFF_RADIUS];

            const double dpd_wall_range =
                0.5 * rc;

            const double wall_candidate_radius =
                std::max(
                    GetInteractionRadius(),
                    dpd_wall_range);

            double local_coord_system[3][3] = {};
            double distance_to_wall = 0.0;
            int contact_type = -1;

            array_1d<double, 4> &rWeights =
                mContactConditionWeights[i];

            array_1d<double, 3> wall_delta_displacement;
            wall_delta_displacement.clear();

            array_1d<double, 3> wall_velocity;
            wall_velocity.clear();

            p_wall->ComputeConditionRelativeData(
                static_cast<int>(i),
                this,
                local_coord_system,
                distance_to_wall,
                rWeights,
                wall_delta_displacement,
                wall_velocity,
                contact_type,
                wall_candidate_radius);

// #pragma omp critical(DPDWallGeometryDebug)
//             {
//                 if (dpd_wall_geometry_debug_counter < 40)
//                 {
//                     KRATOS_INFO("DPD-WALL-GEOMETRY")
//                         << "particle=" << Id()
//                         << " wall=" << p_wall->Id()
//                         << " contact_type=" << contact_type
//                         << " d=" << distance_to_wall
//                         << " candidate_radius="
//                         << wall_candidate_radius
//                         << " rc=" << rc
//                         << " rc_over_2="
//                         << dpd_wall_range
//                         << std::endl;

//                     ++dpd_wall_geometry_debug_counter;
//                 }
//             }

            if (contact_type != 1 &&
                contact_type != 2 &&
                contact_type != 3)
            {
                continue;
            }

            // The ghost-wall DPD force is active only when r_ghost = 2d < rc.
            if (distance_to_wall <= 1.0e-12 ||
                distance_to_wall >= dpd_wall_range)
            {
                continue;
            }

            mDiscontinuumConstitutiveLaw =
                pCloneDiscontinuumConstitutiveLawWithFEMNeighbour(
                    p_wall);

            KRATOS_ERROR_IF(mDiscontinuumConstitutiveLaw == nullptr)
                << "Null DPD FEM constitutive law for particle "
                << Id()
                << " and wall "
                << p_wall->Id()
                << std::endl;

// #pragma omp critical(DPDWallLawDebug)
//             {
//                 if (dpd_wall_law_debug_counter < 20)
//                 {
//                     KRATOS_INFO("DPD-WALL-LAW")
//                         << "particle=" << Id()
//                         << " wall=" << p_wall->Id()
//                         << " particle_property="
//                         << GetProperties().Id()
//                         << " wall_property="
//                         << p_wall->GetProperties().Id()
//                         << " law_type="
//                         << mDiscontinuumConstitutiveLaw
//                                ->GetTypeOfLaw()
//                         << std::endl;

//                     ++dpd_wall_law_debug_counter;
//                 }
//             }

            if (mDiscontinuumConstitutiveLaw->GetTypeOfLaw() !=
                "DEM_DPD_SPH_LIKE")
            {
                continue;
            }

            /*
             * Preserve the physical distance relation used by
             * DEM_DPD_SPH_LIKE::CalculateForcesWithFEM:
             *
             * distance_to_wall =
             * GetInteractionRadius() - indentation.
             */
            const double indentation =
                GetInteractionRadius() - distance_to_wall;

            array_1d<double, 3> global_relative_velocity;
            global_relative_velocity.clear();

            global_relative_velocity[0] =
                rParticleVelocity[0] - wall_velocity[0];

            global_relative_velocity[1] =
                rParticleVelocity[1] - wall_velocity[1];

            global_relative_velocity[2] =
                rParticleVelocity[2] - wall_velocity[2];

            double local_relative_velocity[3] =
                {0.0, 0.0, 0.0};

            ProjectGlobalToLocal(
                global_relative_velocity,
                local_coord_system,
                local_relative_velocity);

            double old_local_elastic_contact_force[3] =
                {0.0, 0.0, 0.0};

            double local_elastic_contact_force[3] =
                {0.0, 0.0, 0.0};

            double local_viscous_contact_force[3] =
                {0.0, 0.0, 0.0};

            double local_delta_displacement[3] =
                {0.0, 0.0, 0.0};

            double cohesive_force = 0.0;
            bool sliding = false;

            mDiscontinuumConstitutiveLaw->InitializeContactWithFEM(
                this,
                p_wall,
                indentation,
                0.0);

            /*
             * Temporary deterministic diagnostic.
             *
             * The log already showed this particle-wall pair has:
             * contact_type = 1
             * d = 0.047870
             * 2d = 0.095740 < rc = 0.18.
             *
             * Therefore this message must appear if this source and this code
             * path are actually being executed.
             */
            // if (Id() == 501 && p_wall->Id() == 3)
            // {
            //     KRATOS_INFO("DPD-FEM-CALLSITE-TEST")
            //         << "particle=" << Id()
            //         << " wall=" << p_wall->Id()
            //         << " contact_type=" << contact_type
            //         << " d=" << distance_to_wall
            //         << " indentation=" << indentation
            //         << " r_ghost=" << 2.0 * distance_to_wall
            //         << " rc=" << rc
            //         << " local_v=("
            //         << local_relative_velocity[0] << ", "
            //         << local_relative_velocity[1] << ", "
            //         << local_relative_velocity[2] << ")"
            //         << std::endl;
            // }

            mDiscontinuumConstitutiveLaw->CalculateForcesWithFEM(
                rProcessInfo,
                old_local_elastic_contact_force,
                local_elastic_contact_force,
                local_delta_displacement,
                local_relative_velocity,
                indentation,
                indentation,
                local_viscous_contact_force,
                cohesive_force,
                this,
                p_wall,
                sliding);

            array_1d<double, 3> global_elastic_force;
            global_elastic_force.clear();

            array_1d<double, 3> global_viscous_force;
            global_viscous_force.clear();

            ProjectLocalToGlobal(
                local_elastic_contact_force,
                local_coord_system,
                global_elastic_force);

            ProjectLocalToGlobal(
                local_viscous_contact_force,
                local_coord_system,
                global_viscous_force);

            array_1d<double, 3> global_total_force(
                global_elastic_force);

            noalias(global_total_force) +=
                global_viscous_force;

            noalias(rElasticForce) +=
                global_elastic_force;

            noalias(rContactForce) +=
                global_total_force;

            /*
             * The particle-side rigid-element-force accumulator stores the
             * equal-and-opposite reaction sent toward the FEM wall.
             */
            noalias(rRigidElementForce) -=
                global_total_force;

// #pragma omp critical(DPDWallForceDebug)
//             {
//                 if (wall_debug_counter < 10)
//                 {
//                     KRATOS_INFO("DPD-WALL")
//                         << "particle=" << Id()
//                         << " wall=" << p_wall->Id()
//                         << " d=" << distance_to_wall
//                         << " 2d="
//                         << 2.0 * distance_to_wall
//                         << " rc=" << rc
//                         << " F=("
//                         << global_total_force[0] << ", "
//                         << global_total_force[1] << ", "
//                         << global_total_force[2] << ")"
//                         << std::endl;

//                     ++wall_debug_counter;
//                 }
//             }
        }
    }

} // namespace Kratos