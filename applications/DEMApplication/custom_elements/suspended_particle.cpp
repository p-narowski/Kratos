#include "custom_elements/suspended_particle.h"

#include "custom_elements/dpd_particle.h"
#include "DEM_application_variables.h"
#include "custom_utilities/GeometryFunctions.h"

#include <cmath>

namespace Kratos
{

namespace
{

void BuildDPDLocalBasis(
    const array_1d<double, 3>& rNormalInput,
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

    if (std::abs(normal[0]) < 0.9) {
        reference[0] = 1.0;
    } else {
        reference[1] = 1.0;
    }

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
    const array_1d<double, 3>& rGlobalVector,
    const double rLocalCoordSystem[3][3],
    double rLocalVector[3])
{
    for (int local_direction = 0; local_direction < 3; ++local_direction) {
        rLocalVector[local_direction] = 0.0;

        for (int global_direction = 0;
             global_direction < 3;
             ++global_direction) {
            rLocalVector[local_direction] +=
                rLocalCoordSystem[local_direction][global_direction]
                * rGlobalVector[global_direction];
        }
    }
}

void ProjectLocalToGlobal(
    const double rLocalVector[3],
    const double rLocalCoordSystem[3][3],
    array_1d<double, 3>& rGlobalVector)
{
    rGlobalVector.clear();

    for (int local_direction = 0; local_direction < 3; ++local_direction) {
        for (int global_direction = 0;
             global_direction < 3;
             ++global_direction) {
            rGlobalVector[global_direction] +=
                rLocalCoordSystem[local_direction][global_direction]
                * rLocalVector[local_direction];
        }
    }
}

bool HasDPDProperties(const Properties& rContactProperties)
{
    return rContactProperties.Has(DPD_SMOOTHING_LENGTH)
        && rContactProperties.Has(DPD_CUTOFF_RADIUS)
        && rContactProperties.Has(DPD_CONSERVATIVE_COEFF)
        && rContactProperties.Has(DPD_DISSIPATIVE_COEFF_NORMAL)
        && rContactProperties.Has(DPD_DISSIPATIVE_COEFF_TANGENTIAL);
}

} // unnamed namespace

void SuspendedParticle::ComputeBallToBallContactForceAndMoment(
    ParticleDataBuffer& rDataBuffer,
    const ProcessInfo& rProcessInfo,
    array_1d<double, 3>& rElasticForce,
    array_1d<double, 3>& rContactForce)
{
    /*
     * The normal DEM implementation must still process:
     *
     * - SuspendedParticle--SuspendedParticle,
     * - SuspendedParticle--ordinary SphericParticle,
     * - all other non-DPD spherical DEM elements.
     *
     * It cannot be called directly after processing DPD neighbours,
     * because it would also apply ordinary DEM contact to the DPD neighbours.
     *
     * Therefore this implementation reproduces the DPDParticle range-force
     * loop only for DPD neighbours. The standard DEM loop is preserved by
     * temporarily removing DPD neighbours, calling the base method, then
     * restoring the original neighbour list.
     */

    std::vector<SphericParticle*> dpd_neighbours;
    std::vector<std::size_t> dpd_neighbour_indices;

    dpd_neighbours.reserve(mNeighbourElements.size());
    dpd_neighbour_indices.reserve(mNeighbourElements.size());

    for (std::size_t i = 0; i < mNeighbourElements.size(); ++i) {
        SphericParticle* p_neighbour = mNeighbourElements[i];

        if (p_neighbour == nullptr) {
            continue;
        }

        if (dynamic_cast<DPDParticle*>(p_neighbour) != nullptr) {
            dpd_neighbours.push_back(p_neighbour);
            dpd_neighbour_indices.push_back(i);
            mNeighbourElements[i] = nullptr;
        }
    }

    /*
     * This executes unchanged DEM contact for:
     * Suspended--Suspended and Suspended--SphericParticle.
     */
    SphericParticle::ComputeBallToBallContactForceAndMoment(
        rDataBuffer,
        rProcessInfo,
        rElasticForce,
        rContactForce);

    /*
     * Restore the exact original neighbour vector before returning.
     * Neighbour-list integrity is important for history, search and output.
     */
    for (std::size_t k = 0; k < dpd_neighbours.size(); ++k) {
        mNeighbourElements[dpd_neighbour_indices[k]] = dpd_neighbours[k];
    }

    const array_1d<double, 3>& r_my_position =
        GetGeometry()[0].Coordinates();

    const array_1d<double, 3>& r_my_velocity =
        GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

    /*
     * Apply DPD interaction only to SuspendedParticle--DPDParticle pairs.
     *
     * The force is accumulated only on this suspended particle's
     * rElasticForce/rContactForce. This follows the existing DEM/DPD
     * force assembly convention in your DPDParticle implementation:
     * each element computes its own force contribution during its own RHS.
     */
    for (SphericParticle* p_neighbour : dpd_neighbours) {
        DPDParticle* p_dpd_neighbour =
            dynamic_cast<DPDParticle*>(p_neighbour);

        if (p_dpd_neighbour == nullptr) {
            continue;
        }

        Properties& r_contact_properties =
            GetProperties().GetSubProperties(
                p_dpd_neighbour->GetProperties().Id());

        if (!HasDPDProperties(r_contact_properties)) {
            continue;
        }

        const double rc =
            r_contact_properties[DPD_CUTOFF_RADIUS];

        mDiscontinuumConstitutiveLaw =
            pCloneDiscontinuumConstitutiveLawWithNeighbour(
                p_dpd_neighbour);

        KRATOS_ERROR_IF(mDiscontinuumConstitutiveLaw == nullptr)
            << "Null DPD constitutive law for suspended particle "
            << Id()
            << " and DPD neighbour "
            << p_dpd_neighbour->Id()
            << std::endl;

        if (mDiscontinuumConstitutiveLaw->GetTypeOfLaw()
            != "DEM_DPD_SPH_LIKE") {
            continue;
        }

        array_1d<double, 3> r_ij;
        r_ij.clear();

        const array_1d<double, 3>& r_neighbour_position =
            p_dpd_neighbour->GetGeometry()[0].Coordinates();

        r_ij[0] = r_my_position[0] - r_neighbour_position[0];
        r_ij[1] = r_my_position[1] - r_neighbour_position[1];
        r_ij[2] = r_my_position[2] - r_neighbour_position[2];

        const double distance = norm_2(r_ij);

        if (distance <= 1.0e-12 || distance >= rc) {
            continue;
        }

        array_1d<double, 3> normal;
        normal.clear();

        normal[0] = r_ij[0] / distance;
        normal[1] = r_ij[1] / distance;
        normal[2] = r_ij[2] / distance;

        double local_coord_system[3][3] = {};
        BuildDPDLocalBasis(normal, local_coord_system);

        const array_1d<double, 3>& r_neighbour_velocity =
            p_dpd_neighbour->GetGeometry()[0]
                .FastGetSolutionStepValue(VELOCITY);

        array_1d<double, 3> global_relative_velocity;
        global_relative_velocity.clear();

        global_relative_velocity[0] =
            r_my_velocity[0] - r_neighbour_velocity[0];

        global_relative_velocity[1] =
            r_my_velocity[1] - r_neighbour_velocity[1];

        global_relative_velocity[2] =
            r_my_velocity[2] - r_neighbour_velocity[2];

        double local_relative_velocity[3] = {
            0.0, 0.0, 0.0};

        ProjectGlobalToLocal(
            global_relative_velocity,
            local_coord_system,
            local_relative_velocity);

        /*
         * DPD indentation is not geometric sphere overlap.
         * It is positive over the DPD interaction range.
         */
        const double dpd_indentation = rc - distance;

        double old_local_elastic_contact_force[3] = {
            0.0, 0.0, 0.0};

        double local_elastic_contact_force[3] = {
            0.0, 0.0, 0.0};

        double local_viscous_contact_force[3] = {
            0.0, 0.0, 0.0};

        double local_delta_displacement[3] = {
            0.0, 0.0, 0.0};

        double cohesive_force = 0.0;
        bool sliding = false;

        mDiscontinuumConstitutiveLaw->InitializeContact(
            this,
            p_dpd_neighbour,
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
            p_dpd_neighbour,
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

        array_1d<double, 3> global_pair_force(global_elastic_force);

        noalias(global_pair_force) += global_viscous_force;

        noalias(rElasticForce) += global_elastic_force;

        noalias(rContactForce) += global_pair_force;
    }
}

void SuspendedParticle::ComputeBallToRigidFaceContactForceAndMoment(
    ParticleDataBuffer& rDataBuffer,
    array_1d<double, 3>& rElasticForce,
    array_1d<double, 3>& rContactForce,
    array_1d<double, 3>& rRigidElementForce,
    const ProcessInfo& rProcessInfo)
{
    /*
     * Suspended particles remain ordinary DEM particles against
     * the screw, barrel and all other rigid FEM walls.
     */
    SphericParticle::ComputeBallToRigidFaceContactForceAndMoment(
        rDataBuffer,
        rElasticForce,
        rContactForce,
        rRigidElementForce,
        rProcessInfo);
}

} // namespace Kratos