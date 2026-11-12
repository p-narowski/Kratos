// Surface-based DPD--DEM constitutive law for suspended DEM particles.

#include "custom_constitutive/dem_dpd_drag_ghost.h"

#include "custom_elements/dpd_particle.h"
#include "custom_elements/suspended_particle.h"

// Project includes
#include "DEM_application_variables.h"
#include "custom_elements/spheric_particle.h"

#include <algorithm>
#include <cmath>

namespace Kratos
{

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_DPD_DRAG_GHOST::Clone() const
    {
        return DEMDiscontinuumConstitutiveLaw::Pointer(
            new DEM_DPD_DRAG_GHOST(*this));
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw>
    DEM_DPD_DRAG_GHOST::CloneUnique()
    {
        return Kratos::make_unique<DEM_DPD_DRAG_GHOST>(*this);
    }

    std::string DEM_DPD_DRAG_GHOST::GetTypeOfLaw()
    {
        return "DEM_DPD_DRAG_GHOST";
    }

    void DEM_DPD_DRAG_GHOST::Check(Properties::Pointer pProp) const
    {
        /*
         * The base check requests the standard DEM properties. Keep it because
         * this class inherits the DEM contact-law infrastructure and the solver
         * may expect those properties to exist even though mKn/mKt are not used
         * in the DPD pair force below.
         */
        DEM_D_Linear_viscous_Coulomb::Check(pProp);

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_SMOOTHING_LENGTH))
            << "DEM_DPD_DRAG_GHOST requires DPD_SMOOTHING_LENGTH in properties "
            << pProp->Id() << "." << std::endl;

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_CUTOFF_RADIUS))
            << "DEM_DPD_DRAG_GHOST requires DPD_CUTOFF_RADIUS in properties "
            << pProp->Id() << "." << std::endl;

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_CONSERVATIVE_COEFF))
            << "DEM_DPD_DRAG_GHOST requires DPD_CONSERVATIVE_COEFF in properties "
            << pProp->Id() << "." << std::endl;

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_DISSIPATIVE_COEFF_NORMAL))
            << "DEM_DPD_DRAG_GHOST requires "
            << "DPD_DISSIPATIVE_COEFF_NORMAL in properties "
            << pProp->Id() << "." << std::endl;

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_DISSIPATIVE_COEFF_TANGENTIAL))
            << "DEM_DPD_DRAG_GHOST requires "
            << "DPD_DISSIPATIVE_COEFF_TANGENTIAL in properties "
            << pProp->Id() << "." << std::endl;

        const double h = (*pProp)[DPD_SMOOTHING_LENGTH];
        const double rc = (*pProp)[DPD_CUTOFF_RADIUS];

        KRATOS_ERROR_IF(h <= 0.0)
            << "DPD_SMOOTHING_LENGTH must be positive. Got " << h
            << " in properties " << pProp->Id() << "." << std::endl;

        KRATOS_ERROR_IF(rc <= 0.0)
            << "DPD_CUTOFF_RADIUS must be positive. Got " << rc
            << " in properties " << pProp->Id() << "." << std::endl;

        KRATOS_ERROR_IF(rc > 2.0 * h)
            << "DEM_DPD_DRAG_GHOST uses a compact cubic-spline-like kernel with "
            << "support rc <= 2*h. Got rc = " << rc << ", h = " << h
            << " in properties " << pProp->Id() << "." << std::endl;
    }

    void DEM_DPD_DRAG_GHOST::InitializeContact(
        SphericParticle *const element1,
        SphericParticle *const element2,
        const double indentation)
    {
        /*
         * This is not a spring-overlap DEM contact. The DPD--DEM interaction is
         * entirely calculated from surface distance, so disable inherited
         * tangential/normal spring stiffnesses.
         */
        mKn = 0.0;
        mKt = 0.0;
    }

    void DEM_DPD_DRAG_GHOST::InitializeContactWithFEM(
        SphericParticle *const element,
        Condition *const wall,
        const double indentation,
        const double ini_delta)
    {
        /*
         * DEM_DPD_DRAG_GHOST is intended only for DPD--suspended-particle pairs.
         * A DPD--wall law should be kept separate.
         */
        mKn = 0.0;
        mKt = 0.0;
    }

    double DEM_DPD_DRAG_GHOST::KernelWeight(
        const double r_surface,
        const double h,
        const double rc) const
    {
        if (r_surface >= rc)
        {
            return 0.0;
        }

        /*
         * r_surface < 0 means the DPD-particle center lies inside the DEM sphere.
         * Use q = 0, the maximum finite kernel weight, rather than passing a
         * negative coordinate to the spline.
         */
        const double r_effective = std::max(0.0, r_surface);
        const double q = r_effective / h;

        if (q < 1.0)
        {
            return 1.0 - 1.5 * q * q + 0.75 * q * q * q;
        }

        if (q < 2.0)
        {
            const double a = 2.0 - q;
            return 0.25 * a * a * a;
        }

        return 0.0;
    }

    double DEM_DPD_DRAG_GHOST::DissipativeWeight(
        const double r_surface,
        const double h,
        const double rc) const
    {
        const double w = KernelWeight(r_surface, h, rc);
        return w * w;
    }

    double DEM_DPD_DRAG_GHOST::GetDemRadius(
        SphericParticle *const element1,
        SphericParticle *const element2) const
    {
        return std::max(
            element1->GetRadius(),
            element2->GetRadius());
    }

    double DEM_DPD_DRAG_GHOST::CalculateCenterDistance(
        SphericParticle *const element1,
        SphericParticle *const element2) const
    {
        const auto &r_x1 = element1->GetGeometry()[0].Coordinates();
        const auto &r_x2 = element2->GetGeometry()[0].Coordinates();

        const double dx = r_x1[0] - r_x2[0];
        const double dy = r_x1[1] - r_x2[1];
        const double dz = r_x1[2] - r_x2[2];

        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    double DEM_DPD_DRAG_GHOST::CalculateNormalForce(
        SphericParticle *const element1,
        SphericParticle *const element2,
        const double indentation,
        double LocalCoordSystem[3][3])
    {
        return 0.0;
    }

    double DEM_DPD_DRAG_GHOST::CalculateNormalForce(
        SphericParticle *const element,
        Condition *const wall,
        const double indentation)
    {
        return 0.0;
    }

    double DEM_DPD_DRAG_GHOST::CalculateNormalForce(const double indentation)
    {
        return 0.0;
    }

    double DEM_DPD_DRAG_GHOST::CalculateCohesiveNormalForce(
        SphericParticle *const element1,
        SphericParticle *const element2,
        const double indentation)
    {
        return 0.0;
    }

    double DEM_DPD_DRAG_GHOST::CalculateCohesiveNormalForceWithFEM(
        SphericParticle *const element,
        Condition *const wall,
        const double indentation)
    {
        return 0.0;
    }

    void DEM_DPD_DRAG_GHOST::CalculateForcesWithFEM(
        const ProcessInfo &rProcessInfo,
        const double OldLocalElasticContactForce[3],
        double LocalElasticContactForce[3],
        double LocalDeltDisp[3],
        double LocalRelVel[3],
        double indentation,
        double previous_indentation,
        double ViscoDampingLocalContactForce[3],
        double &rCohesiveForce,
        SphericParticle *const element,
        Condition *const wall,
        bool &rSliding)
    {
        /*
         * This law is intentionally not a DPD--wall law. Return zero force if it
         * is assigned to a wall pair by mistake, rather than silently using
         * center/surface logic that has no physically defined DEM sphere.
         */
        LocalElasticContactForce[0] = 0.0;
        LocalElasticContactForce[1] = 0.0;
        LocalElasticContactForce[2] = 0.0;

        ViscoDampingLocalContactForce[0] = 0.0;
        ViscoDampingLocalContactForce[1] = 0.0;
        ViscoDampingLocalContactForce[2] = 0.0;

        rCohesiveForce = 0.0;
        rSliding = false;
    }

    void DEM_DPD_DRAG_GHOST::CalculateForces(
        const ProcessInfo &rProcessInfo,
        const double OldLocalElasticContactForce[3],
        double LocalElasticContactForce[3],
        double LocalDeltDisp[3],
        double LocalRelVel[3],
        double indentation,
        double previous_indentation,
        double ViscoDampingLocalContactForce[3],
        double &rCohesiveForce,
        SphericParticle *element1,
        SphericParticle *element2,
        bool &rSliding,
        double LocalCoordSystem[3][3])
    {
        KRATOS_TRY

        for (int k = 0; k < 3; ++k)
        {
            LocalElasticContactForce[k] = 0.0;
            ViscoDampingLocalContactForce[k] = 0.0;
        }

        rCohesiveForce = 0.0;
        rSliding = false;

        auto *p_dem_1 = dynamic_cast<SuspendedParticle *>(element1);
        auto *p_dem_2 = dynamic_cast<SuspendedParticle *>(element2);
        auto *p_dpd_1 = dynamic_cast<DPDParticle *>(element1);
        auto *p_dpd_2 = dynamic_cast<DPDParticle *>(element2);

        SuspendedParticle *p_dem = nullptr;

        if (p_dem_1 != nullptr && p_dpd_2 != nullptr)
        {
            p_dem = p_dem_1;
        }
        else if (p_dem_2 != nullptr && p_dpd_1 != nullptr)
        {
            p_dem = p_dem_2;
        }
        else
        {
            KRATOS_ERROR
                << "DEM_DPD_DRAG_GHOST requires one SuspendedParticle "
                << "and one DPDParticle. Pair: "
                << element1->Id() << ", " << element2->Id()
                << std::endl;
        }

        Properties &r_contact_properties =
            element1->GetProperties().GetSubProperties(
                element2->GetProperties().Id());

        const double h =
            r_contact_properties[DPD_SMOOTHING_LENGTH];

        const double rc =
            r_contact_properties[DPD_CUTOFF_RADIUS];

        const double a =
            r_contact_properties[DPD_CONSERVATIVE_COEFF];

        const double gamma_n =
            r_contact_properties[DPD_DISSIPATIVE_COEFF_NORMAL];

        const double gamma_t =
            r_contact_properties[DPD_DISSIPATIVE_COEFF_TANGENTIAL];

        const double center_distance =
            CalculateCenterDistance(element1, element2);

        KRATOS_ERROR_IF(center_distance <= 1.0e-12)
            << "DEM_DPD_DRAG_GHOST: coincident centers for particles "
            << element1->Id() << " and " << element2->Id()
            << std::endl;

        const double surface_distance =
            center_distance - p_dem->GetRadius();

        if (surface_distance >= 0.5 * rc)
        {
            return;
        }

        // Outside the DEM: exact mirror separation.
        // Inside the DEM: finite maximum-weight recovery force.
        const double ghost_distance =
            2.0 * std::max(0.0, surface_distance);

        const double wc =
            KernelWeight(ghost_distance, h, rc);

        const double wd =
            DissipativeWeight(ghost_distance, h, rc);

        LocalElasticContactForce[2] = a * wc;

        ViscoDampingLocalContactForce[0] =
            -2.0 * gamma_t * wd * LocalRelVel[0];

        ViscoDampingLocalContactForce[1] =
            -2.0 * gamma_t * wd * LocalRelVel[1];

        ViscoDampingLocalContactForce[2] =
            -2.0 * gamma_n * wd * LocalRelVel[2];

        KRATOS_CATCH("")
    }

} // namespace Kratos