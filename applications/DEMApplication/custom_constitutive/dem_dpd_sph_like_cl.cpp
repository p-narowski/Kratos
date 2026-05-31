// DPD-like simplified DEM constitutive law for small particles

#include "custom_constitutive/dem_dpd_sph_like_cl.h"
#include "DEM_application_variables.h"
#include "custom_elements/spheric_particle.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace Kratos
{

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_DPD_SPH_LIKE::Clone() const
    {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_DPD_SPH_LIKE(*this));
        return p_clone;
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_DPD_SPH_LIKE::CloneUnique()
    {
        return Kratos::make_unique<DEM_DPD_SPH_LIKE>(*this);
    }

    std::string DEM_DPD_SPH_LIKE::GetTypeOfLaw()
    {
        return "DEM_DPD_SPH_LIKE";
    }

    void DEM_DPD_SPH_LIKE::Check(Properties::Pointer pProp) const
    {
        DEM_D_Linear_viscous_Coulomb::Check(pProp);

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_SMOOTHING_LENGTH))
            << "Missing DPD_SMOOTHING_LENGTH in properties " << pProp->Id() << std::endl;

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_CUTOFF_RADIUS))
            << "Missing DPD_CUTOFF_RADIUS in properties " << pProp->Id() << std::endl;

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_CONSERVATIVE_COEFF))
            << "Missing DPD_CONSERVATIVE_COEFF in properties " << pProp->Id() << std::endl;

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_DISSIPATIVE_COEFF_NORMAL))
            << "Missing DPD_DISSIPATIVE_COEFF_NORMAL in properties " << pProp->Id() << std::endl;

        KRATOS_ERROR_IF_NOT(pProp->Has(DPD_DISSIPATIVE_COEFF_TANGENTIAL))
            << "Missing DPD_DISSIPATIVE_COEFF_TANGENTIAL in properties " << pProp->Id() << std::endl;
    }

    void DEM_DPD_SPH_LIKE::InitializeContact(
        SphericParticle *const element1,
        SphericParticle *const element2,
        const double indentation)
    {
        mKn = 0.0;
        mKt = 0.0;
    }

    void DEM_DPD_SPH_LIKE::InitializeContactWithFEM(
        SphericParticle *const element,
        Condition *const wall,
        const double indentation,
        const double ini_delta)
    {
        mKn = 0.0;
        mKt = 0.0;
    }

    double DEM_DPD_SPH_LIKE::KernelWeight(const double r, const double h, const double rc) const
    {
        if (r >= rc)
            return 0.0;

        const double q = r / h;

        if (q >= 0.0 && q < 1.0)
        {
            return 1.0 - 1.5 * q * q + 0.75 * q * q * q;
        }
        else if (q >= 1.0 && q < 2.0)
        {
            const double a = 2.0 - q;
            return 0.25 * a * a * a;
        }

        return 0.0;
    }

    double DEM_DPD_SPH_LIKE::DissipativeWeight(const double r, const double h, const double rc) const
    {
        const double w = KernelWeight(r, h, rc);
        return w * w;
    }

    double DEM_DPD_SPH_LIKE::CalculateNormalForce(
        SphericParticle *const element1,
        SphericParticle *const element2,
        const double indentation,
        double LocalCoordSystem[3][3])
    {
        return 0.0;
    }

    double DEM_DPD_SPH_LIKE::CalculateNormalForce(
        SphericParticle *const element,
        Condition *const wall,
        const double indentation)
    {
        return 0.0;
    }

    double DEM_DPD_SPH_LIKE::CalculateNormalForce(const double indentation)
    {
        return 0.0;
    }

    double DEM_DPD_SPH_LIKE::CalculateCohesiveNormalForce(
        SphericParticle *const element1,
        SphericParticle *const element2,
        const double indentation)
    {
        return 0.0;
    }

    double DEM_DPD_SPH_LIKE::CalculateCohesiveNormalForceWithFEM(
        SphericParticle *const element,
        Condition *const wall,
        const double indentation)
    {
        return 0.0;
    }

    void DEM_DPD_SPH_LIKE::CalculateForcesWithFEM(
        const ProcessInfo &r_process_info,
        const double OldLocalElasticContactForce[3],
        double LocalElasticContactForce[3],
        double LocalDeltDisp[3],
        double LocalRelVel[3],
        double indentation,
        double previous_indentation,
        double ViscoDampingLocalContactForce[3],
        double &cohesive_force,
        SphericParticle *const element,
        Condition *const wall,
        bool &sliding)
    {
        std::cout << "[DPD-FEM-DBG] ENTERED CalculateForcesWith FEM" << std::endl;

        KRATOS_TRY

        LocalElasticContactForce[0] = 0.0;
        LocalElasticContactForce[1] = 0.0;
        LocalElasticContactForce[2] = 0.0;

        ViscoDampingLocalContactForce[0] = 0.0;
        ViscoDampingLocalContactForce[1] = 0.0;
        ViscoDampingLocalContactForce[2] = 0.0;

        cohesive_force = 0.0;
        sliding = false;

        Properties &properties_of_this_contact =
            element->GetProperties().GetSubProperties(wall->GetProperties().Id());

        const double h = properties_of_this_contact[DPD_SMOOTHING_LENGTH];
        const double rc = properties_of_this_contact[DPD_CUTOFF_RADIUS];
        const double a = properties_of_this_contact[DPD_CONSERVATIVE_COEFF];
        const double gamma_n = properties_of_this_contact[DPD_DISSIPATIVE_COEFF_NORMAL];
        const double gamma_t = properties_of_this_contact[DPD_DISSIPATIVE_COEFF_TANGENTIAL];

        // Mirror-ghost wall distance based on kernel width (no radius)
        const double R_interaction = element->GetInteractionRadius();
        const double dist_to_wall = R_interaction - indentation; // ≈ DistPToB
        const double d = std::max(dist_to_wall, 1.0e-12);
        const double r = 2.0 * d; // mirror-ghost

        if (r < 1.0e-15 || r >= rc)
        {
            return;
        }

        const double wc = KernelWeight(r, h, rc);
        const double wd = DissipativeWeight(r, h, rc);
        // Enforce ghost velocity = -particle velocity:
        // relative velocity for DPD is v_p - v_g = v_p - (-v_p) = 2 v_p
        const double vrel_t0 = 2.0 * LocalRelVel[0];
        const double vrel_t1 = 2.0 * LocalRelVel[1];
        const double vrel_n = 2.0 * LocalRelVel[2];

        LocalElasticContactForce[2] = a * wc;

        // Dissipative wall force: separate tangential / normal damping (particle vs ghost)
        ViscoDampingLocalContactForce[0] = -gamma_t * wd * vrel_t0;
        ViscoDampingLocalContactForce[1] = -gamma_t * wd * vrel_t1;
        ViscoDampingLocalContactForce[2] = -gamma_n * wd * vrel_n;

        static int fem_dbg_count = 0;
        if (fem_dbg_count < 20)
        {
            const double fel =
                std::sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] +
                          LocalElasticContactForce[1] * LocalElasticContactForce[1] +
                          LocalElasticContactForce[2] * LocalElasticContactForce[2]);

            const double fvis =
                std::sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0] +
                          ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1] +
                          ViscoDampingLocalContactForce[2] * ViscoDampingLocalContactForce[2]);

            std::cout << "[DPD-FEM-DBG] "
                      << "indent=" << indentation
                      << " r=" << r
                      << " rc=" << rc
                      << " wc=" << wc
                      << " wd=" << wd
                      << " |Fel|=" << fel
                      << " |Fvis|=" << fvis
                      << " vrel=(" << vrel_t0 << "," << vrel_t1 << "," << vrel_n << ")"
                      << std::endl;

            ++fem_dbg_count;
        }

        KRATOS_CATCH("")
    }

    void DEM_DPD_SPH_LIKE::CalculateForces(
        const ProcessInfo &r_process_info,
        const double OldLocalElasticContactForce[3],
        double LocalElasticContactForce[3],
        double LocalDeltDisp[3],
        double LocalRelVel[3],
        double indentation,
        double previous_indentation,
        double ViscoDampingLocalContactForce[3],
        double &cohesive_force,
        SphericParticle *element1,
        SphericParticle *element2,
        bool &sliding,
        double LocalCoordSystem[3][3])
    {
        std::cout << "[DPD-PP-DBG] ENTERED CalculateForces" << std::endl;
        KRATOS_TRY

        LocalElasticContactForce[0] = 0.0;
        LocalElasticContactForce[1] = 0.0;
        LocalElasticContactForce[2] = 0.0;

        ViscoDampingLocalContactForce[0] = 0.0;
        ViscoDampingLocalContactForce[1] = 0.0;
        ViscoDampingLocalContactForce[2] = 0.0;

        cohesive_force = 0.0;
        sliding = false;

        Properties &properties_of_this_contact =
            element1->GetProperties().GetSubProperties(element2->GetProperties().Id());

        const double h = properties_of_this_contact[DPD_SMOOTHING_LENGTH];
        const double rc = properties_of_this_contact[DPD_CUTOFF_RADIUS];
        const double a = properties_of_this_contact[DPD_CONSERVATIVE_COEFF];
        const double gamma_n = properties_of_this_contact[DPD_DISSIPATIVE_COEFF_NORMAL];
        const double gamma_t = properties_of_this_contact[DPD_DISSIPATIVE_COEFF_TANGENTIAL];

        const auto &x1 = element1->GetGeometry()[0].Coordinates();
        const auto &x2 = element2->GetGeometry()[0].Coordinates();

        const double dx = x1[0] - x2[0];
        const double dy = x1[1] - x2[1];
        const double dz = x1[2] - x2[2];
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (r < 1.0e-15 || r >= rc)
        {
            return;
        }

        const double wc = KernelWeight(r, h, rc);
        const double wd = DissipativeWeight(r, h, rc);

        // Conservative: purely normal
        LocalElasticContactForce[2] = a * wc;

        // Dissipative: split normal / tangential
        ViscoDampingLocalContactForce[0] = -gamma_t * wd * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = -gamma_t * wd * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = -gamma_n * wd * LocalRelVel[2];

        static int pp_dbg_count = 0;
        if (pp_dbg_count < 20)
        {
            const double fel =
                std::sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] +
                          LocalElasticContactForce[1] * LocalElasticContactForce[1] +
                          LocalElasticContactForce[2] * LocalElasticContactForce[2]);

            const double fvis =
                std::sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0] +
                          ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1] +
                          ViscoDampingLocalContactForce[2] * ViscoDampingLocalContactForce[2]);

            std::cout << "[DPD-PP-DBG] "
                      << "r=" << r
                      << " rc=" << rc
                      << " wc=" << wc
                      << " wd=" << wd
                      << " |Fel|=" << fel
                      << " |Fvis|=" << fvis
                      << " vrel=(" << LocalRelVel[0] << "," << LocalRelVel[1] << "," << LocalRelVel[2] << ")"
                      << std::endl;

            ++pp_dbg_count;
        }

        KRATOS_CATCH("")
    }
} // namespace Kratos